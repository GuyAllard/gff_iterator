"""
Utility functions for working with GFF/GTF (version 2) files
"""

from collections import OrderedDict

# the 9 g(t/f)f columns
gff_fields = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attributes"
]

# quick check that I didn't do something stupid
assert len(gff_fields) == 9, "OOPS! - wrong number of gff fields defined ({})".format(len(gff_fields))


def attribute_to_value(attribute):
    """
    Convert a raw attribute string into a native python type
    """
    if attribute[0] == '"':
        return attribute.strip('"')
    elif attribute.isdigit() or (attribute.startswith('-') and attribute[1:].isdigit()):
        return int(attribute)
    else:
        try:
            return float(attribute)
        except ValueError:
            return attribute


def value_to_attribute(value):
    """
    Convert an attribute value from native python type to gff string representation
    """
    # convert numbers to string
    if isinstance(value, int) or isinstance(value, float):
        return str(value)

    # wrap strings in quotes
    return '"{}"'.format(value)


def column_9_dict(col_9):
    """
    Create a dictionary of key-value pairs from the fields present in column 9 of a row from the input GFF/GTF file
    An ordered dict is used here so that the origin attribute field order is maintained

    Unquoted attribute fields are treated as numeric, quoted attribute fields as strings

    :param col_9: string containing the 9th column of a row from the input file
    :return:      ordered dictionary of key-value pairs
    """
    fields = [field.strip() for field in col_9.split(";")]
    attributes = OrderedDict(
        [(pair[0], pair[1]) for pair in [field.split(" ") for field in fields if len(field) > 0]]
    )

    # remove quotes from string values, convert numbers to numeric types
    for key in attributes:
        attributes[key] = attribute_to_value(attributes[key])

    return attributes


def get_fields(line):
    """
    Create a dictionary containing the data from a single line of a g(t/f)f file

    :param line: string containing a single line of a g(t/f)f file
    :return:    dictionary containing the parsed data
    """
    row = [x.strip() for x in line.split("\t")]
    assert len(row) == 9, "Error - row is not valid GFF\n{}".format("\t".join(row))

    row[-1] = column_9_dict(row[-1])
    record = OrderedDict(zip(gff_fields, row))

    # start and end are integers
    record["start"] = int(record["start"])
    record["end"] = int(record["end"])

    # score and frame may be numbers, but may also be the '.' character
    try:
        record["score"] = float(record["score"])
        record["frame"] = int(record["frame"])
    except ValueError:
        pass

    return record


def to_string(data):
    """
    Convert a gff data dictionary back to a line of text for writing to a g(f/t)f file

    :param data: dictionary of gff/gtf data (see get_fields)
    :return:     string for a gff/gtf row
    """
    if len(data["attributes"]) > 0:
        attribute_string = "; ".join(
            ["{} {}".format(key, value_to_attribute(data["attributes"][key])) for key in data["attributes"]]
        ) + ";"
    else:
        attribute_string = ""

    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
        data["seqname"],
        data["source"],
        data["feature"],
        data["start"],
        data["end"],
        data["score"],
        data["strand"],
        data["frame"],
        attribute_string
    )


class Feature(object):
    """
    Base class to represent gff features
    """

    def __init__(self, data):
        """
        Initialize using a data parsed from a row of the gff file
        """
        self._data = data

    def __str__(self):
        """
        Return a string representation suitable for writing back to gff
        """
        return to_string(self._data)

    def __repr__(self):
        """
        Summary information
        """
        return str((self.chromosome, self.strand, self.start, self.end, self.transcript_id))

    @property
    def chromosome(self):
        return self._data["seqname"]

    @property
    def source(self):
        return self._data["source"]

    @property
    def start(self):
        """
        Start coordinate (remember, 1-based!!)
        """
        return self._data["start"]

    @property
    def end(self):
        return self._data["end"]

    @property
    def score(self):
        return self._data["score"]

    @property
    def strand(self):
        return self._data["strand"]

    @property
    def frame(self):
        return self._data["frame"]

    @property
    def attributes(self):
        return self._data["attributes"]

    @property
    def gene_id(self):
        try:
            return self._data["attributes"]["gene_id"]
        except KeyError:
            return None

    @property
    def transcript_id(self):
        try:
            return self._data["attributes"]["transcript_id"]
        except KeyError:
            return None

    @property
    def exon_id(self):
        try:
            return self._data["attributes"]["exon_id"]
        except KeyError:
            return None

    @property
    def extents(self):
        """
        Start and end position (remember, 1-based closed interval!!)
        """
        return [self.start, self.end]


class Exon(Feature):
    """
    Feature to represent exons
    """
    pass


class CDS(Feature):
    """
    Feature to represent CDSs
    """
    pass


class Container(Feature):
    """
    Abstract Container class to hold sub-features
    """

    def __init__(self, data):
        super(Container, self).__init__(data)
        self._children = []

    def __str__(self):
        """
        String representation of container and children
        suitable for writing to gff
        """
        lines = [super(Container, self).__str__()]  # feature line
        lines += [str(x) for x in self._children]  # line for each child feature
        return "\n".join(lines)

    def can_add(self, feature):
        raise NotImplementedError

    def add_child(self, child):
        """
        Add a child feature
        """
        if not self.can_add(child):
            raise TypeError("cannot add child")

        self._children.append(child)

    @property
    def children(self):
        return self._children

    @property
    def exons(self):
        return [x for x in self._children if isinstance(x, Exon)]

    @property
    def cds(self):
        return [x for x in self._children if isinstance(x, CDS)]


class Transcript(Container):
    """
    Container to represent a transcript
    """
    def can_add(self, feature):
        """
        Transcripts can contain CDS and Exon features
        :param feature:
        :return:
        """
        return (isinstance(feature, CDS) or isinstance(feature, Exon)) and \
               feature.transcript_id == self.transcript_id


class Gene(Container):
    """
    Container to represent a gene
    """
    @property
    def transcripts(self):
        return [x for x in self._children if isinstance(x, Transcript)]

    def can_add(self, feature):
        """
        Genes can contain Transcript, CDS and Exon features
        :param feature:
        :return:
        """
        return (isinstance(feature, Transcript) or isinstance(feature, Exon) or isinstance(feature, CDS)) and \
               feature.gene_id == self.gene_id


def make_feature(data):
    """
    make a feature from a parsed gff line
    :param data: a dictionary of gff fields (see 'get_fields')
    :return: A Feature object
    """
    if data["feature"] == "gene":
        return Gene(data)

    if data["feature"] == "transcript":
        return Transcript(data)

    if data["feature"] == "exon":
        return Exon(data)

    if data["feature"] == "CDS":
        return CDS(data)


def gff_iterator(gff_handle):
    """
    Iterate through a gff file, yielding hierarchical features

    :param gff_handle: an open, readable handle to a gff file (or file-like)
    :yield: a Feature object for each parsed feature in the gff file
    """
    hierarchy = []
    for line in gff_handle:
        # skip comment lines
        if line.startswith('#'):
            continue

        # make a Feature object for this line
        feature = make_feature(get_fields(line))

        # figure out where this feature belongs
        parent = None
        while len(hierarchy) > 0:
            parent = hierarchy[-1]
            try:
                parent.add_child(feature)
                parent = None
                break
            except (TypeError, AttributeError):
                hierarchy = hierarchy[:-1]

        if parent is not None:
            yield parent

        hierarchy.append(feature)

    if len(hierarchy) > 0:
        yield hierarchy[0]
