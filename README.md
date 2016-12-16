GFF Iterator
==============
A very simple python iterator for parsing GFF/GTF files


**Installation**  
from within this directory, run  
```
pip install .
```


**Usage**  
```python
from gff_iterator import gff_iterator
with open("somefile.gff", "r") as gff_file:
    for record in gff_iterator(gff_file):
        print(record)
```
