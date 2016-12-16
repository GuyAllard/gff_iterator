GFF Iterator
==============
A very simple python iterator for parsing GFF/GTF files


**Installation**  
To install with pip, run
```
pip install https://github.com/guyallard/gff_iterator/archive/v1.0.0.zip
```


**Usage**  
```python
from gff_iterator import gff_iterator
with open("somefile.gff", "r") as gff_file:
    for record in gff_iterator(gff_file):
        print(record)
```
