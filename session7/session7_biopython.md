# Introduction to Biopython

## Installing and loading the Biopython library

```bash
pip3 install biopython
```

```Python
# load library
import Bio

# check version of Biopython
Bio.__version__
```

## Working with sequence objects (Seq)

```Python
# load (sub) library
# specific loading keeps the namespace clean
from Bio.Seq import Seq


# generate a new minimal seq object
from Bio.Seq import Seq

# only the sequence (as string is required)
my_seq_obj = Seq("AGGTTTCCGTA")

# direct access
my_seq_obj

# access via print returns only string
print(my_seq_obj)

# implicitly convert to string
str(my_seq_obj)

# easily compute the complement
my_seq_obj.complement()

# or the reverse complement
my_seq_obj.reverse_complement()

# how long is the sequence?
print(len(my_seq_obj))

# cutting slices of sequences
my_seq_obj[4:10]

# lets make the above example more dynamic
my_seq_obj[4:len(my_seq_obj)-1]

# what is the GC fraction?
from Bio.SeqUtils import gc_fraction
gc_fraction(my_seq_obj)

# count simple motifs
my_seq_obj.count("TTT")

# concatenate sequences
my_other_seq_obj = Seq("ACGTAACCGG")
new_seq = my_seq_obj + my_other_seq_obj

# in silico transcription
new_seq.transcribe()

# in silico translation (* == stop codon)
new_seq.translate()

# seq objects are immutable
new_seq[5] = "C"

# create a mutable object
from Bio.Seq import MutableSeq
mutable_seq = MutableSeq(str(new_seq))
```

## Loading larger sequences from files

### FASTA file
```Python
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 

file = open("ls_orchid.fasta") 

records = parse(file, "fasta")

for record in records:    
   print("Id: %s" % record.id) 
   print("Name: %s" % record.name) 
   print("Description: %s" % record.description) 
   print("Annotations: %s" % record.annotations) 
   print("Sequence Data: %s" % record.seq) 
   print("Sequence Alphabet: %s" % record.alphabet)
```

### GBK file (GenBank), includes detailed annotations
```Python
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq 

file = open("ls_orchid.gbk") 

records = parse(file, "genbank")

for record in records:    
   print("Id: %s" % record.id) 
   print("Name: %s" % record.name) 
   print("Description: %s" % record.description) 
   print("Annotations: %s" % record.annotations) 
   print("Sequence Data: %s" % record.seq) 
   print("Sequence Alphabet: %s" % record.alphabet)
```


## Working with graphics output

### create a simple genomic diagram
```Python
# we have to install reportlab first!
from reportlab.lib import colors
from reportlab.lib.units import cm


from Bio.Graphics import GenomeDiagram
from Bio import SeqIO

# directly read GenBank file
record = SeqIO.read("NC_005816.gb", "genbank")

# get a first idea of the file contents
print(record)

# create the base diagram
gd_diagram = GenomeDiagram.Diagram("Yersinia pestis biovar Microtus plasmid pPCP1")
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()

# now we add features
for feature in record.features:
    if feature.type != "gene":
        # Exclude this feature
        continue
    if len(gd_feature_set) % 2 == 0:
        color = colors.blue
    else:
        color = colors.lightblue
    gd_feature_set.add_feature(feature, color=color, label=True)

# write out diagram
gd_diagram.draw(
    format="linear",
    orientation="landscape",
    pagesize="A4",
    fragments=4,
    start=0,
    end=len(record),
)
gd_diagram.write("plasmid_linear.pdf", "PDF")
gd_diagram.write("plasmid_linear.eps", "EPS")
gd_diagram.write("plasmid_linear.svg", "SVG")

# now lets make it circular
gd_diagram.draw(
    format="circular",
    circular=True,
    pagesize=(20 * cm, 20 * cm),
    start=0,
    end=len(record),
    circle_core=0.7,
)
gd_diagram.write("plasmid_circular.pdf", "PDF")

```


## Some recipes from the cookbook

### Simple FASTQ quality filter
```Python
from Bio import SeqIO
import gzip

count = 0

with gzip.open("SRR020192.fastq.gz", "rt") as handle:
    for rec in SeqIO.parse(handle, "fastq"):
        count += 1
print("%i reads" % count)

# minimum quality set to 20 below

good_reads = (
    rec
    for rec in SeqIO.parse("SRR020192.fastq", "rt")
        if min(rec.letter_annotations["phred_quality"]) >= 20
)

count = SeqIO.write(good_reads, "good_quality.fastq", "fastq")
print("Saved %i reads" % count)
```


### Read length histogram
```Python
from Bio import SeqIO
sizes = [len(rec) for rec in SeqIO.parse("ls_orchid.fasta", "fasta")]
len(sizes), min(sizes), max(sizes)

# Now that we have the lengths of all the genes (as a list of integers), we can use the matplotlib histogram function to display it.

from Bio import SeqIO

sizes = [len(rec) for rec in SeqIO.parse("ls_orchid.fasta", "fasta")]

import pylab

pylab.hist(sizes, bins=20)
pylab.title(
    "%i orchid sequences\nLengths %i to %i" % (len(sizes), min(sizes), max(sizes))
)
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.show()
```