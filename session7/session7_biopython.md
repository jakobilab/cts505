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