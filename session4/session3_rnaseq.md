# Introduction to R: RNA-seq analysis

## this workflow is based on https://f1000research.com/articles/5-1408

## install software

```R

# bioconductor is installed, now we cann install ggplot2
BiocManager::install("RNAseq123")

# load libraries
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

```

## getting the data

```R

# URL to GEO data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63310

# direct URL to the file stored in GEO
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"

# downloading file
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 

# unpacking file
utils::untar("GSE63310_RAW.tar", exdir = ".")

# creating a list of gzipped txt files
# these contain the gene expression data
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")

# loop over the files to unpack from .gz to .txt
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

# check contents of first file
read.delim(files[1], nrow=5)

```

## load data into edgeR

```R
# read directly into edgeR
x <- readDGE(files, columns=c(1,3))

# check variable class
class(x)

# check dimensions
dim(x)

```

## fix sample names

```R
# only use the first 12 characters, remove file extensions
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))

# update column names
colnames(x) <- samplenames

# create groups - but as factor
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))

# assign groups
x$samples$group <- group

# assign lane information (optional, not always available)
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))

# assign lane to object
x$samples$lane <- lane

# check samples
x$samples

##                              files group lib.size norm.factors lane
## 10_6_5_11 GSM1545535_10_6_5_11.txt    LP 32863052            1 L004
## 9_6_5_11   GSM1545536_9_6_5_11.txt    ML 35335491            1 L004
## purep53     GSM1545538_purep53.txt Basal 57160817            1 L004
## JMS8-2       GSM1545539_JMS8-2.txt Basal 51368625            1 L006
## JMS8-3       GSM1545540_JMS8-3.txt    ML 75795034            1 L006
## JMS8-4       GSM1545541_JMS8-4.txt    LP 60517657            1 L006
## JMS8-5       GSM1545542_JMS8-5.txt Basal 55086324            1 L006
## JMS9-P7c   GSM1545544_JMS9-P7c.txt    ML 21311068            1 L008
## JMS9-P8c   GSM1545545_JMS9-P8c.txt    LP 19958838            1 L008

```


## fix gene names - offline without biomaRt

```R
# get IDs from row names
geneid <- rownames(x)

# use Mus.musculus package to offline conert IDs
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)

```