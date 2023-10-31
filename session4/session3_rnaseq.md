# Introduction to R: RNA-seq analysis

## Halloween special!
``` 
CATCGTGAYTTYATCAAGAACATGATTACTGGCACCTCGCAGGCCGACTGCGCCATCCTCATCATTGCCG
CTGGCACTGGCGAGTTCGAGGCTGGTATCTCCAAGGATGGCCAGACTCGCGAGCACGCCCTCCTTGCCTA
CACCCTTGGTGTTAAGCAGCTCATCGTTGCTATTAACAAGATGGACACCACCAAGTGGTCCGAGTCTCGT
TTCCAGGAAATCATCAAGGAGACGTCCAACTTCATCAAGAAAGTTGGCTATAACCCCAAGACGGTTCCCT
TCGTTCCCATTTCCGGCTTCAACGGAGACAACATGTTGACCGCCTCCACCAACTGCCCCTGGTACAAGGG
CTGGGAGAAGGAGACCAAGGGCGGCAAGGCCACTGGTAAGACCCTCCTCGAGGCCCTCGACGCCATCGAG
CCCCCCAAGCGTCCCGTCGACAAGCCTCTCCGTCTCCCGCTTCAGGATGTGTACAAGATSGGNGGTATCG
```


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

##    ENTREZID  SYMBOL TXCHROM
## 1    497097    Xkr4    chr1
## 2 100503874 Gm19938    <NA>
## 3 100038431 Gm10568    <NA>
## 4     19888     Rp1    chr1
## 5     20671   Sox17    chr1
## 6     27395  Mrpl15    chr1

# only keep first occurence of each gene name; remove duplicates
genes <- genes[!duplicated(genes$ENTREZID),]

```

## convert raw counts to normalized counts

```R
# we use CPM here (counts per million)
# this is very helpful since we do NOT need gene lengths for CPM calculation

# for FPKM (Fragments Per Kilobase of transcript per Million mapped reads) or 
# RPKM (Reads Per Kilobase per Million mapped reads) we require the gene lengths

# RPKM / FPKM:

# RPKM: single end, FPKM: paired-end

# a) Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
# b) Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
# c) Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.

cpm <- cpm(x)
# log CPM with added prior counts (2/(library size))
lcpm <- cpm(x, log=TRUE)

# show table
summary(lcpm)


```

## remove very lowly expressed reads

```R
# i.e. all rows with ZERO reads in 9 conditions (columns)
table(rowSums(x$counts==0)==9)

# let edgeR do the work for us
# keeps genes with about 10 read counts or more in a minimum number of samples, where the number of samples is chosen according to the minimum group sample size.
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

# visualize cleaned data


L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm[,i])
lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")


```