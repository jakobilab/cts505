# Introduction to R: day-to-day data wrangling and biomart

## loading, subsetting, and merging data frames

```R

gene1 <- read.table("genes1.tsv", header = T, sep = "\t")
gene2 <- read.table("genes2.tsv", header = T, sep = "\t")


```




## install biomart

```R

# install bioconductor first
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.17")

# now install biomart
BiocManager::install("biomaRt")

# load biomart

# does not work
library(biomart)

# load correctly
library(biomaRt)



```


## work with biomart

```R

# set up annotation
annotation <- "rnorvegicus_gene_ensembl"

# set up mart object
ensembl <- useMart("ensembl", dataset = annotation)

# query data
result <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description"), values = names(go_background_blank), mart = ensembl)

```