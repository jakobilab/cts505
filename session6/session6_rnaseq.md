# Introduction to R: RNA-seq analysis - part 3


## Quickstart 

```R
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)
x <- readDGE(files, columns=c(1,3))
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
geneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")
genes <- genes[!duplicated(genes$ENTREZID),]
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors
```

## Building the design matrix

```R

# Full guide to design matrices:
# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html

# inspect sample name, group, lane associations
data.frame(colnames(x),group,lane)
 
# colnames.x. group lane
# 1   10_6_5_11    LP L004
# 2    9_6_5_11    ML L004
# 3     purep53 Basal L004
# 4      JMS8-2 Basal L006
# 5      JMS8-3    ML L006
# 6      JMS8-4    LP L006
# 7      JMS8-5 Basal L006
# 8    JMS9-P7c    ML L008
# 9    JMS9-P8c    LP L008

# we are interested in difference within group
# we want to adjust for differences between lanes

design <- model.matrix(~0+group+lane)
rownames(design) <- colnames(x)
colnames(design) <- gsub("group", "", colnames(design))

# set up contrasts
contr.matrix <- makeContrasts(
   BasalvsLP = Basal-LP, 
   BasalvsML = Basal - ML, 
   LPvsML = LP - ML, 
   levels = colnames(design))
contr.matrix

```

## Differential gene expression

```R

#The square root of the common dispersion gives the
# coefficient of variation of biological variation
y <- estimateDisp(x, design, robust=TRUE)
y$common.dispersion

# fit the generalized linear model (glms)
fit <- glmFit(y, design)

# conduct likelihood ratio tests for LP vs. Basal using our contrast
 lrt <- glmLRT(fit, contrast = contr.matrix[, 'BasalvsLP'])

# get the differentially expressd genes
 topTags(lrt)
```

## GO Enrichment

```R

RPKM <- rpkm(x, x$genes$Length, prior.count = 0)

  RPKM <- RPKM[, c(grep(baseline, colnames(RPKM)), grep(focus, colnames(RPKM)))]
  labels_current <- labels[ c(grep(baseline, labels), grep(focus, labels))]

  # print(colnames(RPKM))
# get list of gene names with rpkm > 1 in all samples

go_background_num <- nrow(RPKM[apply(RPKM[, -1], MARGIN = 1, function(x) all(x > 1)),])
go_background_blank <- as.vector(rep(0, go_background_num))

names(go_background_blank) <- row.names(RPKM[apply(RPKM[, -1], MARGIN = 1, function(x) all(x > 1)),])

message(paste("background gene list contains: ", go_background_num, " genes", sep = ""))
```