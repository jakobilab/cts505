# Introduction to R: RNA-seq analysis - part 3


## Quickstart 

```R
library(limma)
# library(Glimma)
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
# x <- readDGE(files, columns=c(1,3))

count_data <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(count_data) <-  c("EntrezID", "GeneLength")

for (file in files){
  tmp <- read.delim(file, header=T, check.names=F)
  sample_name <- substring(file, 12, nchar(file))
  sample_name <- gsub(".txt","",sample_name)
  message(sample_name)
  colnames(tmp) <- c("EntrezID", "GeneLength", sample_name)
  count_data <- merge(count_data, tmp, by = c("EntrezID", "GeneLength"), all=T)
}
geneid <- count_data$EntrezID
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")
genes <- genes[!duplicated(genes$ENTREZID),]
count_data$Symbol <- genes$SYMBOL
rownames(count_data) <- count_data$EntrezID
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x <- DGEList(group=group, counts=count_data[,3:11], genes=count_data[,c("EntrezID", "Symbol", "GeneLength")])
x$samples$lane <- lane
x$samples$group <- group
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
x <- calcNormFactors(x, method = "TMM")
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

# or use a different contrast
lrt2 <- glmLRT(fit, contrast = contr.matrix[, 'LPvsML'])

# get the differentially expressd genes
topTags(lrt)

# generate a simple decision table
isDE <- as.logical(decideTestsDGE(lrt, p.value = 0.05))
DEnames <- rownames(x)[isDE]
summary(isDE)

edgerTable <- topTags(lrt, n = nrow(x))$table

```

## GO Enrichment

```R

# generate RPKM table
RPKM <- rpkm(x, x$genes$Length, prior.count = 0)

# generate list of gene names with rpkm > 1 in all samples
go_background_num <- nrow(RPKM[apply(RPKM[, -1], MARGIN = 1, function(x) all(x > 1)),])

# set this list to 0
go_background_blank <- as.vector(rep(0, go_background_num))

# assign gene gene IDs
names(go_background_blank) <- row.names(RPKM[apply(RPKM[, -1], MARGIN = 1, function(x) all(x > 1)),])
print(paste("background gene list contains: ", go_background_num, " genes", sep = ""))

tmp <- edgerTable[edgerTable$FDR <= 0.05,]

go_list_both <- go_background_blank
tmp$EntrezID <- as.character(tmp$EntrezID)

go_list_both[tmp$EntrezID] <- tmp$logFC


topDiffGenesBoth <- function(allScore) {
  return(allScore > 0 | allScore < 0)
}
test.stat <- new("elimCount", testStatistic = GOFisherTest, cutOff = 0.05)

allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Mm.eg.db", ID="Entrez")

go_bp <- new("topGOdata", 
             nodeSize = 10,
             ontology = "BP",
             geneSel = topDiffGenesBoth, 
             allGenes = go_list_both, 
             GO2genes=allGO2genes,
             annot = annFUN.GO2genes)

resultElimFisBP <- getSigGroups(go_bp, test.stat)
allRes_bp <- GenTable(go_bp, elimFisher = resultElimFisBP, topNodes = 200, numChar = 100)

```