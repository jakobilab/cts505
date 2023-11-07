# Introduction to R: RNA-seq analysis - part 2


## Fix gene names - offline without biomaRt

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

## Convert raw counts to normalized counts

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

## Remove very lowly expressed reads

```R
# i.e. all rows with ZERO reads in 9 conditions (columns)
table(rowSums(x$counts==0)==9)

# let edgeR do the work for us
# keeps genes with about 10 read counts or more in a minimum number of samples, where the number of samples is chosen according to the minimum group sample size.
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

# visualize cleaned data

# get mean and median values for library sizes, multiply by 1M
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

# By default, the function keeps genes with about 10 read counts or more in
# a minimum number of samples, where the number of samples is chosen according 
# to the minimum group sample size. The actual filtering uses CPM values rather 
# than counts in order to avoid giving preference to samples with large library
# sizes.
 
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


## Normalizing expression distributions

```R

# claculate normalization factors for each library
# use trimmed mean of M-values (TMM) (Robinson and Oshlack 2010) 
x <- calcNormFactors(x, method = "TMM")

# display norm factors for each library
x$samples$norm.factors

# duplicate object
x2 <- x

# manually assign 1 as normalization factor
x2$samples$norm.factors <- 1

# reduce read counts by 5% in first library 
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)

# increase read counts in library x5
x2$counts[,2] <- x2$counts[,2]*5

# plot data to visualize library sizes
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
```

## Unsupervised clustering - principal component analysis

```R

# compute log-CPM values
lcpm <- cpm(x, log=TRUE)

# set up two sets
# set 1: sample group
# set 2: lanes

par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)


plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")

plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

# interactive HTML plot!
glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)

```

## Building the design matrix

```R
# inspect sample name, group, lane associations
data.frame(colnames(x),group,lane)

# we are interested in difference within group
# we want to adjust for differences between lanes

design <- model.matrix(~0+group+lane)
rownames(design) <- colnames(x)
colnames(design) <- gsub("group", "", colnames(design))

# check design


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