# Introduction to R base: fundamental statistical analyses

## Read input files

### Excel files

```R
# first we need to install a suitable package
install.packages("openxlsx", dependencies = TRUE)

# lets look at the function call to read an XLSX file

read.xlsx(
  xlsxFile,
  sheet,
  startRow = 1,
  colNames = TRUE,
  rowNames = FALSE,
  detectDates = FALSE,
  skipEmptyRows = TRUE,
  skipEmptyCols = TRUE,
  rows = NULL,
  cols = NULL,
  check.names = FALSE,
  sep.names = ".",
  namedRegion = NULL,
  na.strings = "NA",
  fillMergedCells = FALSE
)

# sample call
sample_data <- readWorkbook("sample.xlsx", sheet=1)
```

### CSV (comma-separated value) files

```R
# Read tabular data into R
read.table(file, header = FALSE, sep = "", dec = ".")

# Read "comma separated value" files (".csv")
read.csv(file, header = TRUE, sep = ",", dec = ".", ...)

# Or use read.csv2: variant used in countries that 
# use a comma as decimal point and a semicolon as field separator.
read.csv2(file, header = TRUE, sep = ";", dec = ",", ...)

# Read TAB delimited files
read.delim(file, header = TRUE, sep = "\t", dec = ".", ...)
read.delim2(file, header = TRUE, sep = "\t", dec = ",", ...)

data <- read.csv("input.csv") # read contents of input.csv
print(data)
```


## The iris dataset

Fisher (1936). "The use of multiple measurements in taxonomic problems". Annals of Eugenics. 7 (2): 179–188. doi: https://doi.org/10.1111%2Fj.1469-1809.1936.tb02137.x

Full story on Wikipedia:  https://en.wikipedia.org/wiki/Iris_flower_data_set

## T-test

### Unpaired example

```R


attach(iris) #Attaches data example dataset

head(iris) # print first rows of table
colnames(iris) # print first rows of table

hist(iris)  #how does the data look?

hist(Sepal.Width)  # try again?

shapiro.test(iris) # does not work

shapiro.test(Sepal.Width) # works

t.test(x = setosa$Petal.Length) # One-sample t-test, two-sided
t.test(x = setosa$Petal.Length, y = versicolor$Petal.Length) # Welch two sample t-test (default: unpaired)

# write into variable
value <- t.test(x = setosa$Petal.Length, y = versicolor$Petal.Length) 

set.seed(0) # important for reproducible random numbers

ClevelandSpending <- rnorm(50, mean = 250, sd = 75)
NYSpending <- rnorm(50, mean = 300, sd = 80)

# generate e.g. blood pressure data for paired example
sample_data1 <- c(rnorm(1000, mean = 145, sd = 9)) # generate test data
sample_data2 <- c(rnorm(1000, mean = 138, sd = 8)) # generate more test data

t.test(sample_data1, sample_data1, paired = TRUE) # test random data


```

### Paired example

```R

t.test(x = setosa$Petal.Length, y = versicolor$Petal.Length, paired = T) # Paired t-test, two-sided (default: unpaired)

# generate e.g. blood pressure data for paired example
sample_data1 <- c(rnorm(50, mean = 145, sd = 9)) # generate test data
sample_data2 <- c(rnorm(50, mean = 138, sd = 8)) # generate more test data

binary_data <- c(rep("high", 50), rep("low", 50)) # generate non-numerical test data

all_data <- c(sample_data1, sample_data2) # merge two lists

t.test(sample_data1, sample_data2, paired = TRUE) # test random data

t.test(all_data ~ binary_data, var.equal = TRUE) # test with formula and binary data

var.test(sample_data1, sample_data2) # only test variances
```

## ANOVA

### One-way ANOVA

```R

attach(iris) #Attaches data example dataset

lm.iris <- lm(Sepal.Width~Species, data=iris) #makes a 'linear model' object 

library(car) # lets try to load this package

# install.packages("car") # if not installed yet?

leveneTest(lm.iris) # test for equal variances

anova(lm.iris) # run anova

summary(lm.iris1) # generate summary table

# we know: one group is different, but which one (if more than 2 groups provided)

iris.aov<-aov(Sepal.Width~Species,data=iris)  # calculate the test statistic for ANOVA AND determine whether there is significant variation among the groups 

# bonus: works with tukey test

# TukeyHSD: Compute Tukey Honest Significant Differences
TukeyHSD(iris.aov)
```

### Two-way ANOVA

```R

# predict if sepal width differs among species and community

v <-c ("high","low") # generate category 

iris$community<-v # set new data

# always test assumptions before running test:
# adding community variable to formula
# testing for an interaction between Species and community.

lm.iris2<-lm(Sepal.Width~Species*community,data=iris) #makes a 'linear model' object 

leveneTest(lm.iris2) # test variances

# we will find community is not sig. different, so can leave it out next

iris2.aov<-aov(Sepal.Width~Species,data=iris)
TukeyHSD(iris2.aov)

# however, lets play this out with community actually being significant

iris3.aov<-aov(Sepal.Width~Species*community,data=iris)
TukeyHSD(iris3.aov)

```



## Chi-square

```R

attach(iris)

iris <- iris[iris$Species != "virginica",] # subselect only 2 species

iris$Species<- factor(iris$Species) # make a factor

head(iris)

# break the variable into 2 categories, below and above the median value.
iris$Petal.Width.Cat <- cut(iris$Petal.Width, breaks = quantile(iris$Petal.Width, probs = seq(0, 1, 0.5)), include.lowest = TRUE)
levels(iris$Petal.Width.Cat) <- c("below", "above")

head(iris)

# Drop Petal.Width column
iris <- iris[,!(names(iris) %in% "Petal.Width")]

head(iris)

# make table
tab <- table(iris$Petal.Width.Cat, iris$Species)
tab

# To perform the chi-square test we will assume the null hypothesis as below:
# 
# H0 : The Petal.Width.Cat has no affect on the Species
#
# Consequently, the alternate hypotheses will be defined as below:
#
# Ha : The Petal.Width.Cat has some affect on the Species
# 
# We will perform the test using chisq.test function in R.

chi <- chisq.test(tab)
chi

```


## Plotting


```R

# install required packages
install.packages(ggplot2)
install.packages(ggsignif)

# load packages
library(ggplot2)
library(ggsignif)

# reload iris data
data(iris)

# simple example with automatic wilcox test
ggplot(iris, aes(x = Species, y = Sepal.Length)) +
  geom_boxplot() + # using `ggsignif` to display comparison of interest
  geom_signif(
    comparisons = list(c("versicolor", "virginica")),
    map_signif_level = TRUE
  )

# advanced example

dat <- data.frame(
  Group = c("S1", "S1", "S2", "S2"),
  Sub = c("A", "B", "A", "B"),
  Value = c(3, 5, 7, 8)
)

ggplot(dat, aes(Group, Value)) +
  geom_bar(aes(fill = Sub), stat = "identity", position = "dodge", width = .5) +
  geom_signif(
    y_position = c(5.3, 8.3), xmin = c(0.8, 1.8), xmax = c(1.2, 2.2),
    annotation = c("**", "NS"), tip_length = 0
  ) +
  geom_signif(
    comparisons = list(c("S1", "S2")),
    y_position = 9.3, tip_length = 0, vjust = 0.2
  ) +
  scale_fill_manual(values = c("grey80", "grey20"))


```