# Introduction to R base: fundamental statistical analyses

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

```R

attach(iris) #Attaches data example dataset

lm.iris1 <- lm(Sepal.Width~Species, data=iris) #makes a 'linear model' object 

library(car) # lets try to load this package

install.packages("car") # if not installed yet?






```

## Hypothesis testing

