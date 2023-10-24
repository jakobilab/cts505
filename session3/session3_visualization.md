# Introduction to R: data visualization

## loading, subsetting, and merging data frames


## install ggplot2

```R

# bioconductor is installed, now we cann install ggplot2
BiocManager::install("ggplot2")

# does not work
library(ggplot2)

```


## data visualization 101

```R
# let's look at the iris dataset once more 

summary(iris)

# out very first plot
ggplot(iris, aes(x = Petal.Length, y = Sepal.Length, col = Species)) +
  geom_point()

# now something more intricate
# add a new title for x and y axis, remove Species from legend
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point() +
  labs(title = "Sepal Length vs. Sepal Width of Iris Species", x = "Length", y = "Width", col = "") 

# add smoothing line
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point() +
  labs(title = "Sepal Length vs. Sepal Width of Iris Species", x = "Length", y = "Width", col = "") +
  geom_smooth()

# remove gray areas
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point() +
  labs(title = "Sepal Length vs. Sepal Width of Iris Species", x = "Length", y = "Width", col = "") +
  geom_smooth(se=FALSE )

# add theming
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point() +
  labs(title = "Sepal Length vs. Sepal Width of Iris Species", x = "Length", y = "Width", col = "") +
  geom_smooth(se=FALSE ) +
  theme_dark()

# but I wanted a box plot!
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_boxplot() +
  labs(title = "Sepal Length vs. Sepal Width of Iris Species", x = "Length", y = "Width", col = "") +
  theme_dark()

# actually, nevermind, I meant histogram
ggplot(iris) +
  geom_histogram(aes(x=Sepal.Length,fill = Species), bins = 60) +
  labs(title = "Sepal Length", x = "Length", y = "Width", col = "") +
  scale_fill_brewer(palette = "BrBG")

# Frequency plot
ggplot(iris) +
  geom_freqpoly(aes(x=Sepal.Length,col = Species), bins = 60) +
  labs(title = "Sepal Length", x = "Length", y = "Width", col = "") +
  scale_fill_brewer(palette = "BrBG")

```