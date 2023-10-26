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

# change look w/o changing theme
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point(size=3) +
  labs(title = "Sepal Length vs. Sepal Width of Iris Species", x = "Length", y = "Width", col = "") +
  geom_smooth(se=FALSE ) +
  scale_color_manual(values=c("orange","purple","black")) +
  theme_classic()

# simplistic B/W version for non-color journal plots
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point(aes(size=Sepal.Width, shape=Species))+
  labs(title = "Sepal Length vs. Sepal Width of Iris Species", x = "Length", y = "Width", col = "") +
  geom_smooth(se=FALSE, color="black") +
  theme_bw()

# change look w/o changing theme
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_point(size=3) +
  labs(title = "Sepal Length vs. Sepal Width of Iris Species", x = "Length", y = "Width", col = "") +
  geom_smooth(se=FALSE ) +
  theme_classic()

# but I wanted a box plot!
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) +
  geom_boxplot() +
  labs(title = "Sepal Length vs. Sepal Width of Iris Species", x = "Length", y = "Width", col = "") +
  theme_dark()

# actually, nevermind, I meant histogram
ggplot(iris) +
  geom_histogram(aes(x=Sepal.Length,fill = Species), bins = 60) +
  labs(title = "Sepal Length", x = "Length", y = "Width", col = "")

# Frequency plot
ggplot(iris) +
  geom_freqpoly(aes(x=Sepal.Length,col = Species), bins = 60) +
  labs(title = "Sepal Length", x = "Length", y = "Width", col = "")

# what about continous color scales?
# we actively map the colors to a continous variable:
ggplot(data=iris,mapping=aes(x=Petal.Length,y=Petal.Width))+
  geom_point(aes(color=Sepal.Width))+
  geom_smooth(method="lm")

# modifiy titles of axis and legend
ggplot(data=iris,mapping=aes(x=Petal.Length,y=Petal.Width))+
  geom_point(aes(color=Sepal.Width))+
  geom_smooth(method="lm")+
  scale_color_continuous(name="New Legend Title")+
  labs(title="This Is A Title",
       subtitle="This is a subtitle",
       x=" Petal Length", 
       y="Petal Width",
       caption="This is a little caption.") 

# add interestings facets to our data
ggplot(data=iris,mapping=aes(x=Petal.Length,y=Petal.Width))+
  geom_point(aes(color=Sepal.Width))+
  geom_smooth(method="lm")+
  scale_color_continuous(name="New Legend Title")+
  scale_x_continuous(breaks=1:8)+
  labs(title="This Is A Title",subtitle="This is a subtitle",x=" Petal Length", 
       y="Petal Width", caption="This is a little caption.")+
  facet_wrap(~Species)

# add text for data points 
ggplot(data=iris,mapping=aes(x=Petal.Length,y=Petal.Width))+
  geom_point(aes(color=Sepal.Width))+
  geom_text(aes(label=Sepal.Width,hjust=0),nudge_x=0.5,size=3) +
  geom_smooth(method="lm")+
  scale_color_continuous(name="New Legend Title")+
  scale_x_continuous(breaks=1:8)+
  labs(title="This Is A Title",subtitle="This is a subtitle",x=" Petal Length", 
       y="Petal Width", caption="This is a little caption.")+
  facet_wrap(~Species)

# add labels for data points 
ggplot(data=iris,mapping=aes(x=Petal.Length,y=Petal.Width))+
  geom_point(aes(color=Sepal.Width))+
  geom_label(aes(label=Sepal.Width,hjust=0),nudge_x=0.5,size=3) +
  geom_smooth(method="lm")+
  scale_color_continuous(name="New Legend Title")+
  scale_x_continuous(breaks=1:8)+
  labs(title="This Is A Title",subtitle="This is a subtitle",x=" Petal Length", 
       y="Petal Width", caption="This is a little caption.")+
  facet_wrap(~Species)

```


## Using ggrepel for text labels

```R

# install package
install.packages("ggrepel")

# the only change here is geom_label_repel() instead of geom_label()

ggplot(data=iris,mapping=aes(x=Petal.Length,y=Petal.Width))+
  geom_point(aes(color=Sepal.Width))+
  geom_label_repel(aes(label=Sepal.Width,hjust=0),nudge_x=0.5,size=3) +
  geom_smooth(method="lm")+
  scale_color_continuous(name="New Legend Title")+
  scale_x_continuous(breaks=1:8)+
  labs(title="This Is A Title",subtitle="This is a subtitle",x=" Petal Length", 
       y="Petal Width", caption="This is a little caption.")+
  facet_wrap(~Species)

```


## Saving a plot as file

```R

# this specifies PDF format
# play around with width and heith arguments to find a good size
# default unit is inch
pdf('plot.pdf', width=5, height=5)

ggplot(data=iris,mapping=aes(x=Petal.Length,y=Petal.Width))+
  geom_point(aes(color=Sepal.Width))+
  geom_label_repel(aes(label=Sepal.Width,hjust=0),nudge_x=0.5,size=3) +
  geom_smooth(method="lm")+
  scale_color_continuous(name="New Legend Title")+
  scale_x_continuous(breaks=1:8)+
  labs(title="This Is A Title",subtitle="This is a subtitle",x=" Petal Length", 
       y="Petal Width", caption="This is a little caption.")+
  facet_wrap(~Species)

# we have to disable the device after plotting
dev.off()


# this specifies PNG format
# play around with width and heith arguments to find a good size
# default unit is pixel for PNG
png('plot.pdf', width=500, height=500)

ggplot(data=iris,mapping=aes(x=Petal.Length,y=Petal.Width))+
  geom_point(aes(color=Sepal.Width))+
  geom_label_repel(aes(label=Sepal.Width,hjust=0),nudge_x=0.5,size=3) +
  geom_smooth(method="lm")+
  scale_color_continuous(name="New Legend Title")+
  scale_x_continuous(breaks=1:8)+
  labs(title="This Is A Title",subtitle="This is a subtitle",x=" Petal Length", 
       y="Petal Width", caption="This is a little caption.")+
  facet_wrap(~Species)

# we have to disable the device after plotting
dev.off()

```