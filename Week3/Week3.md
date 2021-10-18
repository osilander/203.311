**[Return to the Course Home Page](../index.html)**

# Introduction to Visualization
**Professor Murray Cox**

[Purpose](#purpose)<br>
[Introduction](#introduction)<br>
[A Basic Scatter Plot](#a-basic-scatter-plot)<br>
[Other Types of Graphs](#other-types-of-graphs)<br>
[Using Summary Statistics](#using-summary-statistics)<br>
[Choosing A Plot Type](#choosing-a-plot-type)<br>
[Critically Evaluating Your Data](#critically-evaluating-your-data)<br>
[P Hacking](#p-hacking)<br>

## Purpose

To learn how to design and build clear, simple and informative graphs.


## Introduction

Conveying quantitative information in graphical form sounds easy.  Perhaps surprisingly, it isn't.  A large number of studies dating from the 1960s show that *how* you present data really matters.  The way people interpret graphs can vary greatly based on superficial choices like axis ranges, plot type, color schemes and other seemingly unimportant features.  In this practical, we will explore some of the basics of good graph design.  We will make a series of graphs, and in the process explore some examples of poor pracice.


## A Basic Scatter Plot

Let's start by building a simple scatter plot for a set of 100 random data points drawn from a normal distribution.

```R
x.axis.data <- rnorm(100)
y.axis.data <- rnorm(100)
```

Making a basic scatter plot in R is straightforward.

```R
plot(x.axis.data, y.axis.data)
```

The data you generated is random, but running this command in R should produce a graph that looks something like this.

<img src="graphics/scatterplot_example1.jpg" width="200"/>

This graph is a good start, but there is a lot more we can do with it.  For instance, what exactly is plotted on those axes, and what units are they in?  We can add axis labels to give this information.

```R
plot(x.axis.data, y.axis.data, xlab="Variable 1 (Units)", ylab="Variable 2 (Units)")
```

You can also add a title if you want to.  Graphs used in papers and reports usually don't have titles (they often have a written caption instead), but titles can be useful when you are exploring a dataset and making lots of plots.

```R
plot(x.axis.data, y.axis.data, main="My Title", xlab="Variable 1 (Units)", ylab="Variable 2 (Units)")
```

You can do many other fancy things as well, such as changing the point character using the *pch* option...

```R
plot(x.axis.data, y.axis.data, main="My Title", xlab="Variable 1 (Units)", ylab="Variable 2 (Units)", pch=4)
```

...or changing the color of the points...

```R
plot(x.axis.data, y.axis.data, main="My Title", xlab="Variable 1 (Units)", ylab="Variable 2 (Units)", pch=4, col="blue")
```

... or removing the box around the plot.

```R
plot(x.axis.data, y.axis.data, main="My Title", xlab="Variable 1 (Units)", ylab="Variable 2 (Units)", pch=4, col="blue", bty="n")
```

<img src="graphics/scatterplot_example2.jpg" width="200"/>

As you can see, by making some relatively small changes, you can change how the graph looks in some fairly striking ways.

A good way to get a feel for what designs look good (and how to make them) is by searching for example plots online.  The [R Graph Gallery](https://www.r-graph-gallery.com) is a particularly nice site with lots of worked examples.

> **EXERCISE 1**
> Here is code to generate another dataset.
>
> x <- runif(100)
> y <- 2 + 3 \* x^2 + rnorm(100, 0, 0.25)
> 
> Let's assume this is production of a metabolite in mg/L (x axis data) relative to the cellular expression of a gene in read counts per million (y axis data).
> Can you make a scatter plot that clearly and simply shows the trend between the metabolite and gene expression?


## Other Types of Graphs

Alright, so that's all well and good if you want to make a scatter plot.  But what if you want some other type of plot, perhaps a bar chart?

Well, the R command for a bar chart is also fairly straightforward.  

```R
items <- c("A", "B", "C", "D", "E")
values <- sample(seq(0,100), 5)
barplot(values, names=items)
```

> **EXERCISE 2**
> Can you add *x* and *y* axis labels and a title to the bar chart?

And what if you want a pie chart?  That's also easy.

```R
items <- c("A", "B", "C", "D", "E")
values <- sample(seq(0,100), 5)
pie(values, labels=items)
```

There are a huge number of ways to plot data, including many plot types you have probably never heard of before.  Regardless of what sort of data you have, there will be multiple ways to plot it.  Often, the plot commands in R are very simple – at least for basic styles.

> **EXERCISE 3**
> Take another look through the [R Graph Gallery](https://www.r-graph-gallery.com).
> Choose a plot style that interests you and use the commands given online to replicate the plot.  For this exercise, stick to plots described as 'base R' rather than 'ggplot2'.


## Using Summary Statistics

Plotting data is time consuming, so before we do that, it is best practice to calculate summary statistics, such means (averages) and standard deviations.

Let's consider a dataset showing the gene expression of 142 genes together with their associated protein levels as determined by mass spectrometry.  There are four sets of data generated under four different environmental conditions.  

The first thing we need to do is load this dataset.

```R
load(url("https://github.com/mpcox/203.311/raw/main/Week3/files/expression.Rdata"))
```

Gene expression data should be plotted on the *x* axis and is therefore stored in variables called 'set1.x', 'set2.x', etc.  Protein level data should be plotted on the *y* axis and is stored in variables called 'set1.y', 'set2.y', etc.

Means (averages) are easy to calculate in R.

```R
mean(set1.x)
mean(set2.x)
mean(set3.x)
mean(set4.x)
```

> **EXERCISE 4**
> Calculate the mean (average) values of both the gene expression (*x*) and protein (*y*) levels for the four datasets.  Do they appear to differ under the four environmental conditions?

Standard deviations are also easy to calculate in R.

```R
sd(set1.x)
```

> **EXERCISE 5**
> Calculate the standard deviations of both the gene expression (*x*) and protein (*y*) levels for the four datasets.  Do they appear to differ under the four environmental conditions?

Finally, it can be helpful to calculate correlations between different variables.  A correlation analysis will tell you whether genes with high expression also have high levels of the corresponding protein.

```R
cor.test(set1.x, set1.y)
```

> **EXERCISE 6**
> Calculate correlations of gene expression (*x*) against protein (*y*) levels for the four datasets.  Do they appear to differ under the four environmental conditions?

At this point, you have probably identified that the four datasets are very similar, both in terms of their average expression and protein levels, regardless of which environmental conditions they were generated under. Interestingly, the correlation values are also small – there is little evidence that genes with high expression also produce high levels of the corresponding protein.

Just to confirm this, it's always a good idea to plot your data, just to confirm that there are no real differences between the datasets.

> **EXERCISE 7**
> Make scatter plots for the four datasets.


## Choosing A Plot Type

It's often a good idea to look at your data visually.  The challenge is: what plot type should you use?  After a while you begin to learn what sorts of plots are particularly good at representing certain types of data. Nevertheless, often you just have to try different plot types and see what works.

Let's consider genes grouped into functional categories, such as 'immune genes', 'enzymes', 'cell wall genes' and the like.  Classifying genes in this way is very common, and can be highly effective for distinguishing how certain types of genes change under different conditions (say, in cancer cells vs normal tissue).

Here is a dataset of five immune cell types and their percentage frequency under three conditions: normal tissue, the primary tumor, and a secondary tumor ('metastasis').

```R
cell.type <- c("eosinophils", "mast cells", "lymphocytes", "basophils", "neutrophils")
normal <- c(17.5, 21.5, 20.0, 17.0, 24.0)
primary <- c(20.0, 21.0, 19.0, 20.0, 20.0)
secondary <- c(21.0, 18.5, 20.0, 23.5, 17.0)
```  

> **EXERCISE 8**
> We learned earlier how to make pie charts.  Make pie charts for the normal tissue, primary tumor and secondary tumor.  Are there any clear differences in the proportions of the five immune cell types?

Because it's hard to know what plot type produces the clearest visualation of your data, it is often helpful to try other options and see if you like them better.  For group data, like that given here, bar charts are often used.

It can also sometimes be convenient to re-order your dataset – this can make it easier for your readers to understand the results.  Here, let's move the immune cell types around so they're in alphabetical order.

```R
cell.type <- c("basophils", "eosinophils", "lymphocytes", "mast cells", "neutrophils")
normal <- c(17.0, 17.5, 20.0, 21.5, 24.0)
primary <- c(20.0, 20.0, 19.0, 21.0, 20.0)
secondary <- c(23.5, 21.0, 20.0, 18.5, 17.0)
```

> **EXERCISE 9**
> We learned earlier how to make bar charts.  Make bar charts for the normal tissue, primary tumor and secondary tumor.  Do you prefer them to the pie charts?


## Looking Critically At Your Data

In 1980, Robert Jackman wrote an influential [paper](https://www.jstor.org/stable/2095134) looking at the link between social factors and income for a global range of countries.  Although this study is over 40 years old and is not based on genetic data, the analysis is widely known as a case study of how to look critically at data.  It is worth considering for that reason alone.

Let's start by loading the dataset.

```R
load(url("https://github.com/mpcox/203.311/raw/main/Week3/files/income.Rdata"))
```

This dataset contains two variables for 18 countries: the percentage turnout at national elections ('turnout') and the average level of inequality in the incomes of people ('income.inequality').

> **EXERCISE 10**
> We learned earlier how to calculate correlations.  Take this dataset and calculate the correlation between voter turnout and income inequality.

Hopefully you found, as did many many researchers before Jackman, that there is a strong negative correlation between voter turnout and income inequality (r = –0.78, p = 0.00013).  Because the P value is low, we would say that this is a statistically significant result.  Specifically, countries with a high voter turnout have less income inequality between people.

> **EXERCISE 11**
> Now do what Jackman did: make a scatter plot of voter turnout (*x* axis) vs income inequality (*y* axis).  What do you see?

In this instance, one country – South Africa – stands out as being very different to all the other countries.  In statistics, this process of looking for unusual data points is called [anomaly detection](https://en.wikipedia.org/wiki/Anomaly_detection) and South Africa would be said to be an *outlier*.  Although it is important to think critically before removing data from any analysis, this is a good example of where it would be reasonable to do that and then see if the analysis changes.

So what does the correlation between voter turnout and income inequality look like if we exclude South Africa?

First, we have to take out the South Africa data point.  Conveniently, it is the first entry in each variable.  Because each variable is just a vector of numbers, we can simply ask R to ignore the first country in the list and only consider countries 2 to 18.

```R
cor.test(turnout[2:18], income.inequality[2:18])
```

> **EXERCISE 12**
> Calculate the correlation between voter turnout and income inequality.  How has the correlation value changed?  What is the statistical significance (i.e., what is the probability value)?


## P Hacking








