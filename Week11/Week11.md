**[Return to the Course Home Page](../index.html)**

# Transcriptomic Visualization
**Dr Olin Silander**

## Learning Objectives

1. Understand the steps involved in analysing RNA-seq data
2. Understand how small read count numbers can lead to uncertainty
3. Calculate differences in gene expression between samples (DGE)
4. Use four methods to visualise differences in gene expression: volcano plots, heatmaps, PCA, and UMAP 


## Introduction

Up to this point, we have covered Methods for visualising RNA-seq results via dimensional reduction - specifically, Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP). You have seen that in using these methods that certain "characteristics" (e.g. cocktail ingredients) of certain "things" (e.g. cocktails) can be reduced  


### Why are small numbers unreliable?

I have discussed in class that when only a small number of RNA-seq reads map to a gene, we cannot be as certain of its expression level than when large numbers of reads map to a gene. Let's first investigate why this is so.

Many of you have probably heard of different probability distributions, for example the [Normal distribution](https://en.wikipedia.org/wiki/Normal_distribution "Pi, why are you in this formula?") (or Gaussian distribution), the [exponential distribution](https://en.wikipedia.org/wiki/Exponential_distribution "It's got no memory!"), the [binomial distribution](https://en.wikipedia.org/wiki/Binomial_distribution "Will I succeed or not??"). Each of these is associated with different types of processes or samples, for example, [human height](https://tasks.illustrativemathematics.org/content-standards/HSS/ID/A/4/tasks/1020#:~:text=The%20heights%20of%20adult%20men,standard%20deviation%20of%202.5%20inches. "how tall are you?"), the [lifetime of a car battery](https://opentextbc.ca/introstatopenstax/chapter/the-exponential-distribution/ "I'm not talking about Tesla here"), or the [number of heads you'll get in ten coin flips](https://onlinestatbook.com/2/probability/binomial.html "Heads I win, tails you lose!"), respectively.

Closely related to the binomial distribution is The Poisson Distribution, which is the distribution one would expect in almost any case we are sampling a countable number of things. For example, after a very light rain, we could count the number of raindrops on different sidewalk squares. [These would be Poisson distributed](https://en.wikipedia.org/wiki/Poisson_scatter_theorem#:~:text=The%20expected%20number%20of%20raindrops,with%20intensity%20parameter%202%2F5. "But they would be hard to count"). Maybe we are interested in [how many Prussian cavalry](http://rstudio-pubs-static.s3.amazonaws.com/567089_c15d14f3d35b4edcbf13f33bbe775d4c.html "Not interested, thanks") are likely to be [killed by horse kicks in any given year](https://www.randomservices.org/random/data/HorseKicks.html "Where is Prussia anyway?"). Or maybe we're interested in the [number of calls we can expect at a call centre](https://www.statology.org/poisson-distribution-real-life-examples/ "Not answering my phone"). All of these are Poisson distributed.<br>

<img src="graphics/prob-dists.jpeg" width="800" title="C'mon this looks too complicated"/><br>
**There are lots of distributions and they're all related**<br><br>


The number of RNA-seq reads that map to a gene [is also Poisson distributed](https://www.biostars.org/p/84445/ "I knew it!") (largely speaking). As expected then, most packages for analysing RNA-seq data *model* the data as being Poisson distributed. Let us then see what this means for genes with high and low number of reads mapped to them.

<img src="graphics/poisson.jpeg" width="600" title="Early days"/><br>
**I literally do not understand this at all**. Credit: [xkcd](https://www.explainxkcd.com/wiki/index.php/12:_Poisson "But here's an explanation")<br><br>


Okay, let's make some pretend RNA-seq data. First, we will make a dataset with a very small number of read counts per gene. Maybe it was from a bad library prep. Or maybe from an under-represented barcode in what was otherwise a good sequencing run. To do this, we will sample our read counts as if they were Poisson. Navigate to your `R` console.

```R
# Let's keep this simple
# we'll use R's built-in Poisson random number generator
# In this case, the first argument is the number of random numbers
# and the second is the *mean* of the Poisson distribution
# We don't need to specify the variance, as the variance of 
# a Poisson distribution is EQUAL to the mean!
# And we'll get results for two "samples"
read.counts <- data.frame(
	sample_1=rpois(1000,3), 
	sample_2=rpois(1000,3),
	row.names=paste0("Gene_",1:1000))

# Did it work?
head(read.counts)
``` 



### Differential Gene Expression Analysis

As we will be using the DESeq2 package we will first need to install it. Navigate to your `R` console.

```R
# Get DESeq2 from the bioconductor website
library(BiocManager)
BiocManager::install("DESeq2")
```
We will use the `DESeq2` plotting tools, which require `ggplot2`. Apparently many of you are becoming quite familiar with this package for your portfolio assessments (well done!).

```R
# Install ggplot2 if necessary
# load the libraries for ggplot2 and DESeq2
# you know how to do this
```

The input data 
