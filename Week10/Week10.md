**[Return to the Course Home Page](../index.html)**

# Single Cell Transcriptomics
**Dr Olin Silander**

## Learning Objectives

1. Understand the purpose of dimensional reduction techniques
2. Understand why dimensional reduction is useful for analysing large datasets
3. List three common methods for visualising RNA-seq data (volcano plot, heatmap, and dimensional reduction - PCA/UMAP/tSNE)
4. Explain the differences between each of the above visulaisation methods
5. Explain the insights that dimensional reduction can give for RNA-seq data
6. Perform the steps necessary to implement dimensional reduction on a dataset
7. Interpret dimensional reduction plots
8. Discuss the advantages and disadvantages of two common dimensional reduction techniques, PCA and UMAP.


## Introduction

### Dimensional Reduction

As we have learned throughout the Semester, a key aspect of data analysis is data visualisation. However, working in genomics, we often have extremely complicated data, and coming up with ways to visualise it in an intuitive yet objective manner is hard.

For example, last week you visualised differences in microbiome content based on hundreds of bacteria across tens of microbiome communites. How should we do this if we have tens of thousands of genes (i.e. variables) and tens of thousands of samples? We need to have an effective manner of *reducing* the number of variables so that we can visualise only one or two. How should we reduce the number of variables?

We will use **dimensional reduction** techniques. In this way we can objectively reduce tens of thousands of variables into *combinations of variables* so we can focus only on the one or two (or three) *most important (combinations) of variables* and determine which of our samples are most similar or different on the basis of these combinations.

Dimensional reduction is an important technique. In fact when you have any biological samples that have a large number of variables, e.g.
- microbiome samples with hundreds of different bacteria
- genotype data for individuals with hundreds of different SNPs
- gene expression data from cancer samples for hundreds of different genes
- phenotypic data from hundreds of dogs with information on tens of different phenotypes (e.g. height, weight, disposition, leg length, tail length, hair length, coat colour, and eye colour)

**I would argue that after dimensional reduction is the single most important technique you can apply for visulaisation of the data.** Here, we will focus on two main methods: Principal Component Analysis (PCA) and UMAP.

Before reading further, please take five minutes and read [this quick intro to PCA](https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues "eigen-who?") before continuing.

Becuase this is such an important concept, we are going to spend some time on this.
First some examples that have *nothing* to do with RNA or cells or sequencing. 
But hopefully they give us some insight into how dimensional reduction works and why it's important.

### The Meat and Potatoes

To gain some initial insight we will consider a food dataset from the UK. This is shown below. The dataset shows the consumption (in grams) per person per week of each of the foodstuffs.

<img src="graphics/meat-potatoes.png" width="400" title="That's a lotta potatoes N. Ireland"/><br>
**Yummy**<br><br>

Our aim here is to find out which of these countries :grimacing: differ the most in their diets. But of course diets are not one food or two foods, they are combinations of all foods. So which of these countries differ the most in the combination of all these foods?

We can already see that consumption of some types of foods differs more than others. For example, cereal consumption varies by about 5% between all countries. However, Welsh people drink more than 3.5 times as much alcohol than Irish people (*Northern Irish*).

We can also visualise this as a heatmap, which plots the same information, but more compactly. At the top of the heatmap is a dendrogram, which indicates how similar the countries are using [Ward's method](https://python-data-science.readthedocs.io/en/latest/unsupervised.html#agglomerative-clustering "it's a bit complicated"). N. Ireland appears the most different, while England and Wales appear the most similar. Note, importantly, that it does not say much about *how* similar.

<img src="graphics/diet_heat.png" width="300" title="cookin"/><br>
**It's getting hot in here**<br><br>

But we can also figure out which countries are the most similar or different in their combined diet. For this, we can perform a PCA. This finds the combinations of diet items (components) that vary the most between countries. We can then take these components and plot them. Below, I show the first two components (Dim1 and Dim2) - these are the two most important components. Clearly, Wales and N. Ireland differ the most in the combinations of items in their diets. I have made the x-axis (pc1) approximately three times longer than the y-axis (Dim2), as Dim1 accounts for approximately three times more variance (68%) than Dim2 (25%). 

<img src="graphics/diet_pca.png" width="700" title="N. Ireland is a different place"/><br>
**England is central to it all**<br><br>

Not only that, we can visualise which diet items *contribute* to those components. This is shown below.

<img src="graphics/diet_comps.png" width="400" title="aha it's the vegetables"/><br>
**What are "other veg", Wales?**<br><br>

Now we can see that Dimension (Component) 1 consists primarily of sugars and other_veg, all of which the Welsh consume more of, especially compared to N. Ireland. Dimension 2 consists primarily of the Irish tendency to eat a lot of potatoes (with some avoidance of alcohol). But most importantly, we have shrunk our 17-dimensional dataset to two dimensions that account for 68.3 + 24.9 = 93.2% (!) of the variance in the original 17 dimensions.

Okay, let's repeat this ourselves, with a new dataset.

### Hold my beer - Increasing Sample Size and Dimensions
We will move on to a cocktail dataset and a tutorial derived from [here](https://juliasilge.com/blog/cocktail-recipes-umap/ "Cocktails how are they different") and [here](https://github.com/rfordatascience/tidytuesday/blob/master/data/2020/2020-05-26/readme.md "Cocktails lots of data").


**At this point, open your `terminal`**.

Next, download the data from [here](data/all_cocktails.tab). If you have forgotten how to do that, ask your neighbour.

Navigate to your `RStudio` tab and read this file into `R`. Use the `read.table()` function to do this. Ensure that you use the `header=T` argument and assign it to a reasonably named variable (you can choose, but note that this is a dataset on cocktails. Or, for simplicity you can name it `cocktails_df` as that will match the code below).

We now have a dataset of cocktails and their ingredients. Take a look at the dataset, for example with `head` or `summary`.

Next we need to load a few libraries before we do our first analysis

```R
# it's the tidyverse!
install.packages("tidymodels")
library(tidymodels)

# it's for cats!
install.packages("forcats")
library(forcats)

# it's an obscure stats package!
install.packages("embed")
library(embed)

# it's a famous plotting package!
install.packages("ggplot2")
library(ggplot2)
```

Now we start on the path toward PCA.

```R
# This is a recipe
# We don't really sweat the details
# We just paste the code (*all* of it)
# But I put comments in if you're curious

# Tell the recipe what's happening but have no model ( ~. )
pca_rec <- recipe(~., data = cocktails_df) %>%
  # Note that the name and category of cocktail 
  # are just labels not predictors
  update_role(name, category, new_role = "id") %>%
  # Normalise means so all the 
  # variables have mean 0 and stdev 1
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors())
```

```R
# Actually do the PCA by "preparing"
# the "recipe"
pca_prep <- prep(pca_rec)
# and a smidge of tidying
tidied_pca <- tidy(pca_prep, 2)
```
Phew.

Now we can visualise the results. First, let's take a look at
the first two principal components. Remember, these are the
*combinations* of ingredients that contain the most variance
(in other words, what combinations of ingredients differ
the most between cocktail drinks).

Below we use the [ggplot](https://ggplot2.tidyverse.org/index.html "Thanks, Hadley!") plotting package.
This uses the idea of a *grammar* of graphics
and is among the most popular plotting methods in R

```R
# juice gets the results of the recipe
# and feeds it using %>% to the plotting function
juice(pca_prep) %>%
# the plotting, include the cocktail name
  ggplot(aes(PC1, PC2, label = name)) +
  # make the points colored by category
  geom_point(aes(color = category), alpha = 0.7, size = 2) +
  # add text 
  geom_text(check_overlap = TRUE, hjust = "inward", size = 2) +
  # and don't add more colour anywhere
  labs(color = NULL)
```

Wow, a few cocktails are quite different from others. What's in an Applejack punch?

```R
# we have a cocktail of interest
my.cocktail <- "Applejack Punch"
# Let's find the ingredients and assign it to a variable, "ingredients"
# You should be able to see what the code below is doing

ingredients <- cocktails_df[cocktails_df$name==my.cocktail,]
# Now we can see the ingredients
# What is this code doing? It has a new method, which()
# that we use to only report the ingredients that are 
# greater than 0 (i.e. they're in the cocktail)
cocktails_df[cocktails_df$name==my.cocktail,which(ingredients>0)]

# repeat the above steps but with a 
# cocktail of your choice, or, for example, this one:
my.cocktail <- "Sphinx Cocktail"
```

Not only can we say which cocktails are most different, we can see which are most similar.
This would help us suggest new but similar drinks to customers, for example
if you were bartending.

```R
# Here, we choose a couple of cocktails to look at
# You can choose these or different ontes
my.cocktails <- c("Silver Fizz", "Peach Blow Fizz")
# Another new method, the for loop
# we repeat the same as above, but
# "loop" over all values of the my.cocktails vector above
# of course here, coi means "cocktail of interest"
for(coi in my.cocktails) {
  ingredients <- cocktails_df[cocktails_df$name==coi,]
  print(cocktails_df[cocktails_df$name==coi,which(ingredients>0)])
}
```
So similar yet so different.

So what have we discovered? We have found that dimensional reduction is a powerful method to let us determine what variables (or, *combinations* of variables, e.g. diet items or cocktail ingredients) differentiate samples (e.g. countries or cocktails). We can use this to objectively determine which samples are the most similar, and which are the most different. We can also determine which (combinations of) variables are most responsible for *making* these samples different.

But enough of that, onwards and upwards.

### Who map? UMAP
A second commonly used method for dimensional reduction is UMAP (Uniform Manifold Approximation). UMAP is not as easy as PCA to understand from an algorithmic point of view. It is, however, an extremely powerful method for reducing dimensions while preserving the original structure of the data (i.e. the relative relationships and distances between samples). Please take a couple of minutes to browse [this site](https://pair-code.github.io/understanding-umap/ "Wooly Mammoth!!"). Scroll down to the Woolly Mammoth and adjust the parameters. Specifically, try n_neighbors = 100 and min_dist = 0.25.

Okay, let's go through this quickly just so we can compare to our previous results. We make almost exact the same recipe as before:

```R
umap_rec <- recipe(~., data = cocktails_df) %>%
  update_role(name, category, new_role = "id") %>%
  step_normalize(all_predictors()) %>%
  # this is the different step, where we use the step_umap function
  step_umap(all_predictors())

umap_prep <- prep(umap_rec)

umap_prep
```

And juice our results to plot it.

```R
juice(umap_prep) %>%
  ggplot(aes(UMAP1, UMAP2, label = name)) +
  geom_point(aes(color = category), alpha = 0.7, size = 2) +
  geom_text(check_overlap = TRUE, hjust = "inward", size = 2) +
  labs(color = NULL)
```

Woah. Compare this to the previous PCA result. What is different? Although *both of these methods have the same goal - dimensional reduction - you can see that there are very different results.* Here we can see that UMAP does not aim to find what variables differentiate samples the most (thereby *stretching* some dimensions considerably and *shrinking* others). Rather, UMAP aims to find ways to reduce dimensions while maintaining groupings. If we consider the Woolly Mammoth example from the link above, PCA would find that the variable with the most variation is (largely speaking) length and width. It would then project onto these, leaving differences between the left and right side nearly non-existent. You canimagine, for example, that the two tusks would thus become indistinguishable. However, this is not at all true for UMAP. It groups the tusks (as they are near) but keeps them separate. Similar for the left and right legs.

## RNA-seq

Now we can begin our RNA-seq journey. To do this, we will begin at the beginning, with some RNA-seq reads from human samples. These are from [here](data/data.tar "THE TAR FILE"). This is a relatively large file (20 Mb or so). Let's make a fresh directory for this analysis, perhaps `rnaseq`. Do that, change into that directory, and please download the RNA-seq reads now (`wget`). Note that much of this part of the lab is based off the tutorial [here](https://github.com/griffithlab/rnaseq_tutorial/wiki/RNAseq-Data "Awesome tutorial").

 Let's first untar the [tarball](https://en.wikipedia.org/wiki/Tar_(computing "Sticky!") so that we see the files inside.

 ```bash
# -x extracts -v is verbose -f is the file
tar -xvf data.tar

# upon successful untarring, remove the tarball!
rm data.tar
```

<img src="graphics/tar_2x.png" width="700"/><br>
**You're lucky I told you the command.**<br><br>


Next, check the that `.fastq` files look as we expect. Use your trusty friend, `seqkit`.

Let's do a quick QC step. Before, we used the comprehensive QC tool `fastp`. This is an excellent tool as an all-in-one QC and trimmer. Now we will use `fastqc` and look at a report generated by `multiqc` for rapid QC assessment

We can do a quick install:

```bash
# both at once
# might briefly redline your RAM
mamba install -c bioconda fastqc multiqc
``` 
Now we run the QC steps

```bash
# yes this might redline the RAM again
fastqc *fastq.gz

# it can't be this easy, can it?
mu1tiqc .

```

Go ahead and click on the multiqc report file (`.html`). (Open in your browser)For each of the `.fastq` files we can see a summary of its statistics. Note that there is aclickable menu on the left, and a toolbox available on the right (click the "toolbox" tab). The toolboax allows you to do things like colour smaples by group or hide specific samples. We will not worry about that. However, one important statistic we can see is that there a lot of sequence duplicates.

We are not going to worry about the adaptor trimming step of QC, as *I have already done this for you*. However, under normal circumstances this could be fatal for your pipeline.



<img src="graphics/edger-deseq2.jpeg" width="500"/><br>
**So many choices**<br>
source: https://kaopubear.top/blog/2022-03-20-donot-use-deseq2-edger-in-human-population-samples/

### Next Time

Next time: Why do bulk RNA-seq when we can do single-cell RNA-seq? What is single-cell RNA-seq? And more.

<img src="graphics/single-cell.jpeg" width="500"/><br>
**Problem solved.**<br><br>


<br><br><br>


