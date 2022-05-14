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
- gene expression data from cancer samples for hundreds of differemt genes
- phenotypic data from hundreds of dogs with information on tens of different phenotypes (e.g. height, weight, disposition, leg length, tail length, hair length, coat colour, eye colour, etc.

**I would argue that after dimensional reduction is the single most important technique you can apply for visulaisation of the data.** Here, we will focus on two main methods: UMAP and Principal Component Analysis (PCA).

Before reading further, please take five minutes and read [this quick intro to PCA](https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues "eigen-who?") before continuing.

Becuase this is such an important concept, we are going to spend some time on this.
First some examples that have *nothing* to do with RNA or cells or sequencing. 
But hopefully they give us some insight into how dimensional reduction works and why it's important.

### The Meat and Potatoes

To gain some initial insight we will consider a food dataset from the UK. This is shown below.

<img src="graphics/meat-potatoes.png" width="400" title="That's a lotta potatoes N. Ireland"/><br>
**Yummy**<br><br>

Our aim here is to find out which of these countries :grimacing: differ the most in their diets. But of course diets are not one food or two foods, they are combinations of all foods. So which of these differ the most in the combination of all these foods?

We can already see that consumption of some types of foods differs more than others. For example, cereal consumption varies by about 5% between all countries. However, Welsh people drink more than 3.5 times as much alcohol than Irish people (*Northern Irish*).

We can also visualise this as a heatmap, which plots the same information, but more compactly.

<img src="graphics/diet_heat.png" width="300" title="cookin"/><br>
**It's getting hot in here**<br><br>

But we can also figure out which countries are the most similar or different in their combined diet. For this, we can perform a PCA. This finds the combinations of diet items (components) that vary the most between countries. We can then take these components and plot them. Below, I show the first two components (pc1 and pc2) - these are the two most important components

<img src="graphics/diet_pca.png" width="700" title="N. Ireland is a different place"/><br>
**Nearly matches the geography**<br><br>

Not only that, we can visualise which diet items *contribute* to those components. This is shown below.

<img src="graphics/diet_comps.png" width="400" title="aha it's the vegetables"/><br>
**What are "other veg", Wales?**<br><br>

Now we can see that Dimension (Component) 1 consists primarily of sugars, other_veg, and alcohol, all of which Scots, English, and Welsh consume more of. Dimension 2 consists primarily of the Irish tendency to eat a lot of potatoes. But most importantly, we have shrunk our 17-dimensional dataset to two dimensions that contain 68.3 + 24.9 = 93.2% (!) of the variance in the original 17 dimensions.

Okay, let's repeat this ourselves, with a new dataset.

### Hold my beer - Increasing Sample Size and Dimensions
We will move on to a cocktail dataset and a tutorial derived from [here](https://juliasilge.com/blog/cocktail-recipes-umap/ "Cocktails how are they different") and [here](https://github.com/rfordatascience/tidytuesday/blob/master/data/2020/2020-05-26/readme.md "Cocktails lots of data").


**At this point, open your `terminal`**.

Next, download the data from [here](data/all_cocktails.tab). If you have forgotten how to do that, ask your neighbour.

Navigate to your `RStudio` tab and read this file into `R`. Use the `read.table()` function to do this. Ensure that you use the `header=T` argument and assign it to a reasonably named variable (you can choose, but note that this is a dataset on cocktails. Or, for simplicity you can name it "cocktails_df" as that will match the code below).

We now have a dataset of cocktails and their ingredients. Take a look at the dataset, for example with `head` or `summary`.

We need to load a library before we do our first analysis

```R
# it's the tidyverse!
install.packages("tidymodels")
library(tidymodels)

install.packages("forcats")
library(forcats)
```

Now we do the PCA

```R
# This is a recipe
# We don't really sweat the details
# We just paste the code (*all* of it)
# But I put comments in if you're curious

# Tell the recipe what's happening with no model ( ~. )
pca_rec <- recipe(~., data = cocktails_df) %>%
  # Note that the name and category of cocktail 
  # are just labels not predictors
  update_role(name, category, new_role = "id") %>%
  # Normalise means so all the 
  # variables have mean 0 and stdev 1
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors())

# Actually do the magic PCA by "preparing"
# the "recipe"
pca_prep <- prep(pca_rec)
# and a smidge of tidying
tidied_pca <- tidy(pca_prep, 2)
```
Phew.

Now we can visualise the results. First, let's take a look at
the first two principal components. Remember, these are the
combinations of ingredients that contain the most variance
(in other words, what combinations of ingredients differ
the most between cocktail drinks).

```R
# Now we use the ggplot plotting package
# This uses the idea of a *grammar* of graphics
# and is among the most popular plotting methods in R

# juice gets the results of the recipe
# and feed it %>% to the plotting function
juice(pca_prep) %>%
# the plotting
  ggplot(aes(PC1, PC2, label = name)) +
  # make the points colored by category
  geom_point(aes(color = category), alpha = 0.7, size = 2) +
  # add text 
  geom_text(check_overlap = TRUE, hjust = "inward", family = "Helvetica", size = 8) +
  # but don't colour the text
  labs(color = NULL)
```

Wow, a few cocktails are quite different from others. What's in an Applejack punch?

```R
# we have a cocktail of interest
my.cocktail <- "Applejack Punch"
# lets find the ingredients and assign it to a variable, "ingredients"
ingredients <- cocktails_df[cocktails_df$name==my.cocktail,]
# Now we can see the ingredients
# What is this code doing? It has a new method, which()
cocktails_df[cocktails_df$name==my.cocktail,which(ingredients>0)]

# repeat the above steps but with a new cocktail
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
# "loop" over all values of the vector above
for(coi in my.cocktails) {
  ingredients <- cocktails_df[cocktails_df$name==coi,]
  print(cocktails_df[cocktails_df$name==coi,which(ingredients>0)])
}
```
So similar yet so different.

Enough of that.


### RNA-seq

<img src="graphics/edger-deseq2.jpeg" width="500"/><br>
**So many choices**<br>
source: https://kaopubear.top/blog/2022-03-20-donot-use-deseq2-edger-in-human-population-samples/

### Next Time

Next time: Why do bulk RNA-seq when we can do single-cell RNA-seq? What is single-cell RNA-seq? And more.

<img src="graphics/single-cell.jpeg" width="500"/><br>
**Problem solved.**<br><br>


<br><br><br>


