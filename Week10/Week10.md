**[Return to the Course Home Page](../index.html)**

# Single Cell Transcriptomics
**Dr Olin Silander**

## Learning Objectives

1. Understand the purpose of dimensional reduction techniques
2. Understand why dimensional reduction is useful for analysing large datasets
3. Explain the insights that dimensional reduction can give for RNA-seq data
4. Perform the steps necessary to implement dimensional reduction on a dataset
5. Interpret dimensional reduction plots
6. Disucss the advantages and disadvantages of two common dimensional reduction techniques, PCA and UMAP.


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

**I would argue that after dimensional reduction is the single most important technique you can apply for visulaisation of the data.**

Before reading further, please take three minutes and read [this quick intro](https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues "eigen-who?") before continuing.

Let us see how this works.

### The Meat and Potatoes

To gain some initial insight we will consider the food dataset from the UK. This is shown below.

<img src="graphics/meat-potatoes.png" width="400" title="That's a lotta potatoes N. Ireland"/><br>
**Yummy**

We will be using several new methods for data visalisation. Some of them are more intuitive than others. To help you gain some intuition, we will begin with a cocktail dataset and a tutorial derived from [here](https://juliasilge.com/blog/cocktail-recipes-umap/ "Cocktails how are they different") and [here](https://github.com/rfordatascience/tidytuesday/blob/master/data/2020/2020-05-26/readme.md "Cocktails lots of data").



<img src="graphics/edger-deseq2.jpeg" width="500"/><br>
**So many choices**<br>
source: https://kaopubear.top/blog/2022-03-20-donot-use-deseq2-edger-in-human-population-samples/

Next time: Why do bulk RNA-seq when we can do single-cell RNA-seq? What is single-cell RNA-seq? And more.

<img src="graphics/single-cell.jpeg" width="500"/><br>
**Problem solved.**<br><br>


<br><br><br>


