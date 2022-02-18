**[Return to the Course Home Page](../index.html)**

# Introduction to R
**Professor Patrick Biggs**

[Purpose](#purpose)<br>
[Introduction](#introduction)<br>


## Purpose

To learn how to use a new software environment with the aim of understanding some principles of different data types and their analysis.


## Introduction

The biological world has got much more quantitative and analytical as the way biology is now performed has changed.  Data volumes have increased markedly, and so has our requirement and ability to analyse such data.  The first part of “Bootcamp” will introduce you to two tools that will be used throughout the course:  the statistical software package R and the UNIX command line.  This will be followed up with a look into data visualisation in the third week.

Week 1 focusses on the statistical software package, R. It will also introduce you to the main software environment for the course – RStudio.  We will look at accessing this environment, using it and learning some basic concepts about data handling through the use of R, focussing on data tables.  

<!---In this lab session, you will experience use of R within the RStudio environment as an introduction to the main analyses we will perform for Modules 2 and 3 in the next part of the course.--->


## Lecture Overview

As a brief summary of the lecture, we learnt that the statistical language R, its history, and that it has to be installed on a computer, but it is available for all the major platforms (Windows, OS X and Linux). We talked about its use in this course as an environment in which to explore two aspects of bioinformatics -- transcriptomics and metagenomics -- in Modules 2 and 3 respectively.

We discussed the extensive help and support that these around R and its use, all the packages that can be downloaded to help you. We then talked about the interface to R through the terminal and then with a GUI version, and finally with full-blown environments such as RStudio, through which we also had a quick overview. It is through the use of RStudio, that we will learn how to use R, and that is the focus of today's practical.


<table><tr><td>
<strong>Note for 2022:</strong> Due to Covid-19 and the uncertainties this can still cause us, the practicals for 203328 this year are being performed inside a new cloud-based RStudio environment called "RStudio Cloud". We will discuss the use of this shortly.
</td></tr></table>


We looked at getting data in and out of R, its question and answer model of operation, and then at the five data types within R. After looking at those, we went on to data formats, again a major focus of today's lab practical. We also touched on the installation of packages either through R itself for general package installation, as well as with Bioconductor for more specific bioinformatics applications.

So, now on to the practical...

First thing, you have an electronic copy of this practical -- which you will see via the Stream and this site -- please don't copy and paste in the practical today! Why? Typing these commands is part of the hard wiring required to make them stick! By typing commands it gives you another opportunity to think about the command before you execute it and typing is more likely to trigger a big red flag in dangerous situations! So, think twice, type once, we will discuss this more next week.


**This module is an introduction to the worlds of Linux and R, with a mastery test as an assessment. Please use this practical lab as a learning experience for a new set of tools, which in this case might go beyond this course, and into future study.**


## Conventions used for this RStudio practical

In what follows, we shall use the following typographical conventions for working with R:





So, for example,

```R
 > x <- c(1, 5, 7, 14) `Enter`
```
means \"at the R console prompt `>`  '**\>**', type the code as written to make a vector `x` made up of the values 1, 5, 7 and 14, then press the key marked Enter\".

Don\'t forget to press the **\[Enter\]** key: commands are not sent to the R console until this is done. There will not be **\[Enter\]** at the end of any of the lines of code in this practical, it is taken that this is to be typed at the end of the line _unless otherwise stated_. It is also assumed that the first `>` is the R console prompt.



```R
> x = 1 ### comment
> x
[1] 1
### comment
```

## Practical sections start here

New material goes here...

<br><br><br>

[GitHub Markdown cheat sheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)


This is how you put in an image file:
```
<img src="graphics/giphy.gif" width="300"/>
```
<img src="graphics/giphy.gif" width="300"/>


