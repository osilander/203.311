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
-   Characters written in `this programming style` are commands to be typed into the computer as they stand. In a difference to practical 1, this will be for both functions and filenames, *[and all these commands should be typed in for this afternoon's practical]*.
-   Characters written after something like `###` indicate comments within the code that will help you.  They do not run code. 
-   Characters written in without `>` at the start of the line **(e.g. lines starting with "\[1\]")** indicate responses back from the R console inside RStudio that will help you too. 
-   Words inserted within square brackets **\[Enter\]** indicate keys to be pressed.

So, for example,

```R
 > x <- c(1, 5, 7, 14) [Enter]
```
means \"at the R console prompt `>`, type the code as written to make a vector `x` made up of the values 1, 5, 7 and 14, then press the key marked Enter\".

Don\'t forget to press the **\[Enter\]** key: commands are not sent to the R console until this is done. There will not be **\[Enter\]** at the end of any of the lines of code in this practical, it is taken that this is to be typed at the end of the line _unless otherwise stated_. It is also assumed that the first `>` is the R console prompt.

# Need to tidy up this bit

***Logging into the iMacs in C5-10***

The machines we are using for the course -- Apple iMacs -- have been updated for 2021, and there is a new login method for them. Usernames and passwords will be supplied to you in the lab, but please remember to ignore (i.e. cancel) the dialogue box about the network when you log in. We will be working within web browsers, and Firefox and Chrome are installed on the machines, or Safari if you are an Apple user. We will then login to RStudio Cloud using your personalised account.

***Opening up RStudio Cloud***

***[Setting up an account on RStudio Cloud]{.underline}***

We have organised a student account for you at <https://rstudio.cloud/> to be valid until the end of October 2021 (i.e. for the whole of the Semester). The information for your personalised link has been sent to you separately via a Stream announcement.

***[Important Note]{.underline}***

**Please see the separate document -- Access\_to\_RStudio\_Cloud2021 -- that has been provided both as a Stream announcement and in the "Taiwhanga \| Labs" section on the Stream site for setting up your RStudio Cloud environment. You must do this before you can do any work in this course, so being familiar with working with these projects when links and projects get sent to you by lecturers is extremely important.**

***[Our work today]{.underline}***

We are going to be using RStudio today as a convenient introduction in using R, so we are going to focus mostly on the console today (bottom left pane from the lecture). However, we will be using -- or rather you will need to be aware of -- the whole environment for this morning, and also for both next week, and Modules 2 and 3.

# End of tidying up this bit


## Basic operations

### Assigning variables

Let's start off by making the assigning the variable `x` with a value of 1:

```R
> x = 1
> x
[1] 1
```

We can also use the \<- operator to do the same thing:

```R
> x <- 1
> x
[1] 1
```

I prefer the`<-` as to me it is more demonstrative of what is going on i.e. we are making assigning `x` the value of `1`. Using `=` suggests equivalence, which is not strictly true.


### Functions

We have seen this in the code conventions above. Let's repeat that, but this time call it **vec**:

```R
> vec <- c(1, 5, 7, 14)
> vec
[1] 1 5 7 14
```


### Comments

We have also seen this in the code conventions above. Say we want to make a comment on some addition:

```R
> 1 + 1 ### this is a comment***
> [1] 2
```

This is our first example of using the 3 colours for this practical, and for all subsequent work within R. In the console the comments will be the same colour, but if you wrote this in the script box within RStudio, you would see the below, which are screenshots from my session. This is due to something called syntax highlighting, which is not relevant for today (but useful nonetheless). Of course, the script does not give you an answer, that is what the console is for.

# add in images for R comments here

An alternative way to help you with these practicals would be to rewrite the code like this, and this is the convention I may also use:

```R
### this is a comment and it can be at the top of a code bloc
> 1 + 1 
> [1] 2
```

Which one gets used depends on what we are trying to find out, or the principle we are currently exploring.


### Introduction to packages

As a quick introduction to installing packages, we will go into this is more detail later, this is as easy as:

```R
> install.packages("XXXXXX")
```
where `"XXXXXX"` is the package name you want to install. This requires an active internet connection to retrieve the package.






## Practical sections start here

```R
x = 34 ###comment
x

###comment
```

New material goes here...

<br><br><br>

[GitHub Markdown cheat sheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)


This is how you put in an image file:
```
<img src="graphics/giphy.gif" width="300"/>
```
<img src="graphics/giphy.gif" width="300"/>


