**[Return to the Course Home Page](../index.html)**

# 26-Apr 2022 - DO NOT MAKE A COPY OF THIS FILE UNTIL THIS LABEL IS REMOVED

# Introduction to R

**Professor Patrick Biggs**

[Purpose](#purpose)<br>
[Introduction](#introduction)<br>
[Lecture Overview](#lecture-overview)<br>
[Conventions used for this RStudio practical](#conventions-used-for-this-rstudio-practical)<br>
[Accessing the resources needed](#accessing-the-resources-needed)<br>
<!---[Basic operations](#basic-operations)<br>
[R data types](#r-data-types)<br>
[R data formats](#r-data-formats)<br>
[Functions within R](#functions-within-r)<br>
[Libraries, also known as Packages](#libraries,-also-known-as-packages)<br>
[Bioconductor](#bioconductor)<br>
[A small complete example](#a-small-complete-example)<br> --->
[Conclusion of the Practical](#conclusion-of-the-practical)<br>
[Assessment](#assessment)<br>
[Contact](#contact)<br>


## Purpose

To learn how to use a new software environment with the aim of understanding some principles of different data types and their analysis.


## Introduction

The biological world has got much more quantitative and analytical as the way biology is now performed has changed.  Data volumes have increased markedly, and so has our requirement and ability to analyse such data.  The first part of “Bootcamp” will introduce you to two tools that will be used throughout the course:  the statistical software package R and the UNIX command line.  This will be followed up with a look into data visualisation in the third week.

Week 1 focusses on the statistical software package, R. It will also introduce you to the main software environment for the course -– RStudio.  We will look at accessing this environment, using it and learning some basic concepts about data handling through the use of R, focussing on data tables.  


## Lecture Overview

As a brief summary of the lecture, we learnt that the statistical language R, its history, and that it has to be installed on a computer, but it is available for all the major platforms (Windows, OS X and Linux). We talked about its use in this course as an environment in which to explore two aspects of bioinformatics -- transcriptomics and metagenomics -- in the later parts of this course.

We discussed the extensive help and support that these around R and its use, all the packages that can be downloaded to help you. We then talked about the interface to R through the terminal and then with a GUI version, and finally with full-blown environments such as RStudio, through which we also had a quick overview. It is through the use of RStudio, that we will learn how to use R, and that is the focus of today's practical.

We looked at getting data in and out of R, its question and answer model of operation, and then at the five data types within R. After looking at those, we went on to data formats, again a major focus of today's lab practical. We also touched on the installation of packages either through R itself for general package installation, as well as with Bioconductor for more specific bioinformatics applications.

So, now on to the practical...

First thing, you have an electronic copy of this practical -- which you will see via the Stream and this site -- please don't copy and paste in the practical today! Why? Typing these commands is part of the hard wiring required to make them stick! By typing commands it gives you another opportunity to think about the command before you execute it and typing is more likely to trigger a big red flag in dangerous situations! So, think twice, type once, we will discuss this more next week.


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



## Accessing the resources needed

### Computing

#### general

We will be working within web browsers, and Firefox and Chrome are installed on the machines, or Safari if you are an Apple user. We will then login to RStudio Cloud using your personalised account.

#### Manawatu (iMacs)

The machines we are using for the course -- Apple iMacs -- have been updated for 2022, and there is a new login method for them. Usernames and passwords will be supplied to you in the lab, but please remember to ignore (i.e. cancel) the dialogue box about the network when you log in.

#### Albany (PCs)

The machines are Windows PC's so logging into them should be the same as any other Massey Windows PC.

#### outside Massey 

Most likely this is your own machine, so this should not be an issue.

### Our work today

We are going to be using RStudio today as a convenient introduction in using R, so we are going to focus mostly on the console today (bottom left pane from the lecture). However, we will be using -- or rather you will need to be aware of -- the whole environment for this morning, and also for both next week, and the rest of the course.











## Assessment

To reiterate, for this, the first lab of this course, there is no direct assessment today.  What is required however, is an understanding of the principles we have learnt today, as these will be required for two things:

1. For the mastery test which accounts for 15% of the course.  This will take place between Thursday 24-Mar-2022 and Friday 25-Mar-2022 online.
2. For a very large part of this course in which you will use RStudio to explore examples of genomic data.

The mastery test will test both Linux and R, so half of the marks will come from each part.  We have not done much biology today, but as a hint have a think about how we can work with biological data that can come in a variety of formats and what we can potentially do with those given what we have learnt today.  There will be a worksheet with test questions for you to work through as preparation for the test.



## Contact

I have two offices on the Manawatu campus (as I work for both SNS and SoVS), so I am not always in my Science Tower office D5-30. If you want to discuss anything, it's best to email me beforehand.

Prof Patrick Biggs,

Molecular Biosciences Group,

School of Natural Sciences

```
-. .-.   .-. .-.   .-. .-.   .
||\|||\ /|||\|||\ /|||\|||\ /|
|/ \|||\|||/ \|||\|||/ \|||\||
~   `-~ `-`   `-~ `-`   `-~ `-
```

<p.biggs@massey.ac.nz>
