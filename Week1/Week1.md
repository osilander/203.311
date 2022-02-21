**[Return to the Course Home Page](../index.html)**

# Introduction to R

**Professor Patrick Biggs**

[Purpose](#purpose)<br>
[Introduction](#introduction)<br>
[Lecture Overview](#lectureoverview)<br>


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


## Basic operations[^1]


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
[1] 2
```

This is our first example of using the 3 colours for this practical, and for all subsequent work within R. In the console the comments will be the same colour, but if you wrote this in the script box within RStudio, you would see the below, which are screenshots from my session. This is due to something called syntax highlighting, which is not relevant for today (but useful nonetheless). Of course, the script does not give you an answer, that is what the console is for.

# add in images for R comments here

An alternative way to help you with these practicals would be to rewrite the code like this, and this is the convention I may also use:

```R
### this is a comment and it can be at the top of a code block
> 1 + 1 
[1] 2
```

Which one gets used depends on what we are trying to find out, or the principle we are currently exploring.


### Introduction to packages

As a quick introduction to installing packages, we will go into this is more detail later, this is as easy as:

```R
> install.packages("XXXXXX")
```
where `"XXXXXX"` is the package name you want to install. This requires an active internet connection to retrieve the package.

>**Exercise A:**
>
>Install the package called `alluvial` into your current R session. Write the version column from your package window (lower right quadrant: window #4 for the lectures) in the box below (or in a notebook, text file etc.)
>
> <table><tr><td>
> ____________________
> </td></tr></table>

If we then want to use or load the package we type the following into the console:

```R
> library("XXXXXX")
```


### Getting help

As stated in the lecture, there are many ways to get help within R. Looking for answers on Google can also help -- your question is most likely to have been searched for before! R provides extensive documentation. For example, entering `?c` or `help(c)` at the prompt gives documentation of the function `c()` in R.

```R
> help(c)
```

Nothing returns to the console, but you get help in your help window (lower right quadrant: window #4 for the lectures). This is the difference that is found by using RStudio, but it is nothing to worry about.

If you are not sure about the name of the function you are looking for, you can perform a fuzzy search with the `apropos()` function.

```R
> ### your response might be different depending on the packages installed on your system
> apropos("nova")
> [1] "anova" "manova" "power.anova.test" "stat.anova" "summary.manova"
```

We can also use the `example` command to find out what the package can do.

>**Exercise B:**
> Load the package called `alluvial` into your current R session (there is a hint below the Exercise A answer box). Write the code in the box below to run as example from the `alluvial` package:
> <table><tr><td>
> ____________________
> </td></tr></table>

What it returns to the console and to your \#4 window as a plot is not relevant for now, and we can move on (looks pretty cool though!). Keep on pressing [Enter] (around 9 times) until you get an R command prompt (`>`) back on your console.


### Some basic mathematical operations

We can do basic mathematics on variables and numbers.

```R
> x + x ### examples -- make up your own here
> x * pi
> vec + 46
> 16 + 4567
> vec * vec
```

Try some ideas out for yourself, and place some code examples in the box below. Go back to the lecture notes (pages 10 and 11) to try some things out like addition, multiplication, subtraction etc. If it does not work, can you think why?

> <table><tr><td>
> ____________________
>
> </td></tr></table>
>
> What does `vec * vec` do to the values of `vec`? 

## R data types[^2]

There are 5 data types, and we can find out about them using the `class()` function. As is the lectures that is what we shall do here.


### Numeric

This is the default computational data type. If we assign a decimal value to a variable `a` as follows, `a` will be of the numeric type.

```R
> a <- 10.5   ### assign a decimal value* **
> a            ### print the value of a* **
[1] 10.5
> class(a)     ### print the class name of a
[1] "numeric"
```
Furthermore, even if we assign an integer to a variable `k`, it is still being saved as a numeric value.

```R
> k <- 1
> k           ### print the value of k* **
[1] 1
> class(k)    ### print the class name of k* **
[1] "numeric"
```

>**Exercise C:**
> You have variables `k` and `x` from the above. Check that both are of the numeric type and then write down what the numeric answer for `x + k - 1` is?
> <table><tr><td>
> ____________________
> </td></tr></table>


### Integer

To make integers in R, we need to specific them as such. We can be assured that y is indeed an integer by applying the **is.integer()** function.

```R
> y \<- as.integer(3)
> y               ### print the value of y
[1] 3 
> class(y)        ### print the class name of y
[1] "integer" 
> is.integer(y)   ### is y an integer?
> [1] TRUE
```

We can also coerce both numeric value and decimal values into an integer with the same **as.integer()** function.

```R
> as.integer(3.14) ### coerce a numeric value
[1]
> as.integer("5.27") ### coerce a decimal string
[1] 5
```

On the other hand, it is erroneous trying to parse a non-decimal string.

```R
> as.integer(\"Joe\")   ### coerce a non-decimal string
[1] NA
Warning message:
NAs introduced by coercion
```

### Complex

A complex value in R is defined via the pure imaginary value i.

```R
> z <- 1 + 2i     ### create a complex number
> z               ### print the value of z
[1] 1+2i
> class(z)        ### print the class name of z
[1] "complex"
```

The following gives an error as -1 is not a complex value.

```R
> sqrt(-1)   ### square root of -1
[1] NaN
Warning message:
In sqrt(-1) : NaNs produced
```

Instead, we have to use the complex value -1 + 0i, which we can do by coercing it.

```R
> sqrt(as.complex(-1))  #### square root of -1+0i
[1] 0+1i
```

### Logical

A logical value is often created via comparison between variables.

```R
> x <- 1; y <- 2    ### sample values, two on one line
> z <- x > y       ### is x larger than y?
> z                ### print the logical value
[1] FALSE
> class(z)         #### print the class name of z
[1] "logical"
```

### Character

One of the more important data types, which we can start to manipulate using various functions. A character object is used to represent string values in R. We convert objects into character values with the `as.character()` function:

```R
> x <- as.character(3.14)
> x                        ### print the character string
[1] "3.14"
> class(x)                 ### print the class name of x
> [1] "character"
```

Two character values can be concatenated with the paste function.

```R
> fname = "Joe"; lname = "Bloggs"
> paste(fname, lname)
[1] "Joe Bloggs"
```

However, it is often more convenient to create a readable string with the `sprint()` function, which has a C language syntax (hopefully you can see what is going on in this line of code).

```R
> sprintf("%s has %d dollars", "Dave", 1000)
[1] "Dave has 1000 dollars"
```

To extract a substring, we apply the **substr()** function. Here is an example showing how to extract the substring between the third and twelfth characters in a string.

```R
> substr("Mary has a little lamb.", start = 3, stop = 12)
[1] "ry has a l"
```

And to replace the first occurrence of the word "little" by another word "big" in the string, we apply the sub function.

```R
> sub("little", "big", "Mary has a little lamb.")
[1] "Mary has a big lamb."
```

More functions for string manipulation can be found in the R documentation.

```R
> help("sub")
```

>**Exercise D:**
> Using what you have learned in the section above, write the code to extract the first to seventh characters in the `Mary has a big lamb` string.
> <table><tr><td>
> ____________________
> </td></tr></table>
























### Contact

I have two offices on the campus (as I work for both SFS and SoVS), so I am not always in my Science Tower office D5-30. If you want to discuss anything, it's best to email me beforehand.

Prof Patrick Biggs,

Molecular Biosciences Group,

School of Natural Sciences

<p.biggs@massey.ac.nz>

[^1]: Concepts for this part of the practical came from a tutorial found at: <http://www.r-tutor.com/r-introduction>
[^2]: Concepts for this part of the practical came from a tutorial found at: <http://www.r-tutor.com/r-introduction/basic-data-types>
[^3]: Concepts for this part of the practical came from a tutorial found at: <http://www.r-tutor.com/r-introduction/>
[^4]: You might have to hit **\[Enter\]** in the console to see the plot -- this is a new change for 2019, due to the version of R.



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


