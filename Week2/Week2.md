**[Return to the Course Home Page](../index.html)**

# Introduction to UNIX
**Professor Patrick Biggs**

[Purpose](#purpose)<br>
[Introduction](#introduction)<br>



## Purpose

To learn how to use a new software environment with the aim of understanding some principles of the use of command line tools



## Introduction

The biological world has got much more quantitative and analytical as the way biology is now performed has changed.  Data volumes have increased markedly, and so has our requirement and ability to analyse such data.  The first part of “Bootcamp” will introduce you to two tools that will be used throughout the course:  the statistical software package R and the UNIX command line.  This will be followed up with a look into data visualisation in the third week.

Week 2 focusses on an environment used by many computational scientists, those studying biology and many other disciplines.  This environment is the UNIX command line, and takes us back to the days before everything was graphical and ‘point and click’ based.  This is a very powerful and useful way to analyse biological data, and we will start learning the basics of how to work in this different environment.



## Lecture overview

As a brief summary of the lecture, we learnt that Unix has been around for 50 years or so, and in some ways it is still quite similar to those initial operating systems (OS).  It could be argued it has suffered from being very basic in its visual richness, and in that regard it has not advanced much.  Under the hood, how it works, and its usability is a completely different matter.  If you are interested -- and yes, you can do this is this lab -- have a look at the following URL that gives a quick view of the evolution of Microsoft Windows as a comparison:

https://www.theguardian.com/technology/2014/oct/02/from-windows-1-to-windows-10-29-years-of-windows-evolution

Windows has clearly come on a long way with a shorter history than Unix. You will also remember from the lecture that the current version of the Apple OS -- OS X -- like all recent versions is also based on Unix.  That is why we could very easily get the basics of using the command line using a terminal on our iMacs we use for this course at Manawatu.  There is no need to install any other way to perform command line applications.  Windows machines (like those at Albany) cannot easily run Linux/Unix without additional software.  


<table><tr><td>
<strong>Note for 2022:</strong> Due to Covid-19 and the uncertainties this can still cause us, the practicals for 203311 this year are being performed inside a new cloud-based RStudio environment called "RStudio Cloud". We will discuss the use of this shortly.
</td></tr></table>


This lab practical is on basic usage of Linux/Unix.  From now on we will call this Linux.  We will explore some commands and try to work with a set of small files. Your aim is to get familiar with this environment, so when we use them for the rest of the course, you have the knowledge to perform the tasks asked for, as well as others we shall learn. 

This module is an introduction to the worlds of Linux and R, with a mastery test as an assessment.  Please use this practical lab as a learning experience for a new set of tools, which in this case might go beyond this course, and into future study.

Windows and macOS are very powerful OS and thousands of people working on their development to make them the polished products they are.  However, the time it takes new features and programs to be incorporated takes a very long time because of the huge corporations they have become.  

This is why a large part of scientific research computing uses Unix/Linux as an environment within which to work.  We can access new programs, try them out, and use them in our research if we want to.  We have greater control over the parameters and options that come with the programs.  In other words, many of these programs do not suffer from the "black box" nature of Windows and macOS, i.e. we do not know how they fully work under the surface.  The number of new programs for genome analysis that you can download is growing all the time.  LetÕs first go through some of the programs that come with a standard installation of Linux so that you can work with files and directories.

First thing, don't copy and paste in the practical today!  Why?  Typing these commands is part of the hard wiring required to make them stick!  Also, "with great power comes great responsibility", and the way you are working today has the power to wipe a computer -- no trash recovery -- completely gone forever (except for that system-wide backup you have made right?)!  By typing commands, it gives you another opportunity to think about the command before you execute it and typing is more likely to trigger a big red flag in dangerous situations!  **_Think twice, type once._**

**This module is an introduction to the worlds of Linux and R, with a mastery test as an assessment. Please use this practical lab as a learning experience for a new set of tools, which in this case might go beyond this course, and into future study.**



## Conventions used for this Linux practical

In what follows, we shall use the following typographical conventions:

- Characters written in `this programming style`, e.g. `ls` are commands to be typed into the computer as they stand.
- Characters written in italics like _`this programming style`_, e.g. _`file.txt`_ indicate non-specific file or directory names.
- Words within double arrows such as \<\<this\>\> are examples.
- Words inserted within square brackets [Ctrl] indicate keys to be pressed.

So, for example,

```bash
$ ls <<anydirectory>> [Enter]	- {NOTE: do not do this, it is an example.}
```
  
means "at the UNIX prompt $, type `ls` followed by the name of some directory, then press the key marked Enter."

Don't forget to press the [Enter] key: commands are not sent to the computer until this is done.



## Accessing the resources needed

### Computing

#### Manawatu (iMacs)

The machines we are using for the course -- Apple iMacs -- have been updated for 2022, and there is a new login method for them. Usernames and passwords will be supplied to you in the lab, but please remember to ignore (i.e. cancel) the dialogue box about the network when you log in. We will be working within web browsers, and Firefox and Chrome are installed on the machines, or Safari if you are an Apple user. We will then login to RStudio Cloud using your personalised account.

#### Albany (PCs)

check with Olin tomorrow


### Use of RStudio Cloud for this Linux practical

To start off with log into your RStudio Cloud account, and go to the `week1and2` project in the "Bootcamp_2022" workspace we used last week. In the lower left-hand side -- the area you used the console for last week -- you will see another tab called "Terminal".  This is where we are going to be working today, so click on that, and you should see something like the below.   

# add in image of terminal

The line that says `rstudio-user@application-2573314-deployment-6685328-mjnjr:/cloud/project$` is the current line that can accept input for you to do things with.  There are three things of note here:

| Part of the command  |  The response |
| ------------- | ------------- |
| Our username in the environment:  | `studio-user@application-2571986-deployment-6682510-766qp:` |
| Our location in our environment:  | `/cloud/project`  |
| The command prompt:               | `$`  |
   
I have used a testing environment for developing these resources, so my command prompt is different when it is shown in the below.  

On the right-hand side we have "window #4", which again needs to be switched to be in the "Files" tab, as in the below.  Folders of use for today give us the following view:

# add in image of folder structure

This is the base folder for your work today, and all files for this practical have been loaded into the project for you.

At this point we will start using the University of Surrey Unix notes[^2] to go through some of the commands I used in the lecture.  I will run the practical in the same way as that one.



## Exploring our folder

Because of the RStudio Cloud setup we are using, we can start exploring the files we have (we will discuss moving around the system a little more shortly), so don't worry about the command yet, type:

```bash
$ cd Lab2UnixAndR/
```

To find out what is in this directory, type

```bash
$ ls
```
The `ls` command (lowercase 'L' and lowercase 'S') lists the contents of your current working directory, which in this case is:

```bash
rstudio-user@application-2571986-deployment-6682510 766qp:/cloud/project/Lab2UnixAndR$ ls
starting.txt
```

NB: Because this location name is so long, for the sake of clarity, I will abbreviate the prompt to the end of our name and location only, thus the below is exactly the same:

```bash
766qp:/cloud/project/Lab2UnixAndR$ ls
starting.txt
```

You see a file called _`starting.txt`_.  Is that all that is there?  To find out type

```bash
$ ls -a
```

```bash
766qp:/cloud/project/Lab2UnixAndR$ ls -a
.  ..  .hiddenfile.txt  starting.txt    
```

We can now see another file called _`.hiddenfile.txt`_.  What is going on here?  Well, ls does not, in fact, cause all the files in your current directory to be listed, but only those ones whose name does not begin with a dot (.) Files beginning with a dot (.) are known as hidden files and usually contain important program configuration information. They are hidden because you should not change them unless you are very familiar with UNIX!!



## Assessment

To reiterate, for this, the second lab of this course, there is no direct assessment today.  What is required however, is an understanding of the principles we have learnt today, as these will be required for two things:

1. For the mastery test which accounts for 15% of the course.  This will take place on Friday 18-Mar-2022.
2. For a part of this course in which you will use RStudio to explore examples of genomic data.

The mastery test will test both Linux and R, so half of the marks will come from each part.  We have not done much biology today, but as a hint have a think about how we can work with biological data that can come in a variety of formats (e.g. sequences files in say fasta format) and what we can potentially do with those given what we have learnt today.  There will be a worksheet with test questions for you to work through as preparation for the test.



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

[^1]: Concepts for this part of the practical came from a tutorial found at: <http://www.ee.surrey.ac.uk/Teaching/Unix/index.html> 


<!--- ## Practical sections start here

New material goes here...

<br><br><br>

[GitHub Markdown cheat sheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)


This is how you put in an image file:
```
<img src="graphics/giphy.gif" width="300"/>
```
<img src="graphics/giphy.gif" width="300"/> --->


