**[Return to the Course Home Page](../index.html)**

### **10-May-2022 - PLEASE DO NOT MAKE A COPY OF THIS FILE UNTIL THIS LABEL IS REMOVED**

# Metagenomic and metabarcoding visualisation

**Professor Patrick Biggs**

[Purpose](#purpose)<br>
[Introduction](#introduction)<br>
[Lecture Overview](#lecture-overview)<br>
[Conventions used for this RStudio practical](#conventions-used-for-this-rstudio-practical)<br>
[Accessing the resources needed](#accessing-the-resources-needed)<br>

[Portfolio analysis](#portfolio-analysis)<br>
[Assessment](#assessment)<br>
[Contact](#contact)<br>


## Purpose

To return to a more in-depth analysis of metabarcoding data primarily using the `phyloseq` package within `R`.


## Introduction

Microbial diversity has been revolutionised by the rapid advances in sequencing technology, given us new insights into the role of the microbial world in all environments on Earth. Sequencing of metagenomes (all the DNA in a given sample etc) or of specific markers (metabarcodes e.g. 16S rRNA amplicons) have been crucial in this regard. Knowing what is there, and what it is potentially doing are two of the main research areas of interest.  How we then analyse the data, visualise and interpret the results is thus of importance.

Week 9 focusses on building from the lecture material (see below) to gain familiarity with the `phyloseq` package within `R`.


## Lecture Overview

XXXXXXXXXX

So, now on to the practical...

First thing, you have an electronic copy of this practical -- which you will see via the Stream and this site -- **please don't copy and paste in the practical today!** Why? Typing these commands is part of the hard wiring required to make them stick! By typing commands it gives you another opportunity to think about the command before you execute it and typing is more likely to trigger a big red flag in dangerous situations! So, think twice, type once.


---

## Conventions used for this RStudio practical

As a reminder, in what follows, we shall use the following typographical conventions for working with `R`:

-   Characters written in `this programming style` are commands to be typed into the computer as they stand.
-   Characters written in _`this programming style`_ are objects and file names.
-   Characters written after something like `###` indicate comments within the code that will help you.  They do not run code. 
-   Characters written in without `>` at the start of the line **(e.g. lines starting with "\[1\]")** indicate responses back from the R console inside RStudio that will help you too. 
-   Words inserted within square brackets **\[Enter\]** indicate keys to be pressed.

So, for example,

```R
 > x <- c(1, 5, 7, 14) [Enter]
```
means \"at the `R` console prompt `>`, type the code as written to make a vector `x` made up of the values 1, 5, 7 and 14, then press the key marked Enter\".

Don\'t forget to press the **\[Enter\]** key: commands are not sent to the `R` console until this is done. There will not be **\[Enter\]** at the end of any of the lines of code in this practical, it is taken that this is to be typed at the end of the line _unless otherwise stated_. It is also assumed that the first `>` is the `R` console prompt.

---

## Accessing the resources needed

### Computing

#### General

We will be working within web browsers, and Firefox and Chrome are installed on the machines, or Safari if you are an Apple user. We will then login to RStudio Cloud using your personalised account.

#### Manawatu (iMacs)

The machines we are using for the course -- Apple iMacs -- have been updated for 2022, and there is a new login method for them. Usernames and passwords will be supplied to you in the lab, but please remember to ignore (i.e. cancel) the dialogue box about the network when you log in.

#### Albany (PCs)

The machines are Windows PC's so logging into them should be the same as any other Massey Windows PC.

#### Outside Massey 

Most likely this is your own machine, so this should not be an issue.


### Our work today

We will return to our RStudio cloud environment to analyse the data we started working on in Week 7, and then further that analysis by following some tutorial examples on the pre-loaded `phyloseq` datasets.


### What we are going to do today?

XXXXXXXXXX

The commands below have been checked and should work fine.


---

## Exercise 1: completion of the DADA2 tutorial

The DADA2 pipeline produced a sequence table and a taxonomy table which is appropriate for further analysis in phyloseq. We'll also include the small amount of metadata we have – the samples are named by the gender (G), mouse subject number (X) and the day post-weaning (Y) it was sampled (e.g. GXDY is the format of the sample names).  So, now we shall build up this metadata file, but we shall do it inside of `R` using source files available to us:

```R
### load our required packages
> library(phyloseq)
> packageVersion("phyloseq")
[1] ‘1.30.0’

### modifying our seqtab.nochim object to remove the mock community
> seqtab <- makeSequenceTable(mergers[names(mergers) != "Mock"])
> seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)

### starting to make a dataframe for the samples by getting their names
> samples.out <- rownames(seqtab.nochim)
> samples.out
 [1] "F3D0"   "F3D1"   "F3D141" "F3D142" "F3D143" "F3D144" "F3D145" "F3D146" "F3D147" "F3D148" "F3D149"
[12] "F3D150" "F3D2"   "F3D3"   "F3D5"   "F3D6"   "F3D7"   "F3D8"   "F3D9" 
```

We will now make our other columns of data to add to the metadata:

```R
### set up the subject and gender
> subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
> subject
 [1] "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3" "F3"
> gender <- substr(subject,1,1)
> gender
 [1] "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F" "F"
> subject <- substr(subject,2,999)
> subject
 [1] "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3"
```

Now to define the data post weaning:

```R
### days post weaning
> day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
> day
 [1]   0   1 141 142 143 144 145 146 147 148 149 150   2   3   5   6   7   8   9

### setting up
> samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
> head(samdf)
  Subject Gender Day
1       3      F   0
2       3      F   1
3       3      F 141
4       3      F 142
5       3      F 143
6       3      F 144

### define early and late categories dependent on days post weaning to add to our data
> samdf$When <- "Early"
> samdf$When[samdf$Day>100] <- "Late"
> rownames(samdf) <- samples.out
> head(samdf)
       Subject Gender Day  When
F3D0         3      F   0 Early
F3D1         3      F   1 Early
F3D141       3      F 141  Late
F3D142       3      F 142  Late
F3D143       3      F 143  Late
F3D144       3      F 144  Late
```

Now we construct our initial phyloseq object as a combination to bring all these things together:

```R### this is all one line of code
> ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
+                sample_data(samdf), 
+                tax_table(taxa))
> ps
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 218 taxa and 19 samples ]
sample_data() Sample Data:       [ 19 samples by 4 sample variables ]
tax_table()   Taxonomy Table:    [ 218 taxa by 6 taxonomic ranks ]
```

**If you do not see this structure and output, please let a demonstrator know as soon as possible.**

We can start to do many analyses at this point.  Again as a starting point, we will have a look at the counts, alpha diversity, and an initial ordination to see if there are any differences in our metadata, and whether it has an effect on the results.

We shall start off by looking at the alpha diversity:

```R
### plotting diversity by day using the Shannon and Simpson measures
> plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
```

You should see something like the following:



---


## Portfolio analysis

Once again, this is a two part analysis for the week 9 Portfolio analysis. 

XXXXXXXXXXX

---


## Assessment

To reiterate this week's lab is only assessed via the Portfolio analysis as described above.  **No code shown here that has not been at least briefly described in ther Week 9 lecture will not the subject of assessment in the Mastery test for this Module that is taking place between Thursday 19-May-2022 and Friday 20-May-2022 online.**

The mastery test will test the contents of weeks 7 to 9, more information will follow soon.



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




<!-- <br><br><br>

[GitHub Markdown cheat sheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)


This is how you put in an image file:
```
<img src="graphics/giphy.gif" width="300"/>
```
<img src="graphics/giphy.gif" width="300"/> -->


