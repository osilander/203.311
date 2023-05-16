**[Return to the Course Home Page](../index.html)**

# Portfolio Assistance Week 4
**A/Prof Olin Silander**

[Purpose](#purpose)<br>

## Purpose

1. Understand possible ways to progress through a portfolio assessment.
2. Find basic ways to visualise data from the Week 4 portfolio assessment.
3. Understand resources for visualising different types of data. 

## Introduction

The Portfolio Assessment from Week 4 is the first part of this Assessment. For many of you it will be the first time you have tried to visualise data, or certainly to visualise data outside of Excel. The purpose here is to show possible ways to move step-by-step through the Assessment.

## Background

The Portfolio Assessment is written as follows:

1. The GC (guanine and cytosine) content of a genome (and read) is known to affect how easy it is to sequence. For example, regions of genomes that are GC-rich are often under-represented in sequencing data. One question that arises from this is whether GC content is correlated with quality scores for either Illumina or Oxford Nanopore sequencing reads. **Please address this question.**<br>

    You will have to generate new data on GC content and quality using the `seqkit fx2tab` subcommand (make sure you use the help command, so: `seqkit fx2tab --help`), and combine it with `R` plotting methods to answer this question. Please do this, and state whether you think there is, or is not, a link between GC content and read quality in this data, and provide your graphical analysis as evidence (you are free to present additional evidence).

    You *may* need to access specific columns of your matrix to do this analysis. If you want to do this, there are different methods. Return to the *R Bootcamp* lab to remind yourself.

    If you decide to tell `seqkit` to add a header line to your file (see option `-H`) then you will **have** to go into the file itself and edit the top line, which will probably contain a `#`. Remove this hash, re-save the file, and then `R` should read it properly.

Let's take this step by step

### Please address this question
This is the first and most important part of the assessment. *How will you test "whether GC content is correlated with quality scores for either Illumina or Oxford Nanopore sequencing reads"*?
