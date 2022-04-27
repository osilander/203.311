**[Return to the Course Home Page](../index.html)**

# 27-Apr 2022 - PLEASE DO NOT MAKE A COPY OF THIS FILE UNTIL THIS LABEL IS REMOVED

# Barcodes and diversity - Fastq read analysis for 16S rRNA metabarcoding using the R package DADA2

**Professor Patrick Biggs**

[Purpose](#purpose)<br>
[Introduction](#introduction)<br>
[Lecture Overview](#lecture-overview)<br>
[Conventions used for this RStudio practical](#conventions-used-for-this-rstudio-practical)<br>
[Accessing the resources needed](#accessing-the-resources-needed)<br>
[Theoretical overview](#theoretical-overview)<br>
[Exercise 1: Getting everything ready](#exercise-1-getting-everything-ready)<br>
[Exercise 2: Collecting our data](#exercise-2-collecting-our-data)<br>
[Exercise 3: Examining the quality profiles of forward and reverse reads](#exercise-3-examining-the-quality-profiles-of-forward-and-reverse-reads)<br>
[Exercise 4: Performing filtering and trimming](#exercise-4-performing-filtering-and-trimming)<br>
[Exercise 5: Dereplication](#exercise-5-dereplication)<br>
[Exercise 6: Investigation of the error rates](#exercise-6-investigation-of-the-error-rates)<br>
[Exercise 7: Sample inference](#exercise-7-sample-inference)<br>
[Exercise 8: Merging the paired reads](#exercise-8-merging-the-paired-reads)<br>
[Exercise 9: Constructing the sequence table and removing chimaeras](#exercise-9-constructing-the-sequence-table-and-removing-chimaeras)<br>
[Exercise 10: Checking our progress](#exercise-10-checking-our-progress)<br>
[Exercise 11: Assigning taxonomy](#exercise-11-assigning-taxonomy)<br>
[Exercise 12: Evaluating the accuracy](#exercise-12-evaluating-the-accuracy)<br>
[Exercise 13: Working with fasta files from our analysis](#exercise-13-working-with-fasta-files-from-our-analysis)<br>
[Portfolio analysis](#portfolio-analysis)<br>
[Assessment](#assessment)<br>
[Contact](#contact)<br>



## Purpose

To learn how to perform an analysis of 16S rRNA sequencing data from raw sequencing data to preliminary visualisation and taxonomic classification.



## Introduction

Microbial diversity has been revolutionised by the rapid advances in sequencing technology, given us new insights into the role of the microbial world in all environments on Earth. Sequencing of metagenomes (all the DNA in a given sample etc) or of specific markers (metabarcodes e.g. 16S rRNA amplicons) have been crucial in this regard. Knowing what is there, and what it is potentially doing are two of the main research areas of interest.  How we then analyse the data, visualise and interpret the results is thus of importance.

Week 7 focusses on introducing ideas around the analysis of microbial diversity - for bacteria - within the 16S rRNA amplicon, and performing a typical workflow analysis on such sequences in R using a piece of software called DADA2. These ideas will be developed further in the next two weeks (weeks 8 and 9) of this module. 
 










## Assessment

To reiterate, there is no direct assessment today.  What is required however, is an understanding of the principles we have learnt today, as these will be required for the mastery test which accounts for 15% of the course.  This will take place between Thursday 19-May-2022 and Friday 20-May-2022 online.

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


[^1]: <https://www.nature.com/nmeth/journal/v13/n7/pdf/nmeth.3869.pdf >
[^2]: Concepts for this part of the practical came from a tutorial found at: <https://benjjneb.github.io/dada2/tutorial.html>
