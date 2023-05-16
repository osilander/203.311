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

You need to understand how to break this question down to its most important parts:
1. GC content
2. Correlated
3. Quality Scores
4. Illumina
5. Oxford Nanopore

These are the components of your answer. We need to interpret each a bit further.

### Quality Scores
Quality scores of *what*?
There are two hints in the assessment. First, it ask: *[is] GC content correlated with quality scores for either Illumina or Oxford Nanopore sequencing READS.*
With some careful reading it becomes clear we are talking about sequencing reads.

There is also a second clue: we should use `seqkit fx2tab`. What is the output of `seqkit fx2tab`?

As noted when you did this in the lab, it is "lengths and quality per sequence," or in the comments of commands, "a file with three columns: the NAME, the LENGTH, and the QUALITY".

With some careful reading we should be able to see that we need to get the **quality scores of individual reads**. Furthermore, we know how to do this: using `seqkit fx2tab`, and we know that this will give us column formatted data.

We should remind ourselves **how** to get this type of data. We can look back at the lab, where we see `seqkit fx2tab -qln myseqs.fastq.gz > mydata.txt` gives a file with that information. However, arguably a better way is through using seqkit itself. Two ways to do this - command line help:

```bash
# what is the output of seqkit help
seqkit --help

SeqKit -- a cross-platform and ultrafast toolkit for FASTA/Q file manipulation

Version: 0.13.2

Author: Wei Shen <shenwei356@gmail.com>

Documents  : http://bioinf.shenwei.me/seqkit
Source code: https://github.com/shenwei356/seqkit
Please cite: https://doi.org/10.1371/journal.pone.0163962

Usage:
  seqkit [command]

Available Commands:
  amplicon        retrieve amplicon (or specific region around it) via primer(s)
  bam             monitoring and online histograms of BAM record features
  common          find common sequences of multiple files by id/name/sequence
  concat          concatenate sequences with same ID from multiple files
  convert         convert FASTQ quality encoding between Sanger, Solexa and Illumina
  duplicate       duplicate sequences N times
  faidx           create FASTA index file and extract subsequence
  fish            look for short sequences in larger sequences using local alignment
  fq2fa           convert FASTQ to FASTA
  fx2tab          convert FASTA/Q to tabular format (with length/GC content/GC skew)
  ```

  OR [Google seqkit help] (https://bioinf.shenwei.me/seqkit/usage/ "Google: seqkit help")

Unfortunataely Googling [seqkit help fx2tab](https://bioinf.shenwei.me/seqkit/usage/) does not give us more useful information. Nevertheless we can use command-f (control-f on Windows) to find the relevant section. Finding `fx2tab` or `GC` brings us to the [fx2tab page](https://bioinf.shenwei.me/seqkit/usage/#fx2tab-tab2fx).

We can either command-f quality, which brings us to the *flag* we need to use to get quality, or we can scan the page to see that the second flag listed gives quality. But we don't actually need to any of that, as this was all noted in the previous lab.

### GC content
We next need to find **per read** GC content. We have not encountered this before. There are two solutions here. First, we have only used a single method to get pre-read characteristics, `fx2tab`. We can make an educated guess that we should use this. Alternatively, reading the Portfolio Assessment carefully, we see that we can get data on GC: "new data on GC content and quality using the `seqkit fx2tab`" using the help function.

### Illumina and Nanopore
This should be apparent - we simply need to make sure that we do this for Illumina and Nanopore.

### Correlation
We need to know what we are correlating. This is apparent in the question: "whether GC content is correlated with quality scores". We can only test this if we have the data. You should now have the data *in column format* by using the correct `seqkit f2tab` command. The only way that we have used to visualise so far is `R`, so we will stick with that.


### Correlation in R
The Portfolio Assessment does not require any statistics (although it does not note this explicitly.)

#### Getting the data into R

This is the first hurdle, and we will have to return to different parts of the labs to do this. In the Section "Plotting the data more deliberately in R" we noted that:
"you will need to use the read.table function."

In Lab 3 on data visualisation, we noted that:

