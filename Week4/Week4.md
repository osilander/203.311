**[Return to the Course Home Page](../index.html)**

# Sequencing and Mapping
**Dr Olin Silander**

[Purpose](#purpose)<br>
[Introduction](#introduction)<br>
[Background](#background)<br>
[SARS-CoV-2 Genome Sequencing](#sars-cov-2-genome-sequencing)<br>
[Illumina](#illumina)<br>
[PacBio](#pacbio)<br>
[Oxford Nanopore](#oxford-nanopore)<br>
[The Data](#data)<br>
[Making Good Use of Summary Statistics](#making-good-use-of-summary-statistics)<br>
[Choosing A Plot Type](#choosing-a-plot-type)<br>
[Critically Evaluating Your Data](#critically-evaluating-your-data)<br>
[Take Home Messages](#take-home-messages)<br>


## Purpose

1. To learn how to use a package manager to install software for use on the command line.
2. To understand the advantages and disadvantages of using two of the three most common types of NGS sequencing data.
3. To be able to visualise the differences in these NGS data types.
4. To discuss what the possible applications are for this type of NGS data.

## Introduction

The data we will investigate today is from publicly available SARS-Cov-2 genome sequences. Over the next two weeks. you will learn how to use this data to find specific mutations that new SARS-CoV-2 strains have, how to place these strains in an evolutionary context, and how to visualise this context.

## Background

Soon after the birth of Next Generation Sequencing in 2005 (or so), the technology rapidly proliferated into a number of 
[different platforms](https://en.wikipedia.org/wiki/Massive_parallel_sequencing "NGS sequencing platforms") (e.g. 454, IonTorrent, Helicos, and others). However, this number has been gradually weaned down, and currently there are three dominant NGS sequencing platforms: [Illumina](https://en.wikipedia.org/wiki/Illumina,_Inc. "Illumina on Wikipedai"), which dominates the market; [PacBio](https://en.wikipedia.org/wiki/Pacific_Biosciences "PacBio on Wikipedia"); and [Oxford Nanopore](https://en.wikipedia.org/wiki/Oxford_Nanopore_Technologies "Oxford Nanopore on Wikipedia").
These three technologies differ considerably in their methodologies. For all three, sequencing output ranges from tens of gigabases (billions of nucleotides) to terabases (trillions of nucleotides), depending on the specific platform (e.g. Illumina MiSeq, Illumina NovaSeq, Oxford Nanopore MinION, Oxford Nanopore PromethION, etc.).


### Illumina
<img src="graphics/Illumina.gif" width="300"/>

Illumina sequencing relies on sequencing-by-synthesis in which millions of single DNA molecules are multiplied into millions "clusters", and each cluster is sequenced by the incorporating fluorescent nucleotides and imaging the cluster. [Review the method here](https://www.illumina.com/science/technology/next-generation-sequencing/sequencing-technology.html "Illumina movie"). Read lengths for Illumina range between 75 bp and 300 bp, and are of very high quality (i.e. the sequence of base pairs is almost certainly correct, with an error rate of approximately 1 in 10,000).

### PacBio
<img src="graphics/pacbio.gif" width="300"/>

PacBio sequencing relies on imaging the incorporation of fluorescent nucleotides into single DNA molecules using zero-mode-waveguides ("the worlds smallest microscope"). [Review the method here](https://www.pacb.com/smrt-science/smrt-sequencing/ "PacBio movie"). Read lengths for PacBio range up to 30 kilobase pairs, and range in quality from very low (15% error rate) to very high (1 in 100,000 or lower error rate).

### Oxford Nanopore
<img src="graphics/ont.gif" width="300"/>

Oxford Nanopore sequencing relies on sensing current changes in a pore as a DNA *or RNA* molecule is passed through a pore (a protein taken from *E. coli*). [Review the method here](https://nanoporetech.com/applications/dna-nanopore-sequencing "ONT movie"). Read lengths for Oxford Nanopore are essentially unlimited (e.g. 1 megabase pair), and are of medium quality, with an error rate of approximately 1%.

Note that below, we will refer to any DNA sequence data from and NGS platform as a "read".

Today we will deal with DNA sequence data from two of the most widely-available technologies, Illumina and Oxford Nanopore. The primary difference between these two technolgies is that Illumina provides short, highly accurate reads using a very expensive machine (~ $1 million USD), while Oxford Nanopore provides long, less accurate reads using a very cheap machine (~ $1000 USD). We will see that these characteristics provide different advantages.

<img src="graphics/ont-ill.png" width="500"/>

### Software Management

We need software to be able to process all of this data. Until now, you have used the pre-installed statistical programming software, `R`. However, we need additional software to process DNA sequence data. However, rather than just "softare", we will refer talk about "software **packages**", as one piece of "software" often *actually* contain several pieces of related software, hence *software packages* (you have seen this word before when installing new software in `R`).

As you are probably aware, software **packages** are sets of tools that have been developed to perform specific jobs, or are used to implement specific methods. Your general view of a software package may be something like Excel or Chrome or TikTok. More fundamentally, software is simply a group of instructions used to perform a specific task. In bioinformatics, for example, this could be a set of instructions telling the computer how to interpret and display the quality scores from a `.fastq` file.

*However*, software packages and tools often have **dependencies**, which are other pieces of software or tools that are necessary to run the software you would like to install. For example, to use Instagram, you also need software that controls your phone's camera. This reliance of Instagram on camera-controlling software is known as a **dependency**. Importantly, software like Instagram is designed to be **user-friendly**, and during installation will usually check that such camera-controlling software exists, and if it does not, may try to install it.

Despite the existence of dependencies, many bioinformatics software programs, most of which is written by inexperienced computer scientists (or worse, biologists) do not check for dependencies. This can create significant issues if you try to run a piece of software but are missing dependencies.

<img src="graphics/dependencies.jpg" width="600"/>

To make sure that we resolve all these dependency issues when we install new software, we will use a package managment system. This management system is called `conda`, and it is perhaps the most common package manager used in bioinformatics.

As with any software, the first thing we need to do is install it. The installation of this tool is perhaps the most complicated installation we will do in this course, as we cannot use `conda` to install itself. However, after the installation of `conda`, your life will become far easier (well, in terms of analysing biological data) and you will be on your way to becoming a seasoned [bioinformatician](https://soundcloud.com/microbinfie "binfie").

First, navigate to the command line tab on your RStudio window ("Terminal"). This is on the top of the `R` window.

Next, I need to post a **reminder** you **must never forget** tab-complete. Also, never forget the up arrow. Good. Now, we download `conda`.

```bash
    # download the latest conda installer
    # we cry because we can't use tab-complete here as 
    # the file does not yet exist on our computers.
    # you should be able to copy this line and paste 
    # the whole thing on the command line
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

This file (with the extension `.sh`) is a bash file, which is usually run using the command line program `bash`. Noting the *extension* of a file can be very helpful in figuring out what is in it, or what it does. For example, you should never end a `bash` file with `.txt` as that suggests it is a simple text file, when in fact it is not. Similarly, you would never end a Microsoft Word file with `.xlsx`, you would end it with `.doc` or .docx.

Let's now actually install `conda` (in our case we install a miniature version of it with less bloat, `miniconda`).

**Warning** Be careful when using `rm` in the following command. (Why? What does `rm` do?)

```bash
    # run the installer
    # note: now you can use tab-complete
    # you will need to press enter and spacebar several
    # times, and type "yes" twice 
    bash Miniconda3-latest-Linux-x86_64.sh
    
    # delete the installer after successful run
    rm Miniconda3-latest-Linux-x86_64.sh
```

**IMPORTANT**
Conda does not behave *quite* as it should in this cloud-based platform. For that reason, *every* time we use `conda` we will need to addjust our `$PATH` variable[^1]. You do this by typing the following at the terminal (simply copy-paste):

```bash
export PATH="$HOME/miniconda3/bin:$PATH"
``` 
**If you keep your terminal session active (i.e. don't close it) then you should be able to use `conda` for the rest of the Semester in the terminal. If you do close the terminal, repeat the above step**.

Now, you should be able to use the `conda` command. One useful way to check that `conda` (*or most other command line programs*) is working is to ask the program for help. This is **almost always** done by typing `--help` or `-h` after the command. For example try:

```bash
conda --help
``` 

This should bring up a list of sub-commands that `conda` can do (assuming you have installed it correctly). If this does not work, ask someone for help (lecturer, demonstrator, or classmate).

### Software Installation

The process of installing a software package via `conda` is called a *recipe*, and these recipes are contained in places called *channels*. Most recipes for bioinformatic software are contained in the [bioconda](https://bioconda.github.io "bioconda") channel, which currently has recipes for more than 7000 software packages.

Let's try to install some software now. One piece of software we will need allows us to figure out where a certain sequence of DNA is in a genome. We will explore why you might want to do this later. The software we will use is [minimap2](https://github.com/lh3/minimap2 "minimap2 github"). This process is relatively simple:

```bash
# Below we specify which channel we would like conda
# to look in to find the recipe to install minimap2.
# This is done using the -c option

conda install -c bioconda minimap2
```

`conda` will sit around looking for the recipe for a minute, and then it should ask you whether you want to install the software. Simply press `enter` or type `Y` and press `enter`. Let's now get to the task at hand for today: DNA sequences from SARS-CoV-2.

### SARS-CoV-2 Genome Sequencing
In the past 2.5 years, SARS-CoV-2 (the causative agent of COVID-19) has become one of the most extensivley sequenced organisms on earth, with more than five million whole genome sequences available. SARS-CoV-2 genome sequencing is performed for two primary reasons: (1) to track the emergence of new and possibly more virulent variants, and (2) to track transmission between people. It is this second application that is (until recently) the primarily application used here in New Zealand, in constrast to most other countries (again, until recently with the emergence of the Omicron variant).

#### QUESTION 1
Why has SARS-CoV-2 genome sequencing been used more commonly for transmission tracking in New Zealand compared to other countries?

Please look over [this paper here](files/sc2_flight_transmission.pdf), especially figures 3 and 4; and [this paper here](files/sc2_realtime_genomics.pdf) (how many of those outbreaks do you remember?) for some applications of SARS-CoV-2 genome sequencing data in New Zealand.

## Today's Data

There are several methods used to sequence SARS-CoV-2, but perhaps the most common are via [amplicon panels](https://sg.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-amplicon-sequencing/predesigned-amplicon-panels), in which PCR is used to amplify the entire genome, which is then sequenced. The four most common methods are listed [here](https://sg.idtdna.com/pages/landing/coronavirus-research-reagents/ngs-assays#offerings "IDT SARS-CoV-2 methods"). Note, specifically, the ["xGen SARS-CoV-2 Midnight Amplicon Panel"](https://sg.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-amplicon-sequencing/predesigned-amplicon-panels/sars-cov-2-midnight-amp-panel#product-details "Midnight method") :stuck_out_tongue_winking_eye: as we will be using data generated with that method.

The data that we will be using today are Illumina and Oxford Nanopore reads from two SARS-CoV-2 genomes. The format of the data is *fastq*, which specifies a name for each sequence, the sequence itself (i.e. order of basepairs), and the quality of each basepair (i.e. how certain the sequencing machine is that it is giving you the correct base). Review [fastq format here](https://en.wikipedia.org/wiki/FASTQ_format "fastq on Wikipedia").

The Illumina data are available here: [read1](./data/kwazulu-natal-2020-06-02_R1_sub.fastq.gz) and [read2](./data/kwazulu-natal-2020-06-02_R2_sub.fastq.gz) (the data are *paired end*, so there are two files). The Oxford Nanopore data are available [here](./data/montana-2021-29-09.fastq.gz).

To download the data, click on one of the links above to reach a page linking to the data. It should look like this:

<img src="graphics/fastq_download.png" width="600"/>

Right click the "Download" button and scroll to *copy link address*. Then navigate back to the `Terminal`, and change into your home directory (how?).

Quick note: here and throughout the lab sessions I will often refer to certain files or directories as "myfile.txt" or "mydir/". This does not mean that you should use this name, or that this file or directory even exists for you. Rather, you should replace this name with a filename that *does* exist or which you *do* want to analyse.

Once you have changed into your home directory, make a new directory called `data` (how?). Finally, change into that directory and download the data using *wget*:

```bash
wget https://the_data_file_address_you_just_copied
```

You should see a rapid animated arrow tracking the download.

**Explanation**: `wget` is a program that is used to transfer data to or from a server on the command line. Thus, this command is simply using this program to find the file at the location indicated and then download it.

Repeat this process for all three of the files above. Now you have all the DNA sequence data that we will use today. If you have done this correctly, you should be able to list the files from the command line. The following command should give information on which files are present and whether they contain anything. Note that the command looks specifically for *fastq.gz* files using a *wildcard* (the asterisk) to *list* (`ls`) only those files that match that pattern. There are two additional options in the command below, `-h`, which lists file sizes in human readable format, and `-l`, which lists in long format. Here we combine them into `-lh`:

```bash
ls -lh *fastq.gz
```

#### QUESTIONS
Are all three files present?

Are you sure they all sitting in the `/data` directory that is sitting within your `/home` directory?

What does `.gz` indicate?

##### Organization
Things are starting to get a little crazy. New directories, new files, lousy `ls` commands. Let's see if we can look inside this maze of files in a more accessible way. How? `Tree`. Let's intall `tree`:

```bash
# here we actually change our channel

conda install -c conda-forge tree
```

Now we can use our tree command to see what is where and how it's organised:

```bash
# only look two *levels* deep
cd
tree -L 2 
```

Wow. Nice.


Finally, let's next look quickly inside the files. However, we don't want to open them up - they're quite large. Instead, we use a simple terminal command. You may have encountered it previously
```bash
# here we have to use zcat as the files are zipped. 
# We then pipe | the result to head (this is the command used to 
# look at the first few lines of a file - by default, ten lines).
# If you're interested in what a pipe does I can explain it.
zcat choose_one_fastq_file_to_look_at.fastq.gz | head
```


### Making Good Use of Summary Statistics

Once we have the data, the first thing we will do is get some summary statistics (all good data science and bioinformatics and, indeed, *any science* should begin with actually *looking* at the data). Luckily, there are a number of other pieces of software that have been written to do this, so we will not need to re-invent the wheel. Today we will use two pieces of software. The first is [seqkit](https://bioinf.shenwei.me/seqkit/ "seqkit site"), a blazingly fast and flexible piece of software. Install:

```bash
# below we use conda (of course) and 
# tell conda which *channel* to look in
# for the recipe using the -c option
conda install -c bioconda seqkit
```

If your command does not work, let a lecturer, demonstrator, or classmate know.

We will *also* use a program called [fastp](https://github.com/OpenGene/fastp "fastp github"). Install:

```bash
# once more we look in bioconda
# but now I won't tell you what to type
# because you know how to install things 
# using conda by looking in a specific channel
```

Let's use `seqkit` first. Type `seqkit --help` to make sure it's working. No errors? If you have an error, ask for help.

```bash
# some simple statistics about your files
seqkit stats *fastq.gz
```

#### QUESTIONS
What does `*fastq.gz` mean?
What is the `*` doing?

```bash
# some much prettier stats about your sequences
# first let's get a new package called csvtk
# install with conda (of course)
conda install -c bioconda csvtk

# then some nicer looking stats but no commas
# all the -T -t are *options*. I can explain
# those too. Here we also see the magical | pipe again.
seqkit stats *.fastq.gz -T | csvtk csv2md -t
```

```bash
# even fancier, a histogram of read lengths
# we use the --bins option to clean up the 
# histograms a tiny bit.
seqkit watch --bins 15 choose_one_fastq_file_to_plot.fastq.gz
```

Oops the histrograms are still a bit messy :grimacing:



### Choosing A Plot Type

### Critically Evaluating Your Data

### Take Home Messages

### Notes
Besides now having |conda| available as your package manager, one additional thing has changed - your `$PATH` variable. |conda| has put the directory `~/miniconda3/bin` first on your `$PATH` variable. (The little `~` (tilde) at the start is short-hand for your home directory.) What is your `$PATH` variable and **why** do we need to prepend it with `~/miniconda3/bin`? Read on:

[^1]The `$PATH` variable contains places (directories) in which your computer looks for  programs. These directories are listed one after the other. The computer will search these in the order they are listed until the program you requested is found (or not, then it will complain). For example, you might have a `$PATH` variable that says: first look in my home directory (`~/`), and then in the `/usr/bin/` directory, and then in my friend's directory (`/friends_dir/sneaky_software_i_saved_there/`). However, those are *the only* places the computer will look. If you want the computer to look in more places, you have to add those locations to the `$PATH` variable. The `$` indicates that it is a *variable*.

The installation of `conda` shaould have told the computer to also look in `~/miniconda3/bin` - so that the program `conda` can be found anytime you open a new shell, and any program that `conda` installs will be used first. *However*, `RStudioCloud` does not recognise the addition of `conda` to the path, so we add it manually.


##### scratch below
each of which are more oriented toward a specific sequencing platform. The first of these is [NanoPlot](https://github.com/wdecoster/NanoPlot "NanoPlot github").

First things first, install `Nanoplot`:

```bash
# again we use conda and look in the bioconda channel
conda install -c bioconda nanoplot
```

We will also use a program called [pycoqc](https://github.com/a-slide/pycoQC "pycoqc github")
