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
[Software Management](#software-management)<br>
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
<img src="graphics/Illumina.gif" title="It's bridge PCR!" width="300"/>

Illumina sequencing relies on sequencing-by-synthesis in which millions of single DNA molecules are multiplied into millions "clusters", and each cluster is sequenced by the incorporating fluorescent nucleotides and imaging the cluster. [Review the method here](https://www.illumina.com/science/technology/next-generation-sequencing/sequencing-technology.html "Illumina movie"). Read lengths for Illumina range between 75 bp and 300 bp, and are of very high quality (i.e. the sequence of base pairs is almost certainly correct, with an error rate of approximately 1 in 10,000).

### PacBio
<img src="graphics/pacbio.gif" title="The world's smallest microscope" width="300"/>

PacBio sequencing relies on imaging the incorporation of fluorescent nucleotides into single DNA molecules using zero-mode-waveguides ("the worlds smallest microscope"). [Review the method here](https://www.pacb.com/smrt-science/smrt-sequencing/ "PacBio movie"). Read lengths for PacBio range up to 30 kilobase pairs, and range in quality from very low (15% error rate) to very high (1 in 100,000 or lower error rate).

### Oxford Nanopore
<img src="graphics/ont.gif" title="So many pores" width="300"/>

Oxford Nanopore sequencing relies on sensing current changes in a pore as a DNA *or RNA* molecule is passed through a pore (a protein taken from *E. coli*). [Review the method here](https://nanoporetech.com/applications/dna-nanopore-sequencing "ONT movie"). Read lengths for Oxford Nanopore are essentially unlimited (e.g. 1 megabase pair), and are of medium quality, with an error rate of approximately 1%.

**Note that below, we will refer to any DNA sequence data from and NGS platform as a "read".**

#### Today

Today we will deal with DNA sequence data from two of the most widely-available technologies, Illumina and Oxford Nanopore. The primary difference between these two technolgies is that Illumina provides short, highly accurate reads using a very expensive machine (~ $1 million USD), while Oxford Nanopore provides long, less accurate reads using a very cheap machine (~ $1000 USD). We will see that these characteristics provide different advantages.

Oxford Nanopore and Illumina differ in some other ways, but we will not discuss those in detail today. Perhaps the primary difference is that Oxford Nanopore sequences the *original* molecules of DNA and RNA with all their varied modifications, whereas Illumia sequences *copies* of DNA only.
 
<img src="graphics/ont-ill.png" title="Scary car ride" width="500"/>

### Software Management

We need software to be able to process all of this data. Until now (other than Linux, the operating system itself), you have used the pre-installed statistical programming software, `R`. However, we need additional software to process DNA sequence data.\*

\* *Here, rather than just "software", we will refer to these computer programs as "software packages", as one piece of "software" contains several pieces of related software, hence software package (you have seen the phrase software packages before when installing new programs in `R`).*

As you are probably aware, software **packages** are sets of tools that have been developed to perform specific jobs, or are used to implement specific methods. Your general view of a software package may be something like Excel or Chrome or TikTok. More fundamentally, a software package is simply a group of computer programs used to perform a specific task. In bioinformatics, for example, this could be a set of programs (individual computer programs written in a language such as `python`) telling the computer how to interpret and display the quality scores from a `.fastq` file.

*However*, software packages and tools often have **dependencies**, which are other pieces of software that are necessary to run the software you would like to install. For example, to use Instagram, you also need software that controls your phone's camera. This reliance of Instagram on camera-controlling software is known as a **dependency**. Importantly, software like Instagram is designed to be **user-friendly**, and during installation will usually check whether such camera-controlling software exists, and if it does not, may try to install it.

<img src="graphics/File_dependency.png" title="About to crumble" width="300"/> <br>
**Software dependencies are real** (credit: *xkcd*) <br>


Despite the existence of dependencies, many bioinformatics software programs (most of which are written by academic-oriented computational biologists -- or worse, plain-old biologists) do not check for dependencies. This can create significant issues if you try to run a piece of software but are missing other software that it needs to run. To make sure that we resolve all these dependency issues when we install new software, we will use a **package management** system. This management system is called [conda](https://en.wikipedia.org/wiki/Conda_(package_manager "Wikipedia link"), and it is perhaps the most common package manager used in bioinformatics. It considerably simplifies the process of installing software, negating the need to find websites, download multiple files, unzip files, find compatible files for your operating system, etc.

<img src="graphics/dependencies.jpg" title="Woody and Buzz in dependency Hell" width="600"/>

### Conda Installation

As with any software, the first thing we need to do is install the package manager itself. The installation of this tool is perhaps the most complicated installation we will do in this course, as we cannot use `conda` to install itself (and I have not pre-installed it on your system). However, after the installation of `conda`, your life will become far easier (well, in terms of analysing biological data) and you will be on your way to becoming a seasoned [bioinformatician](https://soundcloud.com/microbinfie "binfie podcast").

First, I need to post a **reminder** -- as we will be operating mostly on the command line, you **must never forget** tab-complete. Also, never forget the up arrow.


<img src="graphics/uparrow.jpg" title="Not just for begginers" width="600"/>

**Even seasoned bioinformaticians use it.**

Second, try to follow the instructions exactly today. If you get an error or warning of any sort, go back *and make sure you have followed the instructions.* If you continue to get the error, then it *could* be my fault.

Good. Now, we download `conda`.

Navigate to the command line tab on your RStudio window ("Terminal"). This is on the top of the `R` window. **Make sure you are in your `/Cloud/project/` directory.**


```bash
    # Download the latest conda installer
    # we cry because we can't use tab-complete here as 
    # the file does not yet exist on our computers.
    # you should be able to copy the line below and  
    # paste it on the command line
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
**Explanation**: `wget` is a program that is used to transfer data to or from a server. Thus, this command is simply using `wget` program to find the file at the location indicated and then download it.

The file you have downloaded (with the extension `.sh`) is a bash file, which is usually run using the command line program `bash`. As you know, noting the *extension* of a file can be very helpful in figuring out what is in it, or what it does. For example, you should never end a `bash` file with `.txt` as that suggests it is a simple text file, when in fact it is not. Similarly, you would never end a Microsoft Word file with `.xlsx`, you would end it with `.doc` or `.docx`. *And* if you do find a file with the suffix `.sh` you can guess it's a `bash` file and use `bash` to run it.

### Aside: Naming Conventions

One important aspect of organising files and directories (folders) is [naming convention](https://en.wikipedia.org/wiki/Naming_convention_(programming "Wiki page on naming convention"). When working on the command line, your life will become considerably easier if you avoid using spaces in your files and directory names. Thus, **never** name your file `my awesome file.txt`. Instead, name it `my_awesome_file.txt` ("snake case"), or `myAwesomeFile.txt` ("camel case") or `my-awesome-file.txt` ("kebab case") or `my.awesome.file.txt` and probably not `MY_AWESOME_FILE.txt` ("screaming snake case") or `MY-AWESOME-FILE.txt` ("spicy kebab case"). You should pick one of these at the start of the course, and *stick to that format throughout the course* (i.e. camel case, or kebab case, etc. I usually use camel case, but not always - note that kebab case requires one less key to type. And use a `.` means that your file names will only ever have one type of non-word character. Last, you should almost **never** begin a file with a `.` (e.g. `.my-awesome-file.txt`) as this will make it a hidden file.

<img src="graphics/naming.jpg" title="Do as I say, not as I do" width="600"/>

**Please be consistent with your file names**

As I pointed out above and will re-emphasise here, the second thing to pay attention to when naming files is the *extension* or suffix. For example *text files* are usually named with the extension `.txt`. *Often*, but not always, file extensions have three characters. There are some well-known exceptions of course - `.html`, `.docx`, `.xlsx`, and others. In this course, we will run into a wide variety of files with a wide variety of extensions, for example `.fastq`, `.sam`, `.bam`, `.txt`, `.sh`, `.fasta`, `.html`, `.gbk`, `.bai`, `.py`, `.r` (sometimes `.R`, `.gz`, `.aln`, `.tre`, `.phy`, `.vcf`,  `.bcf`, and many more!  Hopefully at the conclusion of this Semester you will be familiar with all of these.

Also: if you are ever writing the date, use the format `YYYYMMDD`. That way, when you sort files by name, they will *also* be sorted by date (e.g. if you know that you made one set of data before another, it will be easier to find).

While we are the topic of [naming conventions](https://en.wikipedia.org/wiki/Naming_convention_(programming) "Wiki naminbg conventions"), there are certain characters that you should **always** avoid when naming files and folders. Besides spaces, these are (not necessarily exhaustive):


```bash
   : ; ` " ' \ / ! @ # $ % ^ & * ( ) + , ? [ ] { } | > <
```

### Back to the Topic at Hand

Let's now actually install `conda` (in our case we install a miniature version of it with less bloat, `miniconda`).

**Warning**: Be careful when using `rm` in the following command. (Why? What does `rm` do?)

```bash
    # Run the installer
    # Note: now you can use tab-complete
    # During installation. You will need to 
    # press enter and the spacebar several
    # times at the --More-- prompt, and 
    # type "yes" three times. It should
    # be readily apparent where to do this. 
    bash Miniconda3-latest-Linux-x86_64.sh
    
    # delete the installer after successful run
    rm Miniconda3-latest-Linux-x86_64.sh
```

**IMPORTANT**
Conda may not behave *quite* as it should in this cloud-based platform. Try typing `conda --help`. If there is an error, then to use `conda` we will need to addjust our `$PATH` variable.[^1] You do this by typing the following at the terminal (simply copy-paste) the *entire* line below:

```bash
export PATH="$HOME/miniconda3/bin:$PATH"
``` 
**If you keep your terminal session active (i.e. don't close it) then you should be able to use `conda` for the rest of the Semester in the terminal. If you do close the terminal, repeat the above step**.

Now, you should be able to use the `conda` command. Again, one useful way to check that `conda` (*or most other command line programs*) is working is to ask the program for help. This is **almost always** done by typing `--help` or `-h` after the command. For example try:

```bash
conda --help
``` 

This should bring up a list of sub-commands that `conda` can do (assuming you have installed it correctly). If this does not work, ask someone for help (lecturer, demonstrator, or classmate).

### Software Installation

The process of installing a software package via `conda` is called a *recipe*, and these recipes are contained in places called *channels*. Most recipes for bioinformatic software are contained in the [bioconda](https://bioconda.github.io "bioconda github") channel, which currently has recipes for more than 7000 software packages.

Let's try to install some software now. One piece of software we will need allows us to figure out where a certain sequence of DNA is in a genome (this type of software can be referred to as a "mapper" because it "maps" sequence reads). We will explore why you might want to do this later. The software we will use is [minimap2](https://github.com/lh3/minimap2 "minimap2 github"). This process is relatively simple:

```bash
# Below we specify which channel we would like conda
# to look in to find the recipe to install minimap2.
# This is done using the -c option

conda install -c bioconda minimap2
```

`conda` will sit around looking for the recipe for a minute, and then it should ask you whether you want to install the software. Simply press `enter` or type `Y` and press `enter`. Let's now get to the task at hand for today: analyzing DNA sequences from SARS-CoV-2.

### SARS-CoV-2 Genome Sequencing
In the past 2.5 years, SARS-CoV-2 (the causative agent of COVID-19) has become one of the most extensivley sequenced organisms on earth, with close to [ten million whole genome sequences available](https://www.nature.com/articles/d41586-021-03698-7 "Nature article on Omicron"), and two million sequence [by the UK alone](https://www.gov.uk/government/news/uk-completes-over-2-million-sars-cov-2-whole-genome-sequences "gov.uk press release"). SARS-CoV-2 genome sequencing is performed for two primary reasons: (1) to track the emergence of new and possibly more virulent variants, and (2) to track transmission between people. It is this second application that is (until recently) the primarily application used in New Zealand, in constrast to most other countries (again, until recently with the emergence of the Omicron variant).

#### QUESTION
1. Why has SARS-CoV-2 genome sequencing been used more commonly for transmission tracking in New Zealand compared to other countries?

Please look over [this paper here](files/sc2_flight_transmission.pdf "Strains on a plane"), especially figures 3 and 4; and [this paper here](files/sc2_realtime_genomics.pdf "COVID outbreaks in NZ") (how many of those outbreaks do you remember?) for some applications of SARS-CoV-2 genome sequencing data in New Zealand.

## Today's Data

There are several methods used to sequence SARS-CoV-2, but perhaps the most common is via [amplicon panels](https://sg.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-amplicon-sequencing/predesigned-amplicon-panels), in which PCR is used to amplify the entire genome in short pieces, which are then sequenced. The four most common methods are listed [here](https://sg.idtdna.com/pages/landing/coronavirus-research-reagents/ngs-assays#offerings "IDT SARS-CoV-2 methods"). Note, specifically, the ["xGen SARS-CoV-2 Midnight Amplicon Panel"](https://sg.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-amplicon-sequencing/predesigned-amplicon-panels/sars-cov-2-midnight-amp-panel#product-details "Midnight method") &#128540; as we will be using data generated with that method.

The sequence data that we will be using today are Illumina and Oxford Nanopore reads from two SARS-CoV-2 genomes. The format of the data is *fastq*, which specifies a name for each sequence, the sequence itself (i.e. order of basepairs), and the quality of each basepair (i.e. how certain the sequencing machine is that it is giving you the correct base). Review [fastq format here](https://en.wikipedia.org/wiki/FASTQ_format "fastq on Wikipedia").

The Illumina data are available here: [read1](./data/kwazulu-natal-2020-06-02_R1_sub.fastq.gz) and [read2](./data/kwazulu-natal-2020-06-02_R2_sub.fastq.gz) (the data are *paired end*, so there are two files). The Oxford Nanopore data are available [here](./data/montana-2021-29-09.fastq.gz).

To download the data, first make sure you are in your `/Cloud/project` directory. Second, make a new directory called `data`, and change into that directory. Third, copy the link address (right click on the link and scroll to *Copy Link Address*). Finally, download the files using `wget`:

```bash
wget link_address_you_just_copied.fastq.gz
```
You should see a rapid animated arrow tracking the download.

Repeat this for *all three sequence files*. Quick note: here and throughout the lab sessions I will often refer to certain files or directories as "myfile.txt" or "mydir/". This does not mean that you should use this name, or that this file or directory even exists for you. Rather, you should replace this name with a filename that *does* exist or which you *do* want to analyse or which is nicely descriptive. For example, intead of a ` data/` directory above, you could make a directory called `sequence_data`. Feel free to change the name now if you like.

After downloading all three files, you have all the DNA sequence data that we will use today. If you have done this correctly, you should be able to list the files from the command line. The following command should give information on which files are present in the directory, and whether they contain anything. Note that the command looks specifically for *fastq.gz* files using a *wildcard* (the asterisk) to *list* (`ls`) only those files that match that pattern. There are two additional options in the command below: `-h`, which lists file sizes in human readable format, and `-l`, which lists in long format. Here we combine them into `-lh`:

```bash
ls -lh *fastq.gz
```

#### QUESTIONS
Are all three files present?

Are you sure they all sitting in the `/data` directory that is sitting within your `/Cloud/project/` directory?

What does the `.gz` at the end of the file name indicate?

#### Organization (optional)
Things will start to get a little crazy as we get more data - new directories, new files, lousy `ls` commands. Let's see if we can look inside this maze of files in a more accessible way. How? `tree`. Let's install `tree`:

```bash
# here we actually change our channel
# i.e. it's no longer bioconda

conda install -c conda-forge tree
```

Now we can use our tree command to see what is where and how it's organised:

```bash
# only look two *levels* deep
# using the -L *option*
# first head back up a directory
cd ..
# then run tree
tree -L 2 
```

Check out the `tree` command using the `--help subcommand`. What does the `-L` option do?

Nice.

### Critically Evaluating Your Data

#### Making Good Use of Summary Statistics

First, let's next look quickly inside the files. However, we don't want to open them up - they're quite large (okay, they're not *that* large). Instead, we use a simple terminal command, `head` (you could also use `less` or `tail`. You have encountered all of them previously):
```bash
# Here we first have to use zcat, not cat, as the files are zipped. 
# We then pipe | the result to head (remember, this is the command used to 
# look at the first few lines of a file - by default, ten lines).
# You've also encountered the pipe before - it takes the output of one command
# and feeds it to another.
zcat choose_one_fastq_file_to_look_at.fastq.gz | head
```

Once we have the data and have decided that it looks as we expect, the first thing we will do is get some summary statistics (all good data science and bioinformatics and, indeed, *any science* should begin with actually *looking* at the data). Luckily, there are a number of other pieces of software that have been written to do this, so we will not need to re-invent the wheel. Today we will use two pieces of software. The first is [seqkit](https://bioinf.shenwei.me/seqkit/ "seqkit site"), a blazingly fast and flexible piece of software. Install:

```bash
# Below we use conda (of course) and 
# tell conda which *channel* to look in
# for the recipe using the -c option
conda install -c bioconda seqkit
```

If your command does not work, let a lecturer, demonstrator, or classmate know.

We will *also* use a program called [fastp](https://github.com/OpenGene/fastp "fastp github"). Install using conda:

```bash
# Once more we look in the bioconda channel
# but now I won't tell you what to type
# because you know how to install things with
# simple names like "fastp" by
# using conda by looking in a specific channel
```

#### Summary Stats with seqkit

Let's use `seqkit` first. Type `seqkit --help` to make sure it's working. No errors? If you have an error, ask for help. Next, some simple statistics about your read files:

```bash
# some simple statistics about your files
seqkit stats *fastq.gz
```

#### QUESTIONS
What does `*fastq.gz` mean in the above command?

i.e. what is the `*` doing?

Let's get the same stats but slightly fancier.

```bash
# Some much prettier stats about your sequences.
# First let's get a new package called csvtk.
# Install with conda (of course)
conda install -c bioconda csvtk

# Then some nicer looking stats but no commas.
# We use the stats subcommand the -T -t are *options*.
# I can explain those if you like. Here we also see the magical | pipe again.
# If you want you can paste this line, it's a bit long.
# We take the output of seqkit stats and format it with csvtk
seqkit stats *.fastq.gz -T | csvtk csv2md -t
```
#### QUESTION
How do the average read lengths differ between your sequencing files?

Let's look at whole distributions of read lengths instead of just the *average* read length for all read:

```bash
# Even fancier, a histogram of read lengths.
# We use the watch subcommand and the --bins option to clean up the 
# histograms a tiny bit. The --bins options tells the program
# how many different histogram bins you want (here, 15)
# If you want, you can leave the --bins 15 part of the command out.
# Note that the fastq.gz file name below is not correct.
seqkit watch --bins 15 choose_one_fastq_file_to_plot.fastq.gz
```
Use this `seqkit watch` command for all of your sequencing files. You can also try the `--log` option if you want (what does this option do?).

#### QUESTION
How do the sequencing files differ in the *distributions* of read lengths?

It is also possible to make a simple plot of the average *quality* of each read. In this case, quality of a base (A,C,G,T) in a read refers to the likelihood that the base is incorrect. See the explanation of quality scores [here](https://en.wikipedia.org/wiki/FASTQ_format#Quality "Wikipedia quality scores"). Recall that Illumina and Oxford Nanopore data differ in their read accuracy, and thus quality. Which technology has higher accuracy? Go ahead and plot the quality scores. Here we also use `seqkit` to plot the mean quality of all base pairs in a read (for all reads):

```bash
# below we add a new argument, --fields, to specify which 
# read aspect we would like to plot. Above it was MeanLength
# Now it's MeanQual. We leave the --bins option in (you don't have to))
# When you leave in --bins you need to specify a number (I chose 15)
seqkit watch --fields MeanQual --bins 15 choose_one_fastq_file_to_plot.fastq.gz
```

Do this for both the Oxford Nanopore and Illumina reads.

#### QUESTIONS
1. If a base has a quality score of 10, what is the likelihood that it is the correct base?
2. What if a base has a quality score of 25?
3. How do the sequencing files differ in the distributions of quality scores?
4. What are the *highest* average read qualities for the Oxford Nanopore reads (approximately)?
5. For an Oxford Nanopore read with a "high" average quality, what fraction of bases do you expect to be correct?

#### Summary Stats with fastp

Let's use the somewhat aesthetically more pleasing `fastp` program next. You should have installed it above, check if it's properly installed:

```bash
# as usual, we use the --help option to see if the
# install worked.
fastp --help
```

```bash
# Also use the fastp --help to see the options. There are four important 
# parameters, the input files (-i), output (-o), and the names of the 
# html and json files (where the summary information goes). We have to
# specify all of those. You probably know what a .html file is. You may not 
# know what a .json file is. Do not worry about it, I barely do but look it 
# up if you want. I add the -A options below (please do as well). What does it do?
fastp -A -i montana-2021-29-09.fastq.gz -o montana-2021-29-09.trim.fastq.gz -h montana.html -j montana.json
```

```bash
# repeat for the Illumina data
fastp -i kwazulu-natal-2020-06-02_R1_sub.fastq.gz \
-I kwazulu-natal-2020-06-02_R2_sub.fastq.gz \
-o kwazulu-natal-2020-06-02_R1_sub.trim.fastq.gz \
-O kwazulu-natal-2020-06-02_R2_sub.trim.fastq.gz
```

Unfortunately, the results of the fastp analysis do not display properly in your browser. However, I think it is important that you see them. For this reason, do the following:
1. Find a simple text editor on your computer. For MacOS users, this could be TextEdit. For windows users, this could be Notepad.
2. Open a new document in the editor.
3. In your RStudio _files_, find the files that you just made using fastp (e.g. montana.html and montana.json). open the `.html` _in the RStudio editor window_.
4. Copy all the text from the editor window,  paste it into the text editor window, and save this file to your desktop with a .html suffix.
5. Repeat this for the `.json` file and save it with _the same name as the html file but with a `.json` suffix_.
6. Open the `.html` file now on your desktop by double clicking on it. This shoudl bring you to a web page displaying the `fastp` results.

#### QUESTIONS

### Getting Tabulated Data to Plot in R

In the fastp report we have several useful statistics, such as the mean read quality for each base over the length of the reads, and the fraction of bases at each position that are A, C, G, or T. However, there are additional important statistics that are not provided. Some of these we saw with the `seqkit` tool kit. However, the plots were not aesthetically pleasing, and in some cases difficult to interpret. For that reason, we are going to replot some data, as well as additional data, using our familiar plotting package, `R`.

To begin, we make some files with the new data we would like. In this case, we will use `seqkit` again. First, some data on the distribution of lengths and quality _per sequence_ (rather than per position). Do this for your untrimmed and trimmed Oxford Nanopore data (Montana):

```bash
# -q and -l options give the length and quality for each sequence.
# -n supresses the output of the sequences and qualities of *each bp* 
# of each read
seqkit fx2tab -qln myseqs.fastq.gz
```

Oops. Did you forget that `seqkit` outputs to standard out? You have to use the redirect `>` to put it into a file.

```bash
# -q and -l options give the length and quality for each sequence.
seqkit fx2tab -qln myseqs.fastq.gz > mydata.txt
```

#### QUESTION
Why have I used `.txt` at the end of the example file above?

### Plotting the Data more Deliberately in R
Now we have to think.

What types of plots are we interested in?

What sort of information would we like to see?

What sort of information *is it possible for us to see*?

**What sort of story would we like our plots to tell**?

### Take Home Messages

### Notes

[^1]: The `$PATH` variable contains places (directories) in which your computer looks for programs. These directories are listed one after the other. The computer will search these in the order they are listed until the program you requested is found (or not, then it will complain). For example, you might have a `$PATH` variable that says: first look in my home directory (`~/`), and then in the `/usr/bin/` directory, and then in my friend's directory (`/friends_dir/sneaky_software_i_saved_there/`). However, those are *the only* places the computer will look. If you want the computer to look in more places, you have to add those locations to the `$PATH` variable. The `$` indicates that it is a *variable*.

    After installation, `conda` would usually adjust your `$PATH` variable by adding the directory `~/miniconda3/bin` to your `$PATH` variable -- so that the program `conda` can be found anytime you open a new shell, and any program that `conda` installs will be used first because the computer will look in this directory first. *However*, `RStudioCloud` does not recognise the addition of `conda` to the path, so we add it manually.


##### scratch below
each of which are more oriented toward a specific sequencing platform. The first of these is [NanoPlot](https://github.com/wdecoster/NanoPlot "NanoPlot github").

First things first, install `Nanoplot`:

```bash
# again we use conda and look in the bioconda channel
conda install -c bioconda nanoplot
```

We will also use a program called [pycoqc](https://github.com/a-slide/pycoQC "pycoqc github")
