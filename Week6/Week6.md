**[Return to the Course Home Page](../index.html)**

# Evolutionary Visualization
**Dr Olin Silander**

## Purpose


### Preface

In this section you will use some software to do phylogenetic reconstruction and visualisation using single genes or whole genomes of SARS-CoV-2

### Introduction
[August outbreak](https://nextstrain.org/community/narratives/ESR-NZ/GenomicsNarrativeSARSCoV2/2020-10-01 "The Story")


### Learning outcomes

After studying this tutorial you should be able to:

1. Use bioinformatics software to find orthologues in the NCBI database.
2. Use bioinformatics software to perform sequence alignment.
3. Use bioinformatics software to perform phylogenetic reconstructions.

### Installing the software

Today we will be building an alignment and phylogeny. We will use several pieces of software today to get information on the open reading frames of SARS-CoV-2, and to find some other similar to perform alignments of different open reading frames and infer phylogenies.

There are a large number of ways to approach this. A standard method would be through `blast`; however, this is disabled on your RStudio system. Instead, we will explore it briefly using the web portal.

The first program we need is one that will annotate our genome. Although we have not discussed annotation in this class, it forms a key part of many bioinfomratic pipelines. Here, we will use a standard annotation method for bacterial genomes, `prokka`. Be careful, installation here may redline your RAM.

```R
# Put in some fancy bioinformatics software
install.packages("ape")
getAnnotationsGenBank(c("MN908947.3"))
```

First we must admit that we have done something terrible
```bash
bedtools genomecov -ibam kwazulu-natal.bam -bg | awk '$4 < 12'
```

Now you can annotate your reference `.fasta` file using this program.

```bash
# very simple command structure
prokka reference.fasta
# then an immediate uninstall, bye bye prokka
mamba uninstall prokka
```

The output of prokka will be a directory with a large number of files, several of which will be useful for us. One of these is the `.gbk` file, which is in *genbank* format. This file lists all the annotated reading frames (as well as tRNA, rRNA, exons, introns, etc. if this were a more complicated genome). Go ahead and navigate into the Prokka directory and use `less` to take a look into the `gbk` file. There are several standard pieces of information in this file type, and it is a far more complicated strucutre than the simpler `fasta`, `fastq`, `sam`, or `vcf` format files. We will not discuss this structure in detail, but note which pieces of information are given.

We will use this file to look at the Spike Protein (annotated in most SARS-CoV-2 genomes as the *surface glycoprotein*).

Align homologous nucleotide sequences. This program is `mafft`,  installable using `mamba`.

Second, we will install `raxml`, a phylogenetic tree inference tool, which uses
maximum-likelihood (ML) optimality criterion. This program can also be installed using`mamba`.


### Selecting a gene to build a phylogeny with the software

The first thing we need to do is select a gene that we will 
use to build a *gene tree* to infer the phylogenetic relatedness
of different *SARS-CoV-2* isolates. Today, we will use the *spike* protein. First, we will get this gene (and only this gene) from your annotation. Let's take a quick look at the location of your genes in your annotation.

```bash
# find all coding sequence locations
# CDS stand for "coding sequence" (some may disagree)
grep "CDS" my_prokka_annotation.gbk
```

This should yield a list of the coding regions. The first two are the long polyproteins in SARS-CoV-2. The third is the Spike protein itself.
We can now use the location of this gene `seqkit` to make a new `fasta` file consisting only of the nucleotide sequence of this gene. We can do this in the following way:

```bash
# make a text file with the name
# of the spike gene
echo "spike"  > spikes.txt
# use seqkit to get the sequence of this gene
seqtk subseq -r 12:13 my_prokka_annotation.ffn > spikes.fasta
```

This should yield a `fasta` file containing your gnd sequence only. Check that the file does contain the expected sequence use `cat`.


### Finding orthologues using BLAST

We would usually use |blast| to find matches to your sequence in the NCBI database by requesting a `remote` search using the *spike* gene as the query. However, the remote |blast| service is not currently available. Instead, I have made a `fasta` file of spike sequences available in the `data` directory. You will need to copy this file into your own directory.


You next need to append the sequence of the *gnd* from your own *E. coli* strain to this file, using whatever set of commands you wish/know. Most likely, this will be `cat` and the redirect `>>` (do *not* use `>` as it will overwrite the whole file).


### Performing an alignment


We will use `mafft` to perform our alignment on all the sequences in the |blast| fasta file.
This syntax is very simple (change the filenames accordingly):

```bash
mafft
```

### Building a phylogeny

We will use `iqtree` to build our phylogeny.
This uses a maximum likelihood method to infer parameters of evolution and the topology of the tree.
Again, the syntax of the command is fairly simple.

The arguments are:

- ``-s``: an alignment file
- ``-m``: a model of evolution. In this case we will use a general time reversible model with gamma distributed rates (GTR+GAMMA)
- ``-n``: outfile-name
- ``-p``: specify a random number seed for the parsimony inferences

  
```bash
raxmlHPC -s your_alignment.aln -m GTRGAMMA -n sars_tree -p 12345
```

### Visualizing the phylogeny

We will first use the online software [Interactive Tree of Life (iTOL)](http://itol.embl.de/upload.cgi "just tree it") to visualize the tree.
Navigate to this homepage. Once you have done that, open the file containing your tree, copy the contents, and paste into the web page (in the Tree text box).

You should then be able to zoom in and out to see stuff.