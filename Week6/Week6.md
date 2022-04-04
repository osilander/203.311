**[Return to the Course Home Page](../index.html)**

# Phylogenetics and Evolutionary Visualization
**Dr Olin Silander**

## Purpose


### Preface

In this section you will use some software to do phylogenetic reconstruction and visualisation using single genes or whole genomes of SARS-CoV-2.

### Introduction



[August outbreak](https://nextstrain.org/community/narratives/ESR-NZ/GenomicsNarrativeSARSCoV2/2020-10-01 "The Story")


### Learning outcomes

After studying this tutorial you should be able to:

1. Use bioinformatics software to find orthologues in the NCBI database.
2. Use bioinformatics software to perform sequence alignment.
3. Use bioinformatics software to perform phylogenetic reconstructions.


Today we will be building an alignment and phylogeny. We will use several pieces of software today to get information on the open reading frames of SARS-CoV-2, and to find some other similar to perform alignments of different open reading frames and infer phylogenies.

There are a large number of ways to approach this. A standard method would be through `blast`; however, this is disabled on your RStudio system. Instead, we will explore it briefly using the web portal.

First we must admit that we have done something terrible
```bash
bedtools genomecov -ibam kwazulu-natal.bam -bg | awk '$4 < 12'
```


### Installing the software

The steps below should be done in the `R` console.

```R
# Put in some fancy bioinformatics software
install.packages("ape")
# an example of ehat it can do for us.
getAnnotationsGenBank(c("MN908947.3"))
```

The output of the above command is a list of the annotations of the ancestral SRAS-CoV-2 genome. Most often, annotated genomes are given in Genbank format, usually suffixed with `.gbk` file, which is in *genbank* format. This file lists all the annotated reading frames (as well as tRNA, rRNA, exons, introns, etc. if this were a more complicated genome). Click on this link [here](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3 "Ancestral Genbank") to see what this format looks like. Note that it is considerably more complicated than any other format we have seen so far (`.sam`, `.fastq`, `.fasta`, `.vcf`, `.sh`, and the associated `.fai`, `.bam`, `.bai`, `.bcf`)

We will use this file to look at the Spike Protein (annotated in most SARS-CoV-2 genomes as the *surface glycoprotein*).

Align homologous nucleotide sequences. This program is `mafft`,  installable using `mamba`.

Second, we will install `iqtree`, a phylogenetic tree inference tool, which uses
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
iqtree 
```

### Visualizing the phylogeny

We will first use the online software [Interactive Tree of Life (iTOL)](http://itol.embl.de/upload.cgi "just tree it") to visualize the tree.
Navigate to this homepage. Once you have done that, open the file containing your tree, copy the contents, and paste into the web page (in the Tree text box).

You should then be able to zoom in and out to see stuff.