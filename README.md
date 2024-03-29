# Metagenomics Pipeline

Author: Katie McCauley

Date: August 23, 2022

# Table of Contents

- [Introduction](#introduction)
  - [Quality Control](#quality-control)
    - [QC Methods](#qc-methods)
    - [Running the QC Script](#running-the-qc-script)
  - [Reads Based Metagenomics](#reads-based-metagenomics)
    - [Methods and Included Software](#methods-and-included-software)
    - [Running Reads-Based Analysis](#running-reads-based-analysis)
    - [Performing Analysis]($performing-analysis)
  - [Assembly: In Progress/In Beta](#assembly-pipeline)
    - [Software Used](#software-used)
- [Bug Reports and Roadmap](#bug-reports-and-roadmap)

# Introduction
This document will describe how the metagenomics pipeline is run and give a sense of what is contained within the pipeline itself. This document will encompass the information I provide in the live workshops I've given to date. It will also be a live document which will change as the pipeline changes (the pipeline is not static).

The pipeline has been built to maximize efficiency as much as possible. This means that you may need to run the pipeline once, review your results, and run the pipeline again. However, I hope that despite the length of this document, you will notice that the actual effort required to do so is quite minimal. It was developed so that minimal, if any, hands-on manipulation or knowledge of bash scripting would be needed on the user side. There are likely to be outlier instances (ie, really low-quality sequencing run) where we need to start thinking about modifying some parameters, but unless you have a known outlier case, no major edits should be needed.

This pipeline has also been validated on RNAseq data. It uses bbmap which is a splice-aware aligner, so both DNA and RNA can be aligned to the human genome and filtered out. Some discussion exists that suggests that bbmap is imperfect, but alas... A future iteration of the pipeline will likely include a second human filtering step that will use bowtie2 for DNA sequences and STAR for RNA. This will likely be another parameter to provide to the script, but we'll get there when we get there.

## Quality Control

### QC Methods

The steps in this section of the pipeline perform the following:

1. Concatenate reads across lanes, if needed

2. Run FastQC to get quality information about your sequenced data. The settings for quality filtering are a tad bit higher than typically recommended, so there shouldn't be too much of a need to

3. Use [bbTools](https://jgi.doe.gov/data-and-tools/bbtools/) for:

  a. Adapter Trimming

bbTools provides a file containing the most common adapters used in metagenomimcs sequencing, so reads are mapped against that database to identify cutpoints for adapter removal.

  b. PhiX Removal
This step checks for the presence of the PhiX genome, which is typically added to the sequencing process in order to increase sequence complexity. It's also a good gut check that the data was sequenced as evenly as expected (if you put in 40% PhiX, you would expect to recover 40% PhiX). In parallel, checking that a sample that is not expected to have PhiX did indeed not have hits to PhiX is great confirmation.

  c. Quality Filtering and Trimming

Current recommended filtering protocols suggest trimming at Q-scores (based on the Phred scale) less than 10. I have modified this to enforce trimming at a quality score of 15 and an average quality of 20. In addition, a read must have a minimum length of 50 after trimming, or else it is discarded. This latter component ensures that the smallest sequence you have is actually rather large (in comparison). Some data I've processed previously resulted in reads with 15 or fewer bases, which ultimately isn't useful.

  d. Human (host) Filtering

Before moving on with microbial profiling, human reads should be discarded (very few metagenomics pipelines suggest maintaining human reads). bbMap is used for this step, which is a splice-aware aligner, meaning that you can apply this step to both metagenomics (DNA) and metatranscriptomics (RNA) data. There is some discussion of how accurate bbMap is in identifying reads to filter, and some reads *do* get by the filter, but a few reads getting by do not impact downstream results.

### Running the QC Script

As previously mentioned, this part is pretty simple. Samples can be read into the pipeline in one of two ways:

(1) A directory of sample-specific directories
(2) A directory of R1 and R2 sample fastqs

  If you go this latter way, you won't be able to concatenate across lanes too easily -- this is on my roadmap -- so ideally you have just one lane of data when it is stored this way.

First copy the code to a directory for the project so you have the current working version. Let's say we're working on a project called `myMetagenomics`. If you are on Wynton, you probably already have your sequence data uploaded, so you can start from there or do something like:

```
## Make a directory for your study
mkdir myMetagenomics
## Copy the current working metagneomics code from Katie's directory to yours
cp /wynton/group/lynch/kmccauley/Metagenomics_Pipeline_KM/ myMetageomics/
```

The script that you will work with first is `array_working_pipeline.sh`. Open it up in your favorite viewer (vim, emacs, nano, etc). The top should look something like this:

```
#!/usr/bin/bash
#$ -cwd
#$ -pe smp 10
#$ -l mem_free=10G
#$ -t 1:30
#$ -l h_rt=300:00:00
#$ -m e
#$ -R y

## Metagenomics Pipeline
## Current Version Written by Katie McCauley, 25 Nov 2019 with Ariane Panzer; last revised 02May2021

## Step 1:
## Check the above information to make sure you have the number of samples you plan to run (-t 1:30 is 30 samples), cores per sample (-pe smp 10 uses 10 cores), etc. Keep run time (h_rt) at 300 hours (it doesn't hurt and mostly helps); memory (mem_free) is also fine -- don't go too much higher

## Step 2:
# Where do your metagenomics fastq files live? This directory either contains (a) all FASTQ files for your metagenomics or (b) directories of sample data (often what is downloaded directly from BaseSpace)
## If (b), you can have other files in the directory of directories, but they cannot be FASTQ files

FASTA_DIRECTORY="/your/fasta/directory/location/here/"

## Where do you want your data to end up?
RESULT_DIRECTORY="/where/you/want/your/files/to/end/up/"

## Step 3:
## Now submit to the queue and cross your fingers (and maybe your toes)!
## qsub array_working_pipeline.sh
```

It's pretty well-documented here, but review the steps you see in the bash script and then make edits as necessary. Once you are happy, save your edits, exit and run `qsub array_working_pipeline.sh` from your `myMetagenomics` directory. That's it! You're done! Woot!

The script will figure out how your metagenomics data is stored and determine how to pull the necessary information so that your samples enter the queue as space becomes available.

As the pipeline is running, there will be a .o and .e file for every sample. If you view these files, you can get a sense of progress and general status. All of the verbosity I've built into the pipeline can be explored in the .o files, and the standard output of tools like bbmap can be viewed in the .e files. The .o files are very useful for troubleshooting.

This whole process can be long or short... Who knows... Depends on lots of stuff (proportion of human sequences, sequence complexity, read depth, etc), so check in on it occasionally and don't be surprised if you're waiting a couple of days to move on to the next step. However, setting up the code as an array *is* still a good thing -- imagine waiting for each of your samples to go through one at a time and one of those samples takes a week alone.

## Reads Based Metagenomics

Congratulations! You're back!

### Methods and Included Software

The steps in this portion of the pipeline are two-fold-fold: (1) two softwares, and (2) two sets of code!

The software:

1. Kraken2

Kraken2 assigns taxonomy to your reads using a k-mer based approach in which the software goes through each k-mer, assigns an NCBI-based taxonomy, and then looks at all of the calls across your forward and reverse reads and essentially takes a consensus vote. It will give you the taxonomy that provides the greatest degree of confidence, even if that means providing lower resolution (ie, at the order level instead of at the genus level). The database currently in use includes bacteria, fungi, viruses, human, and plants, and includes everything with complete genomes on NCBI.

2. HUMAnN 3.0

[HUMAnN 3.0](https://huttenhower.sph.harvard.edu/humann/) is a software designed to predict the conserved function (remember, for most of you this is DNA-based, not RNA-based, so it's measuring functional potential, not actual expression/activity) of your metagenomic data. It uses marker genes and is still...predictive, since it's still looking at a single read essentially in isolation. It will also assign taxonomy to your reads, though it tends to be much more conservative about the types of bacteria present in your samples. Sometimes it can pick up fungi/viruses, most often not.

### Running Reads Based Analysis

This component needs to be undertaken in two steps. The first is an array-based step, much like the code for Sample QC, since there's no sample dependency to carry out Kraken2 or the computationally-intensive HUMAnN, so the gains of running it as an array are worthwhile to capitalize on. However, the step afterwards combines all of the HUMAnN 3 results into an analysis file. The script that you will use for this first step is the `reads_based_data_generation.sh` script.

Running it is very similar to running the QC script. Copy the file to your working directory, and then open it up and change the paths in the first few lines of the script to reflect the location of your QC'd data. Also, make sure that the number of array jobs in the top are accurate to the number of samples you want to process. Other than that, you're good to save and run `qsub reads_based_data_generation.sh`.

The second script `humann_readout_standardization_withKEGG.sh` combines your individual-sample HUMAnN3 data into an analyzable dataset in copies per million. It provides data in the form of Enzyme Commission (EC) numbers as well as KEGG pathways. There's no need to array the data this time, so just make sure the file paths are accurate before running the qsub command again.

### Performing analysis

Before I really say much about performing analysis, I should stress that analyzing your metagenomics data, much like your 16S data, should come with a plan -- what are you hoping to learn from your data? You just spent $$$$ to generate it, so clearly there is a goal in mind. 

Once Kraken2 is complete, you can use a software called [Pavian](https://ccb.jhu.edu/software/pavian/) (see the "Try this instead" section) to visualize your data. I typically install it locally and then copy my kraken files locally before running Pavian. This will get you your table of counts across all samples and some pretty Sankey plots on a per-sample basis. I have yet to perform true statistical analyses with these datasets -- often I've already analyzed immense amounts of 16S rRNA data. Mostly it's to give me a sense of what types of bacteria are contained within the microbiome data, and perhaps perform a comparison (with Procrustes) to my 16S rRNA data. However, a few thoughts I have include:

1) A need to determine a method of standardization for the data -- whether that would be something like multiply-rarefying with our in-house script, or using a copies per million transformation.
2) Also determine at what level you will analyze the data at. Kraken2 provides the data at each level of phylogeny, and some reads are defined at the species level, while others are defined at the order level, and there's the trade off of information and specificity. I am working on testing out an alternative reads-based taxonomic classifier, Kaiju, but have mixed feelings about it so far.

Analysis of the HUMAnN3 data is very similar to analyzing 16S data for the most part. You can generate distance matrices and perform PERMANOVA analyses, or differential function analysis using the many model script. Or generate networks with WGCNA (which is better for the distribution of HUMAnN3 data than SparCC, which is better for 16S rRNA data). This is really all up to you. The one note that I would highly suggest is removing the rows that stratify out the pathways into different bacterial contributors. I often use this information to investigate a functional finding (ie, finding that a Tryptophan pathway is associated with asthma and then determining if there are specific bacteria that appear to be driving that association. 

## Assembly Pipeline

Assembly is a different beast at this point -- it's best use is as a follow-up to your reads-based analysis. What did you learn from that analysis that you think is worth diving into? The script itself is built to "perform assembly" but doesn't yet have a clear analytic post-assembly plan as far as my scripts are concerned. You won't get a nice analysis-ready object at the end, and that's probably going to feel a bit unsatisfying. You'll kind of need to create that for yourself at this point, but the goal of this is to get the major nuts and bolts of assembly out of the way so you can really dive into your question.

And before going further, I again want to stress that it's still very much in beta mode. As we have larger discussions, the script, and thus this documentation, will evolve.

P.S. Several of the steps I describe also work with metatranscriptomics data, but it's not currently available in the pipeline. Talk to me (Katie) if you need something for metatranscriptomics.

### Software Used

Probably the best way to introduce the pipeline is to just go through the software that I have pipelined.

metaSPAdes

metaSPAdes with default parameters and paired reads 

bbMap and basic filter

quast & DeepMAsED

Kraken2

Prokka

Binning

Metabat2

CONCOCT

Bin Checks: CheckM & DAStool

# Bug Reports and Roadmap

To submit a bug report, please use this form: https://airtable.com/shrhBdwkZJAQLJ77p

Roadmap:
- Allow for filtering of non-human hosts (mice, cat, dog).
- For one-dir-of-all-files setup, allow for multiple lanes.
- Add a second host-removal step based on bowtie2 or STAR, based on DNA vs or RNA.
- Find out why sometimes the reads-based pipeline works and sometimes doesn't.
- <s>Rename scripts so that they make more sense to an outside observer.</s>
