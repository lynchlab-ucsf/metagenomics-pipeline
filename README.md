# Metagenomics Pipeline

Author: Katie McCauley

Date: May 19, 2021

# Table of Contents

- [Introduction](#introduction)
  - [Quality Control Steps](#quality-control-steps)
  - [Running the QC Script](#running-the-qc-script)

# Introduction
This document will describe how the metagenomics pipeline is run and give a sense of what is contained within the pipeline itself. This document will encompass the information I provide in the live workshops I've given to date. It will also be a live document which will change as the pipeline changes (the pipeline is not static).

The pipeline has been built to maximize efficiency as much as possible. This means that you may need to run the pipeline once, review your results, and run the pipeline again. However, I hope that despite the length of this document, you will notice that the actual effort required to do so is quite minimal. It was developed so that minimal, if any, hands-on manipulation or knowledge of bash scripting would be needed on the user side. There are likely to be outlier instances (ie, really low-quality sequencing run) where we need to start thinking about modifying some parameters, but unless you have a known outlier case, no major edits should be needed.

This pipeline has also been validated on RNAseq data. It uses bbmap which is a splice-aware aligner, so both DNA and RNA can be aligned to the human genome and filtered out. Some discussion exists that suggests that bbmap is imperfect, but alas... A future iteration of the pipeline will likely include a second human filtering step that will use bowtie2 for DNA sequences and STAR for RNA. This will likely be another parameter to provide to the script, but we'll get there when we get there.

## Quality Control Steps
The steps in this section of the pipeline perform the following:

1. Concatenate reads across lanes

2. Run FastQC

3. Use bbTools for:

{add additional information about the outputs, etc.}

## Running the QC Script

As mentioned above, this part is pretty simple. First copy the code to a directory for the project so you have the current working version. Let's say we're working on a project called `myMetagenomics`. If you are on Wynton, you probably already have your sequence data uploaded, so you can start from there or do something like:

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

This whole process can be long or short... Who knows... Depends on lots of stuff (proportion of human sequences, read depth, etc), so check in on it occasionally and don't be surprised if you're waiting a couple of days to move on to the next step.
