# Metagenomics Pipeline

Author: Katie McCauley
Date: May 19, 2021

# Table of Contents

## Intro

This document will describe how the metagenomics pipeline is run and give a sense of what is contained within the pipeline itself. This document will encompass the information I provide in the live workshops I've given to date. It will also be a live document which will change as the pipeline changes (the pipeline is not static).

The pipeline has been built to maximize efficiency as much as possible. This means that you may need to run the pipeline once, review your results, and run the pipeline again. However, I hope that despite the length of this document, you will notice that the actual effort required to do so is quite minimal. It was developed so that minimal, if any, hands-on manipulation or knowledge of bash scripting would be needed on the user side. There are likely to be outlier instances (ie, really low-quality sequencing run) where we need to start thinking about modifying some parameters, but unless you have a known outlier case, no major edits should be needed.

This pipeline has also been validated on RNAseq data. It uses bbmap which is a splice-aware aligner, so both DNA and RNA can be aligned to the human genome and filtered out. Some discussion exists that suggests that bbmap is imperfect, but alas... A future iteration of the pipeline will likely include a second human filtering step that will use bowtie2 for DNA sequences and STAR for RNA. This will likely be another parameter to provide to the script, but we'll get there when we get there.

## QC Steps

The steps in this section of the pipeline perform the following:

1. Concatenate reads across lanes
2. Run FastQC
3. Use bbTools for:
  a. Adapter Trimming
  b. PhiX Removal
  c. Quality Filtering
  d. Human Removal
  
### How do I actually run this?
  
## Reads-Based Analysis

### And running this step?
