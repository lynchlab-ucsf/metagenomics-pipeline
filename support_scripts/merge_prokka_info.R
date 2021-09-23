## Script to combine the PROKKA results into one file for one sample. Paste didn't cut it for what I needed
library(dplyr)

argv <- commandArgs(TRUE)

file1 <- read.table(argv[1], sep="\t", header=TRUE, quote="")
file2 <- read.table(argv[2], sep="\t", header=TRUE, quote="")
file3 <- read.table(argv[3], sep="\t", header=TRUE, quote="")

combined <- file1 %>%
   full_join(file2, by=c("locus_tag"="SubContigID")) %>%
   full_join(file3, by=c("locus_tag"="Gene")) %>%
   mutate(SampleID = argv[4])

names(combined)[names(combined) %in% "Counts"] <- "ReadCounts"
names(combined)[length(combined)] <- "Coverage"

write.table(combined, paste0(argv[4], "_sample_level_summary.tsv"), sep="\t", quote=F, row.names=F)

## Then, once all the samples are done running, we can bind the columns together and do comparative analyses
