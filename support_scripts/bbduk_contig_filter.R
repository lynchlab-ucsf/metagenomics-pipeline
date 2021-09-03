#!/usr/bin/Rscript

contig.info <- read.table("coveragestats.txt", header=T, check.names=F, sep="\t", comment="") ## Need to come up with more generic name for this file.

subset.contigs <- subset(contig.info, Avg_fold > 5 & Length > 150) ## Reducing to 5x cov and 150bp length for read-recruitment

keep.contigs <- c(subset.contigs[,1])
write.table(keep.contigs, "contigs_to_keep.txt", sep="\t", row.names=F, quote=F, col.names=F)
