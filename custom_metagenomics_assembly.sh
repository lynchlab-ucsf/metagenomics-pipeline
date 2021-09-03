#!/bin/bash
#$ -cwd
#$ -pe smp 10
#$ -t 1:3
#$ -l mem_free=10G
#$ -l h_rt=300:00:00
#$ -R y

files_dir=/wynton/group/lynch/kmccauley/URECA_assembly_testing/
results_dir=raw_files


cd $files_dir
sample=$(echo `ls $results_dir` | cut -d " " -f $SGE_TASK_ID)

echo "Processed Sample is:" $sample

if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
    mkdir -p "$TMPDIR"
    export TMPDIR
fi

echo "Temporary directory is" $TMPDIR ", which will be deleted upon completion of the job"

cp $files_dir/$results_dir/$sample/bbduk_ReadCleaning_Output/HumanFiltering/*R1*clean* $TMPDIR/${sample}_R1.fastq.gz
cp $files_dir/$results_dir/$sample/bbduk_ReadCleaning_Output/HumanFiltering/*R2*clean* $TMPDIR/${sample}_R2.fastq.gz

spades_loc=/wynton/group/lynch/kmccauley/mySoftware/SPAdes-3.11.1-Linux/bin/
software_dir=/wynton/group/lynch/kmccauley/mySoftware/
scripts_loc=/wynton/group/lynch/kmccauley/Metagenomics_Pipeline_KM/support_scripts/
db_dir=/wynton/group/lynch/Shared/
cd $TMPDIR

$spades_loc/spades.py -o meta_spades --meta -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz -t $NSLOTS

awk '/^>/{print ">contig" ++i; next}{print}' < meta_spades/contigs.fasta > meta_spades/contigs_renamed.fasta

singularity exec -B $PWD:$PWD $software_dir/bbtools.img bbmap.sh in=${sample}_R1.fastq.gz in2=${sample}_R2.fastq.gz ref=meta_spades/contigs_renamed.fasta covstats=coveragestats.txt po=t outm=${sample}_results.bam

module load CBI r
Rscript $scripts_loc/bbduk_contig_filter.R

singularity exec -B $PWD:$PWD $software_dir/qiime.img filter_fasta.py -f meta_spades/contigs_renamed.fasta -o meta_spades/filtered_contigs.fasta -s contigs_to_keep.txt

singularity exec -B $PWD:$PWD $software_dir/quast.sif metaquast.py --threads $NSLOTS meta_spades/filtered_contigs.fasta --glimmer

module load samtools bedtools2

samtools sort -o ${sample}_sorted.bam -@ $NSLOTS ${sample}_results.bam
samtools index -@ $NSLOTS ${sample}_sorted.bam

##
## Add kraken2

singularity exec -B $PWD:$PWD,$db_dir/kraken2_complete/:/mnt $software_dir/kraken2.img kraken2 --db /mnt/ --report ${sample}_kraken_report_0.2.txt --use-names --confidence 0.2 --memory-mapping --threads $NSLOTS meta_spades/filtered_contigs.fasta --output ${sample}_kraken_output_0.2.txt --unclassified-out ${sample}_unclassified_0.2.fa

## Run PROKKA (can prokka be used for metatranscriptomics?)
singularity exec -B $PWD:$PWD $software_dir/prokka.img prokka --outdir ${sample}_prokka --metagenome --cpus $NSLOTS meta_spades/filtered_contigs.fasta

## Organize PROKKA files
## basically I had created a file with the counts, coverage, gene IDs and other stuff as a summary that I was going to bring together in a second script.
bash $files_dir/gfftobed.sh ${sample}_prokka/PROKKA_*.gff > ${sample}_prokka/bedfile.bed
bedtools coverage -hist -b ${sample}_results.bam -a ${sample}_prokka/bedfile.bed > ${sample}_prokka/${sample}.map.hist
python3 $files_dir/get_coverage_for_genes.py -i ${sample}_prokka/${sample}.map.hist > ${sample}_prokka/${sample}.coverage

bedtools coverage -counts -b ${sample}_results.bam -a ${sample}_prokka/bedfile.bed > ${sample}_prokka/${sample}.counts
echo -e "ContigID\tTranscripts" > ${sample}_prokka/${sample}_summary.tsv && cut -f1,5 ${sample}_prokka/${sample}.counts >> ${sample}_prokka/${sample}_summary.tsv

paste -d'\t' ${sample}_prokka/PROKKA_*.tsv ${sample}_prokka/${sample}_summary.tsv > test1.tsv
paste -d'\t' test1.tsv ${sample}_prokka/${sample}.coverage > ${sample}_prokka/${sample}_sample_level_summary.tsv
rm test1.tsv

## Make a workable R table:
Rscript -e "test <- read.table('${sample}_prokka/${sample}_sample_level_summary.tsv', sep='\t', nrow=3, header=T); write.table(test,'${sample}_prokka/${sample}_sample_level_summary_clean.tsv')"

## Build MAGs using metabat and CONCOCT
singularity exec -B $PWD:$PWD $software_dir/metabat2.img metabat2 -i meta_spades/filtered_contigs.fasta -o bins/${sample}_metabat -m 1500 --seed 123

singularity exec -B $PWD:$PWD $software_dir/concoct.img cut_up_fasta.py meta_spades/filtered_contigs.fasta -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
singularity exec -B $PWD:$PWD $software_dir/concoct.img concoct_coverage_table.py contigs_10K.bed ${sample}_sorted.bam > coverage_table.tsv
singularity exec -B $PWD:$PWD $software_dir/concoct.img concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output_${sample}/ --threads $NSLOTS
## Need an "if bins were created" argument
if ! grep -q 'Not enough contigs pass the threshold filter' "concoct_output_${sample}/log.txt"; then
singularity exec -B $PWD:$PWD $software_dir/concoct.img merge_cutup_clustering.py concoct_output_${sample}/clustering_gt1000.csv > concoct_output_${sample}/clustering_merged.csv
singularity exec -B $PWD:$PWD $software_dir/concoct.img extract_fasta_bins.py meta_spades/filtered_contigs.fasta concoct_output_${sample}/clustering_merged.csv --output_path bins/
else
echo "Not enough contigs pass the threshold filter"
fi

### check bins with checkM
singularity exec -B $PWD:$PWD $software_dir/checkm.img checkm lineage_wf bins/ checkm/ ## Not sure if the no-bin error is fatal to the point where my files don't get returned.

#Remove the initial fastq files
rm *.fastq.gz

## Creating directory where results from this will get saved
mkdir -p $files_dir/results_assembly_${JOB_ID}/$sample
mv ./* $files_dir/results_assembly_${JOB_ID}/$sample


## This is separate, getting functions hopefully by sample?e
## $anvio anvi-mcg-classifier -p SAMPLES-MERGED/PROFILE.db -c all_contigs.db -O mcg_classifier_res

## Post-run steps (run from command line):
##Rscript  make_internal_genomes_table.R results_assembly_${JOB_ID}
## anvio="singularity exec -B $PWD:$PWD,/wynton/group/lynch/Shared/:/mnt/ /wynton/group/lynch/kmccauley/mySoftware/anvio.sif"
## $anvio anvi-gen-genomes-storage -i internal_genomes.txt -o libsamples-GENOMES.db
## What actually worked:: $anvio anvi-gen-genomes-storage -e internal_genomes.txt -o libsamples-GENOMES.db
## $anvio anvi-pan-genome -g libsamples-GENOMES.db --project-name "COVID" --output-dir single_assembly/ --num-threads 15 --mcl-inflation 10 --use-ncbi-blast

## $anvio anvi-display-pan -p single_assembly/COVID-PAN.db -g libsamples-GENOMES.db -P 8010

## $anvio anvi-compute-functional-enrichment -p single_assembly/COVID-PAN.db -g libsamples-GENOMES.db --category status --annotation-source COG20_FUNCTION -o enriched-functions-status.txt --functional-occurrence-table-output COVID-functions-occurrence-frequency.txt
## $anvio anvi-export-misc-data -p single_assembly/COVID-PAN.db -info.txt
## I'd like to be able to pull out the COG functions as a txt file, but not sure how to do that yet. Might need synteny? Or just the co-assembly....

