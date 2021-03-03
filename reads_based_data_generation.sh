#!/bin/bash
#$ -cwd
#$ -pe smp 5
#$ -l mem_free=10G
#$ -l h_rt=80:00:00
#$ -t 1-5
## #$ -m e
## #$ -M kathryn.mccauley@ucsf.edu

## Where is your QC results directory?
files_dir=/wynton/group/lynch/kmccauley/URECA_metagenomics/
## What is the name of the specific directory?
results_dir=metagenomics_results_113481



######## No edits needed beyond this point
software_location=/wynton/group/lynch/kmccauley/mySoftware/

cd $files_dir
sample1=$(echo `ls $results_dir` | cut -d " " -f $SGE_TASK_ID)

echo "Processed Sample is:" $sample1

## Jump straight into Kraken

db_location=/wynton/group/lynch/Shared/kraken2_db/
current_dir=`pwd`
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
    mkdir -p "$TMPDIR"
    export TMPDIR
fi

## run kraken providing the paired reads
cp $results_dir/$sample1/bbduk_ReadCleaning_Output/HumanFiltering/${sample1}_R1_clean_001.fastq.gz $TMPDIR
cp $results_dir/$sample1/bbduk_ReadCleaning_Output/HumanFiltering/${sample1}_R2_clean_001.fastq.gz $TMPDIR

cd $TMPDIR

## kraken runs, but it's not saving the actual output
singularity exec -B $PWD:$PWD,$db_location:/mnt ${software_location}kraken2.img kraken2 --db /mnt/ --report ${sample1}_kraken_report.txt --use-names --gzip-compressed --confidence 0.95 --memory-mapping --paired --threads $NSLOTS ${sample1}_R1_clean_001.fastq.gz ${sample1}_R2_clean_001.fastq.gz > ${sample1}_kraken_output.txt

## Run Humann3 using catted reads, so cat reads first:
cat ${sample1}_R1_clean_001.fastq.gz ${sample1}_R2_clean_001.fastq.gz > ${sample1}_clean_001.fastq.gz

## then run humann3
singularity exec -B $PWD:$PWD,/wynton/group/lynch/Shared/humann_db/metaphlan:/metaphlan,/wynton/group/lynch/Shared/humann_db/chocophlan:/chocophlan,/wynton/group/lynch/Shared/humann_db/uniref:/uniref $software_location/humann.img humann -i ${sample1}_clean_001.fastq.gz -o humann3_results --metaphlan-options "--bowtie2db /metaphlan" --nucleotide-database /chocophlan --protein-database /uniref --threads $NSLOTS

rm *.gz
rm -r humann3_results/*_temp

mkdir -p $current_dir/humann3_results/${sample1}
cp -r $TMPDIR/* $current_dir/humann3_results/${sample1}

