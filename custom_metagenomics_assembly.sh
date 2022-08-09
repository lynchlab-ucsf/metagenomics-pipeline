#!/bin/bash
#$ -cwd
#$ -pe smp 10
#$ -t 1:30
#$ -l mem_free=10G
#$ -l h_rt=300:00:00
#$ -R y

files_dir=/wherever/your/metagenomics_results/file/lives
results_dir=metagenomics_results_XXXXX

#########################
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
db_dir=/wynton/group/lynch/Shared/
scripts_loc=/wynton/group/lynch/kmccauley/Metagenomics_Pipeline_KM/support_scripts/
cd $TMPDIR

$spades_loc/spades.py -o meta_spades --meta -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz -t $NSLOTS

awk '/^>/{print ">contig" ++i; next}{print}' < meta_spades/contigs.fasta > meta_spades/contigs_renamed.fasta

singularity exec -B $PWD:$PWD $software_dir/bbtools.img bbmap.sh in=${sample}_R1.fastq.gz in2=${sample}_R2.fastq.gz ref=rna_spades/contigs_renamed.fasta covstats=coveragestats.txt po=t outm=${sample}_results.bam

module load CBI r
Rscript $files_dir/bbduk_contig_filter.R

singularity exec -B $PWD:$PWD $software_dir/qiime.img filter_fasta.py -f rna_spades/contigs_renamed.fasta -o rna_spades/filtered_contigs.fasta -s contigs_to_keep.txt

singularity exec -B $PWD:$PWD $software_dir/quast.sif quast.py --threads $NSLOTS rna_spades/filtered_contigs.fasta --glimmer
## Try out DeepMAsED from Ruth Ley's Lab as well -- perhaps switch to regular quast above as well?:
mkdir -p deepmased
echo "bam	fasta" > deepmased/bam_fasta_file.tsv
echo "${sample}_results.bam	rna_spades/contigs_renamed.fasta" >> deepmased/bam_fasta_file.tsv
singularity exec -B $PWD:$PWD $software_dir/deepmased.img DeepMAsED features deepmased/bam_fasta_file.tsv --outdir deepmased/ -p $NSLOTS
singularity exec -B $PWD:$PWD $software_dir/deepmased.img DeepMAsED predict deepmased/feature_file_table.tsv --save-path deepmased/ --cpu-only --n-procs $NSLOTS

module load samtools/1.10 bedtools2

samtools sort -o ${sample}_sorted.bam -@ $NSLOTS ${sample}_results.bam
samtools index -@ $NSLOTS ${sample}_sorted.bam

##
## Add kraken2

singularity exec -B $PWD:$PWD,$db_dir/kraken2_complete/:/mnt $software_dir/kraken2.img kraken2 --db /mnt/ --report ${sample}_kraken_report_0.2.txt --use-names --confidence 0.2 --memory-mapping --threads $NSLOTS rna_spades/filtered_contigs.fasta --output ${sample}_kraken_output_0.2.txt --unclassified-out ${sample}_unclassified_0.2.fa

## Run PROKKA (can prokka be used for metatranscriptomics?)
singularity exec -B $PWD:$PWD $software_dir/prokka.img prokka --outdir ${sample}_prokka --metagenome --cpus $NSLOTS rna_spades/filtered_contigs.fasta

## Organize PROKKA files
## basically I had created a file with the counts, coverage, gene IDs and other stuff as a summary that I was going to bring together in a second script.
bash $scripts_loc/gfftobed.sh ${sample}_prokka/PROKKA_*.gff > ${sample}_prokka/bedfile.bed
bedtools coverage -hist -b ${sample}_results.bam -a ${sample}_prokka/bedfile.bed > ${sample}_prokka/${sample}.map.hist
python3 $scripts_loc/get_coverage_for_genes.py -i ${sample}_prokka/${sample}.map.hist > ${sample}_prokka/${sample}.coverage

bedtools coverage -counts -b ${sample}_results.bam -a ${sample}_prokka/bedfile.bed > ${sample}_prokka/${sample}.counts
echo -e "ContigID\tTranscripts" > ${sample}_prokka/${sample}_summary.tsv && cut -f1,5 ${sample}_prokka/${sample}.counts >> ${sample}_prokka/${sample}_summary.tsv

paste -d'\t' ${sample}_prokka/PROKKA_*.tsv ${sample}_prokka/${sample}_summary.tsv > test1.tsv
paste -d'\t' test1.tsv ${sample}_prokka/${sample}.coverage > ${sample}_prokka/${sample}_sample_level_summary.tsv
rm test1.tsv

## Make a workable R table:
Rscript $scripts_loc/merge_prokka_info.R ${sample}_prokka/PROKKA_*.tsv ${sample}_prokka/${sample}_summary.tsv ${sample}_prokka/${sample}.coverage $sample

mkdir binning

## Build MAGs using metabat and CONCOCT
singularity exec -B $PWD:$PWD $software_dir/metabat2.img metabat2 -i meta_spades/filtered_contigs.fasta -o binning/metabat_bins/${sample} -m 1500 --seed 123 --saveCls -t $NSLOTS -v -s 10000

singularity exec -B $PWD:$PWD $software_dir/concoct.img cut_up_fasta.py meta_spades/filtered_contigs.fasta -c 10000 -o 0 --merge_last -b binning/contigs_10K.bed > binning/contigs_10K.fa
singularity exec -B $PWD:$PWD $software_dir/concoct.img concoct_coverage_table.py binning/contigs_10K.bed ${sample}_sorted.bam > binning/coverage_table.tsv
singularity exec -B $PWD:$PWD $software_dir/concoct.img concoct --composition_file binning/contigs_10K.fa --coverage_file binning/coverage_table.tsv -b binning/concoct_output_${sample}/ --threads $NSLOTS
## Need an "if bins were created" argument
if ! grep -q 'Not enough contigs pass the threshold filter' "binning/concoct_output_${sample}/log.txt"; then
singularity exec -B $PWD:$PWD $software_dir/concoct.img merge_cutup_clustering.py binning/concoct_output_${sample}/clustering_gt1000.csv > binning/concoct_output_${sample}/clustering_merged.csv
mkdir binning/concoct_bins/
singularity exec -B $PWD:$PWD $software_dir/concoct.img extract_fasta_bins.py meta_spades/filtered_contigs.fasta binning/concoct_output_${sample}/clustering_merged.csv --output_path binning/concoct_bins/
else
echo "Not enough contigs pass the threshold filter"
fi

### check bins with checkM
singularity exec -B $PWD:$PWD $software_dir/checkm.img checkm lineage_wf binning/metabat_bins metabat_checkm -x .fa -t $NSLOTS
singularity exec -B $PWD:$PWD $software_dir/checkm.img checkm lineage_wf binning/concoct_bins concoct_checkm -x .fa -t $NSLOTS

### also check/improve bins with DAS Tool
singularity exec -B $PWD:$PWD $software_dir/DASTool_v1_1_4.img Fasta_to_Contig2Bin.sh -e fa -i binning/concoct_bins/ > binning/concoct_contig2bin_${sample}.txt
singularity exec -B $PWD:$PWD $software_dir/DASTool_v1_1_4.img Fasta_to_Contig2Bin.sh -e fa -i binning/metabat_bins/ > binning/metabat_contig2bin_${sample}.txt
singularity exec -B $PWD:$PWD $software_dir/DASTool_v1_1_4.img DAS_Tool -i binning/concoct_contig2bin_${sample}.txt,binning/metabat_contig2bin_${sample}.txt -c meta_spades/filtered_contigs.fasta -o binning/dastool_results

#Remove the initial fastq files
rm *.fastq.gz

## Creating directory where results from this will get saved
mkdir -p $files_dir/results_assembly_${JOB_ID}/$sample
mv ./* $files_dir/results_assembly_${JOB_ID}/$sample

