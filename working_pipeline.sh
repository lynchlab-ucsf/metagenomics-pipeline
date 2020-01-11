#!/usr/bin/bash
#$ -cwd
#$ -pe smp 15
#$ -l mem_free=20G
#$ -l h_rt=48:00:00
#$ -l s_rt=48:00:00
#$ -m bae
#$ -M kathryn.mccauley@ucsf.edu

## Metagenomics Pipeline
## Current Version Written by Katie McCauley, 25 Nov 2019 with Ariane Panzer; last revised 03 Dec 2019

# Where do your metagenomics fastq files live?
FASTA_DIRECTORY="/wynton/group/lynch/kmccauley/01_raw/"

######## All variables beyond this point do not need to be set.
software_location="/wynton/group/lynch/kmccauley/mySoftware" ## Making a variable now in the event that this location changes.

## Define Temporary Directory location (in scratch directory). Not using node-specific scratch space
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

FASTQC_OUT_1="FastQC_Out_1"
FASTQC_RAW_1="$TMPDIR"

echo "Temporary directory is" $TMPDIR ", which will be deleted upon completion of the job"

CURRENT_DIR=`pwd`

echo "Final files will be saved to" $CURRENT_DIR

cd $TMPDIR

# Copy Metagenomics files to Scratch space.
cp -r ${FASTA_DIRECTORY}/* .

# Ensure the output directory always exists (setting up before loading CBI and the fastqc software)
mkdir -p "${FASTQC_OUT_1}"

module load CBI fastqc
## Will print information about the version of fastqc being used.
module list

## Here, I'll want to concatenate any samples spread across lanes (and maybe add a check to see if this is even needed so the time isn't wasted. Also want to delete the initial files so that we don't process them in addition to the concatenated files).
sampnames=`ls | grep '[.]fastq' | awk -F '_' '{print $1"_"$2}' | uniq`
for i in $sampnames;
  do for j in 1 2; 
     do
     cat ${i}_*_R${j}_001.fastq.gz > ${i}_R${j}_001_fastq.gz;
     rm ${i}_*_R${j}_001.fastq.gz ;
  done; 
done;

## Instead of running files ending in fastq, some files will also end in fastq.gz. Therefore, I am creating a variable called `files` that identifies any file in the FASTA_DRECTORY that includes ".fastq", which will include both ".fastq" files and ".fastq.gz" files.
files=`ls | grep "[.]fastq"` # Gets any file in the directory with fastq in the name.

## I'm trying really hard to get this right... something to come back to.
#resources=$(qstat -j "$JOB_ID" | grep -F "hard resource_list:" | sed 's/.*: *//')
#IFS=, read obj1 obj2 obj3 <<< $resources
#value='mem_free'
#for i in $obj1 $obj2 $obj3; do
#	if [[ $i =~ " ${value} " ]]; then
#		mem1=$(echo $i | cut -d= -f2)
#	fi
#done
#num1="${resources//[!0-9]/}" ## This will pull out the number from my mem_free string...

## Now that we have the file names, determine where overlaps in *lanes* exist within samples, which may be easier to do if I'm able to pull in a .txt file....
## I'm still playing around with how I want to do this... My brain isn't working at its best right now....

## Determine how to cat files that are from different lanes.


run_fastqc () {
        for f in $files; do
        echo "Running FastQC on" $f # Print file name to remaine apprised of progress
                fastqc \
                -t $fastqc_threads \
                "${f}" \
                --nogroup \
                --outdir "${FASTQC_OUT_1}" \
                --dir "${TMPDIR}"
        echo "Finished running FastQC on" $f
        done
}

# Definiting Output Directory names and creating said directories
BBDUK_DIR="bbduk_ReadCleaning_Output"
TRIMMED_DIR="AdapterTrimming"
HUMAN_DIR="HumanFiltering"
PhiX_DIR="PhiXRemoval"
QualFilt_DIR="QualitySeqs"

mkdir -p "${BBDUK_DIR}"/"${PhiX_DIR}"
mkdir -p "${BBDUK_DIR}"/"${QualFilt_DIR}"
mkdir -p "${BBDUK_DIR}"/"${HUMAN_DIR}"
mkdir -p "${BBDUK_DIR}"/"${TRIMMED_DIR}"

MEMORY="-Xmx40g"
PARAMETERS_ADAPT_TRIM="ktrim=r k=23 mink=11 hdist=1 tbo tpe" # can also include tbo (trim adapters based on pair overlap detection using BBMerge [which does not require known adapter sequnces]) or tpe (specifies to trim both reads to the same length in the event that an adapter kmer was only detected in one of them)

## Current function uses the $files list generated above to ID only the Read1 files, confirms there's a Read2 file, and then continues.

run_bbtools() {
for f in $files; do
if [[ $f == *"_R1_"* ]] && test -f "${f/_R1/_R2}"; then
echo "Running BBDuk on" $f ";" "${f/_R1/_R2}" "exists."

echo "Adapter trimming for" $f
## Adapter Trimming
singularity exec ${software_location}/bbtools.img bbduk.sh \
"${MEMORY}" \
threads=$bbduk_threads \
in1=$f \
in2="${f/_R1/_R2}" \
out1="${BBDUK_DIR}"/"${TRIMMED_DIR}"/"${f/_R1/_R1_adapt_trim}" \
out2="${BBDUK_DIR}"/"${TRIMMED_DIR}"/"${f/_R1/_R2_adapt_trim}" \
ref=/bbmap/resources/adapters.fa \
"${PARAMETERS_ADAPT_TRIM}"

echo "PhiX removal for" $f
## PhiX Removal
singularity exec ${software_location}/bbtools.img bbduk.sh \
"${MEMORY}" \
threads=$bbduk_threads \
in1="${BBDUK_DIR}"/"${TRIMMED_DIR}"/"${f/_R1/_R1_adapt_trim}" \
in2="${BBDUK_DIR}"/"${TRIMMED_DIR}"/"${f/_R1/_R2_adapt_trim}" \
out1="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R1_remPhiX}" \
out2="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R2_remPhiX}" \
ref=/bbmap/resources/phix174_ill.ref.fa.gz \
k=31 hdist=1 stats="${BBDUK_DIR}"/"${PhiX_DIR}"/phix_stats${f}.txt

echo "Quality Filtering and Final Quality Histograms for" $f
## Quality Filtering/Trimming
singularity exec ${software_location}/bbtools.img bbduk.sh \
"${MEMORY}" \
threads=$bbduk_threads \
in1="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R1_remPhiX}" \
in2="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R2_remPhiX}" \
out1="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R1_qTrim}" \
out2="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R2_qTrim}" \
k=23 hdist=1 qtrim=rl trimq=15 maq=10 \
qhist="${BBDUK_DIR}"/"${QualFilt_DIR}"/qhist_${f}.txt

echo "Removing Human Sequences for" $f
singularity exec -B "${software_location}:/mnt" ${software_location}/bbtools.img bbmap.sh \
"${MEMORY}" \
threads=$bbduk_threads \
in="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R1_qTrim}" \
in2="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R2_qTrim}" \
outm1="${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_hg38map}" \
outu1="${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_clean}" \
outm2="${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R2_hg38map}" \
outu2="${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R2_clean}" \
path=/mnt

fi
done
}

## Move human reads off of Wynton and to Lynchserver2 to pipe into my pipeline for human gender and SNP calling.
## Try to think about how to implement this when I don't want the person to need to interact with the pipeline in any way. Could I write a second script that identifies a directory that the human reads go, and then copies any files occupying that directory over to lynchserver2 maybe using the manager account? Still don't want to make it necessary to enter a password.... Maybe discuss this at lab meeting??
## scp "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_hg38map}" "${lynchserver2_user}"@lynchserver2.ucsf.edu:/data/

## Looking at Ariane's Parameters to determine what she used for each step.
#PARAMETERS_ADAPT_TRIM="ktrim=r k=23 hdist=1 mink=11" 
#PARAMETERS_PHIX="k=31 hdist=1 stats=stats.txt"
#PARAMETERS_QUALITY_TRIM="k=27 hdist=1 qtrim=rl trimq=17 cardinality=t"
## Initial comments are that I am unsure why different kmer lengths were used for each step. Some discussion about ideal parameters for each of these steps may be warranted.

if [[ $NSLOTS -gt 3 ]]; then

# Determines distribution of threads for initial processes, based on number of cores alloted for the submission:
	let cores=$NSLOTS
	let fastqc_threads=cores/3
	let bbduk_threads=cores-fastqc_threads

	echo $NSLOTS "cores detected; running FastQC on" `expr $NSLOTS - $bbduk_threads` "core(s) & quality processing on" ${bbduk_threads} "core(s)"
        run_fastqc &
        run_bbtools &
        wait ${!} ## This will wait for the two functions above to complete before continuing on.
else
	let fastqc_threads=$NSLOTS
	let bbduk_threads=$NSLOTS
        echo $NSLOTS "core(s) is less than or equal to 3. Will not run FastQC and quality processing in parallel, but will use the specified number of cores for each process. If you wish to run these processes in parallel, change the -pe smp option to a number greater than 3."
        run_fastqc
        run_bbtools
fi 
#############################################################################
### FUNCTION TO RUN MIDAS     #####
#############################################################################


## Set up Code for MIDAS, which I think I'm going to try to parallelize in the same way like I did for QC and BBtools
## The question, though, is whether activating the MIDAS virtual environment in some way impacts the metaSPAdes run... I would like to just send the command to the virtual environment, instead of opening the virtual environment to support the running of my commands....

run_midas() {
for f in $files; do

if [[ $f == *"_R1_"* ]] && test -f "${f/_R1/_R2}"; then
echo "Running MIDAS" on "$f"

source ${software_location}/metagenomics_midas2/bin/activate

run_midas.py species ./MIDAS_results/"${f/001.fastq.gz/_MIDAS}" -1 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_clean}" -2 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R2_clean}" -t $midas_cores
run_midas.py genes ./MIDAS_results/"${f/001.fastq.gz/_MIDAS}" -1 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_clean}" -2 "${BBDUK_DIR}"/"${HUMAN_DIR}"/${f/_R1/_R2_clean} -t $midas_cores
run_midas.py snps ./MIDAS_results/"${f/001.fastq.gz/_MIDAS}" -1 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_clean}" -2 "${BBDUK_DIR}"/"${HUMAN_DIR}"/${f/_R1/_R2_clean} -t $midas_cores

deactivate
fi
done
}

###############################################################################
###### FUNCTION TO RUN METASPADES        ##########
###############################################################################



merge_and_make_contigs() {
for f in $files; do
if [[ $f == *"_R1_"* ]] && test -f "${f/_R1/_R2}"; then

echo "Merging Paired Reads for" $f

singularity exec /wynton/group/lynch/kmccauley/mySoftware/vsearch_bc.img vsearch \
--fastq_mergepairs "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_clean}" \
--reverse "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R2_clean}" \
--eetabbedout ${f/.fastq.gz/_}stats.txt \
--fastq_allowmergestagger \
--fastqout "${BBDUK_DIR}"/"${HUMAN_DIR}"/merged_"${f}" \
--threads $spades_cores1


echo "Running metaSPAdes on" $f "with merged reads using" $spades_cores "core(s)."
metaspades.py -k 21,33,55,77 \
--pe1-1 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_clean}" \
--pe1-2 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R2_clean}" \
--pe1-12 "${BBDUK_DIR}"/"${HUMAN_DIR}"/merged_"${f}" \
-o metaspades_results/"${f/.fastq.gz/_contig_dat_withmerge}" \
-t $spades_cores1

fi
done
}

just_make_contigs() {
for f in $files; do
if [[ $f == *"_R1_"* ]] && test -f "${f/_R1/_R2}"; then

echo "Running metaSPAdes on" $f "with merged reads using" $spades_cores "core(s)."
metaspades.py -k 21,33,55,77 \
--pe1-1 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_clean}" \
--pe1-2 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R2_clean}" \
-o metaspades_results/"${f/.fastq.gz/_contig_dat_nomerge}" \
-t $spades_cores2

fi
done
}



############################################
### RUN MIDAS AND METASPADES IN PARALLEL ###     ###Unsure if this is safe to do.....
############################################


## Determine how we want to download our own database.

# This setting of the PATH variable works:
export PATH=$PATH:${software_location}/SPAdes-3.11.1-Linux/bin/

## try spades on the one sample:
mkdir metaspades_results
mkdir MIDAS_results


## After a test run, it may be best to allot more cores to metaSPAdes -- it takes *WAY* longer to run than MIDAS. This is assuming that it is safe to run the virtual environment in parallel with metaSPAdes....

let midas_cores=cores
#let spades_cores1=(cores-midas_cores)/2
#let spades_cores2=cores-midas_cores-spades_cores1



run_midas &
#merge_and_make_contigs &
#just_make_contigs &
wait

echo "PIPELINE COMPLETE! Moving files back to" $CURRENT_DIR

rm $files 

## Adding a job ID and date to the name of the directory that gets saved.
currdate=`date +%m%d%Y`
mkdir "$CURRENT_DIR"/metagenomics_results_${currdate}_${JOB_ID}/
mv ./* "$CURRENT_DIR"/metagenomics_results_${currdate}_${JOB_ID}/
