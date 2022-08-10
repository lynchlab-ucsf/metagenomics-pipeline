#!/usr/bin/bash
#$ -cwd
#$ -pe smp 10
#$ -l mem_free=10G
#$ -t 1:262
#$ -l h_rt=300:00:00
#$ -R y

## Metagenomics Pipeline
## Current Version Written by Katie McCauley, 25 Nov 2019 with Ariane Panzer; last revised 02May2021

## Step 1:
## Check the above information to make sure you have the number of samples you plan to run (-t 1:30 is 30 samples), cores per sample (-pe smp 10 uses 10 cores), etc. Keep run time (h_rt) at 300 hours (it doesn't hurt and mostly helps); memory (mem_free) is also fine -- don't go too much higher

## Step 2:
# Where do your metagenomics fastq files live? This directory either contains (a) all FASTQ files for your metagenomics or (b) directories of sample data (often what is downloaded directly from BaseSpace)
## If (b), you can have other files in the directory of directories, but they cannot be FASTQ files

FASTA_DIRECTORY="/wynton/group/lynch/kmccauley/ITN_Metagenomics/combined_files/"

## Where do you want your data to end up?
RESULT_DIRECTORY="/wynton/group/lynch/kmccauley/ITN_Metagenomics/"

## Step 3:
## Now submit to the queue and cross your fingers (and maybe your toes)!
## qsub sample_qc_pipeline.sh


######## All variables beyond this point do not need to be set.
software_location="/wynton/group/lynch/kmccauley/mySoftware" ## Making a variable now in the event that this location changes.

cd $FASTA_DIRECTORY

## This part determines how the data is stored and, later, processed.
if ls *fastq* &>/dev/null; then

FASTA_LIST=`ls *fastq* | awk 'NR % 2 == 1 {print}'` ## This obtains every other fastq file line in the dataset and gets us what we need.
FASTA_FILE=$(echo $FASTA_LIST | cut -d " " -f $SGE_TASK_ID)
echo $FASTA_FILE;

else

FASTA_LIST=`ls $FASTA_DIRECTORY`
FASTA_FILE_DIR=$(echo $FASTA_LIST | cut -d " " -f $SGE_TASK_ID)
echo "Fasta File is:"
echo $FASTA_FILE_DIR
cd $FASTA_FILE_DIR
FASTA_FILE=`ls *fastq* | awk 'NR==1 {print}'`
echo $FASTA_FILE;
fi

## Define Temporary Directory location (in scratch directory). Uses node-specific scratch space.
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
    mkdir -p "$TMPDIR"
    export TMPDIR
fi

FASTQC_OUT_1="FastQC_Init"
FASTQC_RAW_1="$TMPDIR"

echo "Temporary directory is" $TMPDIR ", which will be deleted upon completion of the job"

echo "Final files will be saved to" $RESULT_DIRECTORY

echo "Now we are in" $PWD

#This part concatenates lanes as needed.
sampnames=`echo $FASTA_FILE | awk 'BEGIN{ FS="_S"} {print $1}' | uniq`
for i in $sampnames;
  do for j in 1 2;
     do
     cp ${i}_*_R${j}*.fast* $TMPDIR
     cat $TMPDIR/${i}_*_R${j}_001.fastq.gz > $TMPDIR/${i}_R${j}_001.fastq.gz
     rm -f $TMPDIR/${i}_*_R${j}_001.fastq.gz
  done
done

cd $TMPDIR

echo "These are the files now in tmpdir"

ls

# Copy Metagenomics files to Scratch space.
##cp -r ${FASTA_DIRECTORY}/${FASTA_FILE} .

# Ensure the output directory always exists (setting up before loading CBI and the fastqc software)
mkdir -p "${FASTQC_OUT_1}"

module load CBI fastqc

module list

## Instead of running files ending in fastq, some files will also end in fastq.gz. Therefore, I am creating a variable called `files` that identifies any file in the FASTA_DRECTORY that includes ".fastq", which will include both ".fastq" files and ".fastq.gz" files.
files=`ls | grep "[.]fastq"` # Gets any file in the directory with fastq in the name.

echo $files

run_fastqc () {
        for f in $files; do
        echo "Running FastQC on" $f
                fastqc \
                -t $fastqc_threads \
                "$f" \
                --nogroup \
                --outdir $FASTQC_OUT_1 \
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

MEMORY="-Xmx50g"
PARAMETERS_ADAPT_TRIM="ktrim=r k=23 mink=11 hdist=1 tbo tpe" # can also include tbo (trim adapters based on pair overlap detection using BBMerge [which does not require known adapter sequnces]) or tpe (specifies to trim both reads to the same length in the event that an adapter kmer was only detected in one of them)

## Current function uses the $files list generated above to ID only the Read1 files, confirms there's a Read2 file, and then continues.

run_bbtools() {
for f in $files; do
if [[ $f == *"_R1"* ]] && test -f "${f/_R1/_R2}"; then
echo "Running BBDuk on" $f ";" "${f/_R1/_R2}" "exists."

echo "Adapter trimming for" $f
## Adapter Trimming
singularity exec --bind $PWD:$PWD ${software_location}/bbtools.img bbduk.sh \
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
singularity exec --bind $PWD:$PWD ${software_location}/bbtools.img bbduk.sh \
"${MEMORY}" \
threads=$bbduk_threads \
in1="${BBDUK_DIR}"/"${TRIMMED_DIR}"/"${f/_R1/_R1_adapt_trim}" \
in2="${BBDUK_DIR}"/"${TRIMMED_DIR}"/"${f/_R1/_R2_adapt_trim}" \
out1="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R1_remPhiX}" \
out2="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R2_remPhiX}" \
ref=/bbmap/resources/phix174_ill.ref.fa.gz \
k=31 hdist=1 stats="${BBDUK_DIR}"/"${PhiX_DIR}"/phix_stats_${f/.fastq/}.txt

echo "Quality Filtering and Final Quality Histograms for" $f
## Quality Filtering/Trimming
singularity exec --bind $PWD:$PWD ${software_location}/bbtools.img bbduk.sh \
"${MEMORY}" \
threads=$bbduk_threads \
in1="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R1_remPhiX}" \
in2="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R2_remPhiX}" \
out1="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R1_qTrim}" \
out2="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R2_qTrim}" \
k=23 hdist=1 qtrim=rl trimq=15 maq=20 minlen=50 \ ## Upped MAQ from 10 to 20 and added minlen of 50, which seems reasonable
qhist="${BBDUK_DIR}"/"${QualFilt_DIR}"/qhist_${f/.fastq/}.txt

echo "Removing Human Sequences for" $f
singularity exec -B $PWD:$PWD,"${software_location}:/mnt" ${software_location}/bbtools.img bbmap.sh \
"${MEMORY}" \
threads=$bbduk_threads \
in="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R1_qTrim}" \
in2="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R2_qTrim}" \
outm1="${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_hg38map}" \
outu1="${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_clean}" \
outm2="${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R2_hg38map}" \
outu2="${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R2_clean}" \
path=/mnt

## Potentially add bowtie2 filtering here as well

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

# Determines distribution of threads for initial processes, based on number of cores alloted for the submission (realistically, in order to get these done in a reasonable amount of time, all of the processing power should be focused on the bbduk commands).
	let cores=$NSLOTS
	let fastqc_threads=1
	let bbduk_threads=cores-fastqc_threads

	echo $NSLOTS "cores detected; running FastQC on" `expr $NSLOTS - $bbduk_threads` "core(s) & quality processing on" ${bbduk_threads} "core(s)"
        run_fastqc &
        run_bbtools &
        wait ${!} ## This will wait for the two functions above to complete before continuing on.
else
        let cores=$NSLOTS
	let fastqc_threads=$NSLOTS
	let bbduk_threads=$NSLOTS
        echo $NSLOTS "core(s) is less than or equal to 3. Will not run FastQC and quality processing in parallel, but will use the specified number of cores for each process. If you wish to run these processes in parallel, change the -pe smp option to a number greater than 3."
        run_fastqc
        run_bbtools
fi 

# Perform FastQC of the human-filtered reads to confirm number and quality going into downstream analysis.
mkdir -p FastQC_Final

for f in bbduk_ReadCleaning_Output/HumanFiltering/*clean*; do
   echo "Running final FastQC on" $f
      fastqc \
         -t $fastqc_threads \
         $f \
         --nogroup \
         --outdir FastQC_Final \
         --dir "$PWD"
   echo "Finished running final FastQC on" $f
done

cd $TMPDIR

echo "PIPELINE COMPLETE! Moving files back to" $RESULT_DIRECTORY

rm $files
## This keeps all of the text output provided by the pipeline, but removes the intermediate fastqs to keep space utilization down.
rm bbduk_ReadCleaning_Output/AdapterTrimming/*fastq*
rm bbduk_ReadCleaning_Output/QualitySeqs/*fastq*
rm bbduk_ReadCleaning_Output/PhiXRemoval/*fastq*

## Adding a job ID to the directory that gets saved.
mkdir -p "$RESULT_DIRECTORY"/metagenomics_results_${JOB_ID}/$sampnames
mv ./* "$RESULT_DIRECTORY"/metagenomics_results_${JOB_ID}/$sampnames
