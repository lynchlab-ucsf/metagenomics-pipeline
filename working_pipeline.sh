#!/usr/bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l mem_free=20G
#$ -l h_rt=08:00:00
#$ -l s_rt=08:00:00
#$ -m bae
#$ -M kathryn.mccauley@ucsf.edu

## This version will be used for testing of code (fasta_directory points to a directory with only one sample in it). Parts that are working will be moved to the "combine_scripts.sh".

## Metagenomics Pipeline
## Current Version Written by Katie McCauley, 25 Nov 2019 with Ariane Panzer; last revised 03 Dec 2019

# Where do your metagenomics fastq files live?
FASTA_DIRECTORY="/wynton/group/lynch/kmccauley/01_raw/test_dir/"



######## All variables beyond this point do not need to be set.


FASTQC_OUT_1="FastQ_Out_1"
FASTQC_RAW_1="$TMPDIR"

## Define Temporary Directory Location (as scratchspace).
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

save_cwd=$PWD
echo "Current working directory is" $save_cwd


cd $TMPDIR

# Copy Metagenomics files to Scratch space.
cp -r ${FASTA_DIRECTORY}/* .

# Ensure the output directory always exists (setting up before loading CBI and the fastqc software)
mkdir -p "${FASTQC_OUT_1}"

module load CBI fastqc
## Will print information about the version of fastqc being used.
module list

## Instead of running files ending in fastq, some files will also end in fastq.gz. Therefore, I am creating a variable called `files` that identifies any file in the FASTA_DRECTORY that includes ".fastq", which will include both ".fastq" files and ".fastq.gz" files.
files=`ls | grep "[.]fastq"` # Gets any file in the directory with fastq in the name.

run_fastqc () {
        for f in $files; do
        echo "Running FastQC on" $f # Print file name to remaine apprised of progress
                fastqc \
                -t $fastqc_threads \
                "${f}" \
                --nogroup \
                --outdir "${TMPDIR}"/"FASTQC_Results/" \
                --dir "${TMPDIR}"
        echo "Finished running FastQC on" $f
        done
}

mv ${TMPDIR}/"FASTQC_Results" "${FASTA_DIRECTORY}/${FASTQC_OUT_1}"

# Definiting Output Directory names and creating said directories
BBDUK_DIR="bbduk_ReadCleaning_Output"
TRIMMED_DIR="AdapterTrimming"
HUMAN_DIR="HumanReads"
PhiX_DIR="PhiXRemoval"
QualFilt_DIR="QualitySeqs"

mkdir -p "${BBDUK_DIR}"/"${PhiX_DIR}"
mkdir -p "${BBDUK_DIR}"/"${QualFilt_DIR}"
mkdir -p "${BBDUK_DIR}"/"${HUMAN_DIR}"
mkdir -p "${BBDUK_DIR}"/"${TRIMMED_DIR}"

## Possibly useful parameters to look into and determine what should go up at the top.
MEMORY="-Xmx40g"
PARAMETERS_ADAPT_TRIM="ktrim=r k=23 mink=11 hdist=1" # can also include tbo (trim adapters based on pair overlap detection using BBMerge [which does not require known adapter sequnces]) or tpe (specifies to trim both reads to the same length in the event that an adapter kmer was only detected in one of them)

## Current function uses the $files list generated above to ID only the Read1 files, confirms there's a Read2 file, and then continues.

run_bbtools() {
for f in $files; do
if [[ $f == *"_R1_"* ]] && test -f "${f/_R1/_R2}"; then
echo "Running BBDuk on" $f ";" "${f/_R1/_R2}" "exists."

echo "Adapter trimming for" $f
## Adapter Trimming
singularity exec /wynton/group/lynch/kmccauley/mySoftware/bbtools.img bbduk.sh \
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
singularity exec /wynton/group/lynch/kmccauley/mySoftware/bbtools.img bbduk.sh \
"${MEMORY}" \
threads=$bbduk_threads \
in1="${BBDUK_DIR}"/"${TRIMMED_DIR}"/"${f/_R1/_R1_adapt_trim}" \
in2="${BBDUK_DIR}"/"${TRIMMED_DIR}"/"${f/_R1/_R2_adapt_trim}" \
out1="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R1_remPhiX}" \
out2="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R2_remPhiX}" \
ref=/bbmap/resources/phix174_ill.ref.fa.gz \
k=31 hdist=1 stats=stats.txt

echo "Quality Filtering for" $f
## Quality Filtering/Trimming
singularity exec /wynton/group/lynch/kmccauley/mySoftware/bbtools.img bbduk.sh \
"${MEMORY}" \
threads=$bbduk_threads \
in1="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R1_remPhiX}" \
in2="${BBDUK_DIR}"/"${PhiX_DIR}"/"${f/_R1/_R2_remPhiX}" \
out1="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R1_qTrim}" \
out2="${BBDUK_DIR}"/"${QualFilt_DIR}"/"${f/_R1/_R2_qTrim}" \
k=27 hdist=1 qtrim=rl trimq=17 cardinality=t

echo "Human Removal for" $f
singularity exec -B "/wynton/group/lynch/kmccauley/mySoftware:/mnt" /wynton/group/lynch/kmccauley/mySoftware/bbtools.img bbmap.sh \
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

# Determines distribution of threads for initial processes, based on number of cores alloted for the process
	let cores=$NSLOTS
	let fastqc_threads=cores/3
	let bbduk_threads=cores-fastqc_threads

	echo $NSLOTS "cores detected; running FastQC on" `expr $NSLOTS - $bbduk_threads` "core(s) & quality processing on" ${bbduk_threads} "core(s)"
        run_fastqc &
        run_bbtools &
        wait ## This will wait for the two functions above to complete before continuing on.
else
	let fastqc_threads=$NSLOTS
	let bbduk_threads=$NSLOTS
        echo $NSLOTS "core(s) is less than or equal to 3. Will not run FastQC and quality processing in parallel, but will use the specified number of cores for each process. If you wish to run these processes in parallel, change the -pe smp option to a number greater than 3."
        run_fastqc
        run_bbtools
fi 



## Generate contigs, using SPAdes. Need to figure out why MIDAS isn't working (and what MIDAS uses as input) -- both are having installation troubles.
# Also, there was discussion of pulling down genomes based on OTU table sequences, and using that as reference. Won't be helpful for my study, but something to figure out how to implement through Wynton and Unix. Would this be better executed by using a fasta file to find matches to all OTUs at 97% identity? Or use the names in the OTU taxonomies to pull, say, all genomes belonging to the tagged genus. Capabilities for the latter are there, but the former is going to require more work. Theoretically, it can be done with BLAST or SILVA using their GUIs, but it would be nice to be able to do this from the command line.

## Just some points for discussion/thought...

#Pulling some sequences using BLAST, but need to find something that works on the command line, and only requires a fasta file for input.

# This setting of the PATH variable works:
export PATH=$PATH:/wynton/group/lynch/kmccauley/mySoftware/SPAdes-3.11.1-Linux/bin/
## spades.py --test ## successfully executes.

## try spades on the one sample:
mkdir spades_contig_output

make_contigs() {
for f in $ files; do
if [[ $f == *"_R1_"* ]] && test -f "${f/_R1/_R2}"; then

metaspades.py -k 21,33,55,77 \   ## Check for something that allows for combination of the R1 and R2.
-1 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R1_clean}" \
-2 "${BBDUK_DIR}"/"${HUMAN_DIR}"/"${f/_R1/_R2_clean}" \
-o spades_contig_output/"${f/_R1/_contig_dat}"

fi
done
}

make_contigs

mkdir -r "$save_cwd"/final_results/
mv $TMPDIR "$save_cwd"/final_results/
