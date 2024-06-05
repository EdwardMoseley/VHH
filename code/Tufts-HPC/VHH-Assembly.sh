#!/bin/bash
#SBATCH -J alignmentrerun #job name
#SBATCH --time=01-02:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 2 #1 nodes
#SBATCH -n 4 #2 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=64g #requesting 2GB of RAM total
#SBATCH --output=AlignmentJob.%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=AlignmentJob.%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=edward.moseley@tufts.edu

## ###
## VHH-Assembly.sh
## Relies heavily on runalignmentpipeline.sh
## (Written by FAJO, of the Cowen Lab)
## ###

# Define environmental variables

PRESTO_BIN=/cluster/tufts/cowenlab/immcantation_suite/presto/bin
DATA_BIN=/cluster/tufts/cowenlab/emosel01/data/assembly


# ${PRESTO_BIN}/FilterSeq.py quality -s /cluster/tufts/cowenlab/emosel01/data/S3_assemble-pass.fastq -q 20 --outname S3

# Iterate over each file in the assembly folder
# Pipe through each unique file name after trimming everyhing after 2nd underscore
for sample in $(ls $DATA_BIN | cut -d'_' -f1,2 | uniq); do
    echo $sample
    # Create all necessary file names
    SEQ_RONE=${DATA_BIN}/${sample}_L001_R1_001.fastq
    SEQ_RTWO=${DATA_BIN}/${sample}_L001_R2_001.fastq
    echo 'Processing: '
    echo $SEQ_RONE
    echo $SEQ_RTWO
    SEQ_QUAL=${DATA_BIN}/${sample}_assemble-pass.fastq
    echo $SEQ_QUAL
   # STEP 1: Alignment
   ${PRESTO_BIN}/AssemblePairs.py align -1 $SEQ_RTWO -2 $SEQ_RONE --coord illumina --rc tail --maxerror 0.2 --outname $sample
   ${PRESTO_BIN}/FilterSeq.py quality -s $SEQ_QUAL -q 20 --outname $sample
done

