#!/bin/bash
#SBATCH -J alignmentrerun #job name
#SBATCH --time=01-02:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 2 #1 nodes
#SBATCH -n 4 #2 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=64g #requesting 2GB of RAM total
#SBATCH --output=19AlignmentJobRerun.%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=19AlignmentJobRerun.%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email optitions
#SBATCH --mail-user=edward.moseley@tufts.edu

PRESTO_BIN=/cluster/tufts/cowenlab/immcantation_suite/presto/bin

# STEP 1: Alignment

# STEP 2: Quality Filtering
              
${PRESTO_BIN}/FilterSeq.py quality -s /cluster/tufts/cowenlab/emosel01/data/S3_assemble-pass.fastq -q 20 --outname S3