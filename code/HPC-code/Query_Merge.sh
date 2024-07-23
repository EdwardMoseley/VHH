#!/bin/bash
#SBATCH -J query_merge #job name
#SBATCH --time=01-02:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 2 #1 nodes
#SBATCH -n 4 #2 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=32g #requesting 2GB of RAM total
#SBATCH --output=query_merge.%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=query_merge.%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=edward.moseley@tufts.edu

#[commands_you_would_like_to_exe_on_the_compute_nodes]

#Run R script
#module purge
module load R/4.3.0

Rscript /cluster/tufts/cowenlab/emosel01/code/Query_Merge.R
