#!/bin/bash
#SBATCH -J NGO0265_hclust #job name
#SBATCH --time=06-23:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 1 #1 nodes
#SBATCH -n 1 #1 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=256g #requesting 2GB of RAM total
#SBATCH --output=NGO0265_hclust.%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=NGO0265_hclust.%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=edward.moseley@tufts.edu

#[commands_you_would_like_to_exe_on_the_compute_nodes]

#Run R script
#module purge
module load R/4.3.0
Rscript /cluster/tufts/cowenlab/emosel01/code/HierarchicalClustering.R
