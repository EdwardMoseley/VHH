#!/bin/bash
#SBATCH -J cleanclusters #job name
#SBATCH --time=06-23:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 1 #1 nodes
#SBATCH -n 1 #2 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=32g #requesting 2GB of RAM total
#SBATCH --output=cleanclust.%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=cleanclust.%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=edward.moseley@tufts.edu

#[commands_you_would_like_to_exe_on_the_compute_nodes]

# INPUT_FILE=GS_d2_18Jun24_min_threshold_clone-pass.csv

#Run R script
#module purge
module load R/4.3.0

#Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterPool.R GS_d2_18Jun24_min_threshold_clone-pass.csv
#Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClustertoAA.R GS_d2_18Jun24_pooled.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterMMseqs2.R GS_d2_18Jun24_pooled_VHH_database.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterRefining.R GS_d2_18Jun24_pooled_subclustered_database.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterRecluster.R GS_d2_18Jun24_pooled_outlier_database.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterUnion.R GS_d2_18Jun24_pooled_reclustered_database.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterPool.R GS_d2_18Jun24_results.csv

