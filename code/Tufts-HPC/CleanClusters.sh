#!/bin/bash
#SBATCH -J CleanClusters #job name
#SBATCH --time=01-02:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 2 #1 nodes
#SBATCH -n 4 #2 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=32g #requesting 2GB of RAM total
#SBATCH --output=%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=edward.moseley@tufts.edu

#[commands_you_would_like_to_exe_on_the_compute_nodes]

#Run R script
#module purge
module load R/4.3.0
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterPool.R 21_VecF4_R3_S21_min_threshold_clone-pass.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClustertoAA.R 21_VecF4_R3_S21-rule-pooled-families.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterMMseqs2.R 21_VecF4_R3_S21-rule-pooled-families_VHH_database.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterRefining.R 21_VecF4_R3_S21-rule-pooled-families_subclustered_database.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterRecluster.R 21_VecF4_R3_S21-pooled-families_subclustered_database.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanClusterUnion.R 21_VecF4_R3_S21-rule-pooled-families_subclustered-reclustered-database.csv
Rscript /cluster/tufts/cowenlab/emosel01/code/CleanCluster-Pool.R 21_VecF4_R3_S21-Pooled-Cleaning-Reclustering-Initial-Results.csv
