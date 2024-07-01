#!/bin/bash
#SBATCH -J tofasta #job name
#SBATCH --time=06-23:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 1 #1 nodes
#SBATCH -n 1 #1 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=32g #requesting 2GB of RAM total
#SBATCH --output=tofasta.%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=tofasta.%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=edward.moseley@tufts.edu

module load R/4.3.0

Rscript tofasta.R

