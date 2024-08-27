#!/bin/bash
#SBATCH -J mash #job name
#SBATCH --time=06-23:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 1 #1 nodes
#SBATCH -n 1 #1 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=64g #requesting 2GB of RAM total
#SBATCH --output=mash.%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=mash.%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=edward.moseley@tufts.edu

/cluster/home/emosel01/mash/mash-Linux64-v2.3/mash sketch -o ../data/sketches/genomes /cluster/tufts/cowenlab/emosel01/data/XIV/*db*.fasta

/cluster/home/emosel01/mash/mash-Linux64-v2.3/mash dist ../data/sketches/genomes.msh ../data/sketches/genomes.msh > MASH_XIV_distances25Aug24.tsv
