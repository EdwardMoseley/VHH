#!/bin/bash
#SBATCH -J makeblastdb #job name
#SBATCH --time=01-02:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 1 #1 nodes
#SBATCH -n 4 #2 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=32g #requesting 2GB of RAM total
#SBATCH --output=%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=edward.moseley@tufts.edu

#PRESTO_BIN=/cluster/tufts/cowenlab/immcantation_suite/presto/bin
CHANGEO_BIN=/cluster/tufts/cowenlab/immcantation_suite/changeo/bin
NCBI_BIN=/cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/bin

# Move to directory to see if that makes it work
cd /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0

# IgBLAST

${NCBI_BIN}/igblastn -germline_db_V /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/internal_data/alpaca/alpaca_V -germline_db_D /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/internal_data/alpaca/alpaca_D -germline_db_J /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/internal_data/alpaca/alpaca_J -domain_system imgt -ig_seqtype Ig -organism alpaca -auxiliary_data /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/optional_file/alpaca_gl.aux -outfmt '7 std qseq sseq btop' -query /cluster/tufts/cowenlab/emosel01/data/assembly/5_S3_atleast-1.fasta -out /cluster/tufts/cowenlab/emosel01/data/assembled/5_S3_atleast-1.fmt7

# STEP 2: MakeDB

${CHANGEO_BIN}/MakeDb.py igblast -i /cluster/tufts/cowenlab/emosel01/data/assembled/5_S3_atleast-1.fmt7 -s /cluster/tufts/cowenlab/emosel01/data/assembly/5_S3_atleast-1.fasta -r /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/database/alpaca/IG_DNA/imgt_gapped_alpaca_IGHV.fasta /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/database/alpaca/IG_DNA/imgt_gapped_alpaca_IGHD.fasta /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/database/alpaca/IG_DNA/imgt_gapped_alpaca_IGHJ.fasta --outdir /cluster/tufts/cowenlab/emosel01/data/assembled/ --outname 5_S3 --failed --log /cluster/tufts/cowenlab/emosel01/code/MakeDb.log --format airr --extended