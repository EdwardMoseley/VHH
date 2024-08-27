#!/bin/bash
#SBATCH -J assembly #job name
#SBATCH --time=06-23:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 1 #1 nodes
#SBATCH -n 1 #1 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=64g #requesting 2GB of RAM total
#SBATCH --output=8125.assembly.%j.%N.out #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=8125.assembly.%j.%N.err #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=edward.moseley@tufts.edu

## ###
## VHH-Assembly.sh
## Relies heavily on runalignmentpipeline.sh
## (Written by FAJO, of the Cowen Lab)
## ###

# Define environmental variables

PRESTO_BIN=/cluster/tufts/cowenlab/immcantation_suite/presto/bin
DATA_BIN=/cluster/tufts/cowenlab/emosel01/data/240520_0044_Charles_Shoemaker-8125/fastq_Lane1
CHANGEO_BIN=/cluster/tufts/cowenlab/immcantation_suite/changeo/bin
NCBI_BIN=/cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/bin
OUT_BIN=/cluster/tufts/cowenlab/emosel01/data/240520_0044_Charles_Shoemaker-8125/fastq_Lane1


# Iterate over each file in the assembly folder
# Pipe through each unique file name after trimming everyhing after the "_L"
for sample in $(ls $DATA_BIN | awk -F'_L0' '{print $1}' | uniq); do
    echo $sample
    # Create all necessary file names
    SEQ_RONE=${DATA_BIN}/${sample}_L001_R1_001.fastq
    SEQ_RTWO=${DATA_BIN}/${sample}_L001_R2_001.fastq
    echo 'Processing the following files: '
    echo $SEQ_RONE
    echo $SEQ_RTWO
    SEQ_QUAL=${DATA_BIN}/${sample}_assemble-pass.fastq
    echo $SEQ_QUAL
    SEQ_MPV=${DATA_BIN}/${sample}_quality-pass.fastq
    echo $SEQ_MPV
    SEQ_MPC=${DATA_BIN}/${sample}_MPV_primers-pass.fastq
    echo $SEQ_MPC
    SEQ_COL=${DATA_BIN}/${sample}_MPC_primers-pass.fastq
    echo $SEQ_COL
    SEQ_DUP=${DATA_BIN}/${sample}_collapse-unique.fastq
    echo $SEQ_DUP
    SEQ_FASTA=${DATA_BIN}/${sample}_atleast-1.fasta
    echo $SEQ_FASTA
    SEQ_FMT7=${OUT_BIN}/${sample}_atleast-1.fmt7
    echo $SEQ_FMT7
    # STEP 1: Alignment
    ${PRESTO_BIN}/AssemblePairs.py align -1 $SEQ_RTWO -2 $SEQ_RONE --coord illumina --rc tail --maxerror 0.2 --outname $sample
    # STEP 2: Quality Filtering
    ${PRESTO_BIN}/FilterSeq.py quality -s $SEQ_QUAL -q 20 --outname $sample
    # STEP 3: Variable Primer Pass
    ${PRESTO_BIN}/MaskPrimers.py score -s $SEQ_MPV -p /cluster/tufts/cowenlab/emosel01/data/AlpVHH_R2_primers_05Feb24.fasta --start 0 --mode cut --pf VPRIMER --outname ${sample}_MPV
    # STEP 4: Constant Primer Pass
    ${PRESTO_BIN}/MaskPrimers.py score -s $SEQ_MPC -p /cluster/tufts/cowenlab/emosel01/data/AlpVHH_R1_primers_05Feb24.fasta --start 0 --mode cut --revpr --pf CPRIMER --outname ${sample}_MPC
    # STEP 5: Deduplicate
    ${PRESTO_BIN}/CollapseSeq.py -s $SEQ_COL -n 20 --inner --uf CPRIMER --cf VPRIMER --act set --outname $sample
    # STEP 5: Split by duplicate count
    ${PRESTO_BIN}/SplitSeq.py group -s $SEQ_DUP -f DUPCOUNT --num 1 --outname $sample
    # Perform this same operation again, but this time we want a .fasta output file
    ${PRESTO_BIN}/SplitSeq.py group -s $SEQ_DUP -f DUPCOUNT --num 1 --fasta --outname $sample
    # Move to directory to see if that makes it work
    cd /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0

    # Now we want to use ncbi
    ${NCBI_BIN}/igblastn -germline_db_V  /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/internal_data/alpaca/alpaca_V -germline_db_D /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/internal_data/alpaca/alpaca_D -germline_db_J /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/internal_data/alpaca/alpaca_J -domain_system imgt -ig_seqtype Ig -organism alpaca -auxiliary_data /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/optional_file/alpaca_gl.aux -outfmt '7 std qseq sseq btop' -query $SEQ_FASTA -out $SEQ_FMT7
    # Now we build the database
    ${CHANGEO_BIN}/MakeDb.py igblast -i $SEQ_FMT7 -s $SEQ_FASTA -r /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/database/alpaca/IG_DNA/imgt_gapped_alpaca_IGHV.fasta  /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/database/alpaca/IG_DNA/imgt_gapped_alpaca_IGHD.fasta /cluster/tufts/cowenlab/emosel01/ncbi-igblast-1.22.0/database/alpaca/IG_DNA/imgt_gapped_alpaca_IGHJ.fasta --outdir $OUT_BIN --outname $sample --failed --log ${DATA_BIN}/MakeDb_${sample}.log --format airr --extended

    # Now we'll move back to the code directory
    cd /cluster/tufts/cowenlab/emosel01/code
done
