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
#SBATCH --mail-user=faith.ocitti@tufts.edu

CURRENT_DIR=$(pwd)
# SAMPLE=13_VecF1_R2_S13
SAMPLE=19_VecF1_R3_S19
PRESTO_BIN=/cluster/tufts/cowenlab/immcantation_suite/presto/bin
CHANGEO_BIN=/cluster/tufts/cowenlab/immcantation_suite/changeo/bin
DATA_BIN=/cluster/tufts/cowenlab/focitt01/sigh/alpaca/phaseI/240124_0629_Charles_Shoemaker-7981/fastq_Lane1
# RESULTS_DIR_30=/cluster/tufts/cowenlab/focitt01/sigh/alpaca/phaseI/231224_0629_Charles_Shoemaker-7981/results/$SAMPLE/30
RESULTS_DIR=/cluster/tufts/cowenlab/focitt01/sigh/alpaca/phaseI/240124_0629_Charles_Shoemaker-7981/results/$SAMPLE
LOG_DIR=/cluster/tufts/cowenlab/focitt01/sigh/alpaca/phaseI/240124_0629_Charles_Shoemaker-7981/logs
FASTA_DIR=${RESULTS_DIR}/fasta
FASTQ_DIR=${RESULTS_DIR}/fastq

if [ ! -d $LOG_DIR ]; then mkdir -p $LOG_DIR; fi
if [ ! -d $RESULTS_DIR ]; then mkdir -p $RESULTS_DIR; fi
if [ ! -d $FASTA_DIR ]; then mkdir -p $FASTA_DIR; fi
if [ ! -d $FASTQ_DIR ]; then mkdir -p $FASTQ_DIR; fi

SEQF=${DATA_BIN}/${SAMPLE}_L001_R1_001.fastq
SEQR=${DATA_BIN}/${SAMPLE}_L001_R2_001.fastq
SEQ1=${RESULTS_DIR}/AssemblePairs${SAMPLE}.fastq
SEQ2=${RESULTS_DIR}/FilterSeq${SAMPLE}.fastq
SEQ3=${RESULTS_DIR}/MaskVPrimer${SAMPLE}.fastq
SEQ4=${RESULTS_DIR}/MaskCPrimer${SAMPLE}.fastq
SEQ5=${RESULTS_DIR}/CollapseSeq${SAMPLE}.fastq
TAB1=${RESULTS_DIR}/ParseHeaders${SAMPLE}.tab

# STEP 1: Alignment

# ${PRESTO_BIN}/AssemblePairs.py align -1 $SEQF -2 $SEQR --maxerror 0.2 --coord illumina  -o ${SEQ1} --rc tail --log ${LOG_DIR}/AssemblePairs${SAMPLE}.log
# throw out low quality reads, anything below Q30 gets tossed
# ${PRESTO_BIN}/FilterSeq.py quality -s $SEQ1 -q 20 -o $SEQ2 --log ${LOG_DIR}/FilterSeq${SAMPLE}.log
# match the forward primers and read the forward half of the barcode
${PRESTO_BIN}/MaskPrimers.py align -s $SEQ2 -p ${DATA_BIN}/L1_AlpVHH.fasta --barcode --maxerror 0.01 --mode cut --pf VPRIMER  -o $SEQ3 --log ${LOG_DIR}/MaskVPrimer${SAMPLE}.log
# match the reverse primers and read the reverse half of the barcode
${PRESTO_BIN}/MaskPrimers.py align -s $SEQ3 -p ${DATA_BIN}/SH_AlpVHH.fasta --barcode --maxerror 0.01 --mode cut --revpr --pf CPRIMER  -o $SEQ4 --log ${LOG_DIR}/MaskCPrimer${SAMPLE}.log
# remove duplicate sequences and annotate with number of identical reads with DUPCOUNT
${PRESTO_BIN}/CollapseSeq.py -s $SEQ4 -n 20 --inner --act set -o $SEQ5 --cf VPRIMER --uf CPRIMER --log ${LOG_DIR}/CollapseSeq${SAMPLE}.log
# split sequences and output as fasta files
${PRESTO_BIN}/SplitSeq.py group -s $SEQ5 -f DUPCOUNT --num 1 --fasta --outname ${SAMPLE} --outdir ${FASTA_DIR}
# split sequences and output as fastq files
${PRESTO_BIN}/SplitSeq.py group -s $SEQ5 -f DUPCOUNT --num 1 --outname ${SAMPLE} --outdir ${FASTQ_DIR}
# output duplicate annotations as table
${PRESTO_BIN}/ParseHeaders.py table -s $SEQ5 -f ID DUPCOUNT CPRIMER VPRIMER -o $TAB1

# STEP 2: Create germlines

SEQ6=${FASTA_DIR}/${SAMPLE}_atleast-1.fasta
SEQ7=${RESULTS_DIR}/${SAMPLE}_db-pass.tsv
FMT7=${RESULTS_DIR}/${SAMPLE}.fmt7

# run IGBlast on the duplicate barcode sequences
DATA_DIR=/cluster/tufts/cowenlab/focitt01/sigh/alpaca/phaseI/data
NCBI_DIR=${DATA_DIR}/ncbi-igblast-1.22.0
# to account for a weird restriction in IGBlast, option to use IGDATA
cd $NCBI_DIR
NCBI_VDJ_DIR=${NCBI_DIR}/internal_data/alpaca
AUX_DATA_FILE=${NCBI_DIR}/optional_file/alpaca_gl.aux
# identitification of VDJ regions using IG blast nucleotide 
${NCBI_DIR}/bin/igblastn -germline_db_V ${NCBI_VDJ_DIR}/alpaca_V -germline_db_D ${NCBI_VDJ_DIR}/alpaca_D -germline_db_J ${NCBI_VDJ_DIR}/alpaca_J -auxiliary_data $AUX_DATA_FILE -domain_system imgt -ig_seqtype Ig -organism alpaca -outfmt '7 std qseq sseq btop' -query $SEQ6 -out $FMT7

# revert to initial working directory for funsies
cd $CURRENT_DIR

# write the IGBlast output to an AIRR-SEQ database file
IGMT_DATA=${DATA_DIR}/igmt
${CHANGEO_BIN}/MakeDb.py igblast -i $FMT7 -s $SEQ6 -r ${IGMT_DATA}/alpaca_IGHV.fasta ${IGMT_DATA}/alpaca_IGHD.fasta ${IGMT_DATA}/alpaca_IGHJ.fasta --outdir $RESULTS_DIR --outname $SAMPLE --failed --log ${LOG_DIR}/MakeDb${SAMPLE}.log --format airr --extended
: '
#filter for productive sequences, sequences that have a high quality v call
${CHANGEO_BIN}/ParseDb.py select -d $SEQ7 -f productive -u T --regex --outdir $RESULTS_DIR --outname $SAMPLE

SEQ8=${RESULTS_DIR}/${SAMPLE}_parse-select.tsv
THR_CALC_SCRIPT=/cluster/tufts/cowenlab/focitt01/sigh/alpaca/phaseI/scripts/calculateClonalThreshold.R 
THR=$(Rscript ${THR_CALC_SCRIPT} $SEQ8)
# define clones
${CHANGEO_BIN}/DefineClones.py -d $SEQ8 --act set --model ham --norm len --dist $THR --outname $SAMPLE --outdir $RESULTS_DIR --log ${LOG_DIR}/DefineClones${SAMPLE}.log --failed 


# create germlines
SEQ9=${RESULTS_DIR}/${SAMPLE}_clone-pass.tsv
${CHANGEO_BIN}/CreateGermlines.py -d $SEQ9 -g dmask --cloned -r ${IGMT_DATA}/alpaca_IGHV.fasta ${IGMT_DATA}/alpaca_IGHD.fasta ${IGMT_DATA}/alpaca_IGHJ.fasta --outdir $RESULTS_DIR --outname $SAMPLE --log ${LOG_DIR}/CreateGermlines${SAMPLE}.log --failed
SEQ10=${RESULTS_DIR}/${SAMPLE}_germ-pass.tsv
# build phylogenetic trees
${CHANGEO_BIN}/BuildTrees.py -d $SEQ10 --outname $SAMPLE --outdir $RESULTS_DIR --log ${LOG_DIR}/BuildTrees${SAMPLE}.log --collapse --igphyml --clean all
