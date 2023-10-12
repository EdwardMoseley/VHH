## ###
## Immcantation Pipeline
## Edward Moseley, Cowen Lab
## Tufts University
## 11Oct23
## ###

## Commands run

# Rscript /Users/emoseley/VHH/code/main.R A2A_S1_L001_R1_001.fastq A2A_S1_L001_R2_001.fastq
# Rscript /Users/emoseley/VHH/code/main.R A2G1_S6_L001_R1_001.fastq A2G1_S6_L001_R2_001.fastq
# Rscript /Users/emoseley/VHH/code/main.R A2G2_S12_L001_R1_001.fastq A2G2_S12_L001_R2_001.fastq

## Cammands running

# Rscript /Users/emoseley/VHH/code/main.R L2A_S2_L001_R1_001.fastq  L2A_S2_L001_R2_001.fastq
# Rscript /Users/emoseley/VHH/code/main.R L2G1_S7_L001_R1_001.fastq  L2G1_S7_L001_R2_001.fastq

## Commands to run

# Rscript /Users/emoseley/VHH/code/main.R L2G2_S13_L001_R1_001.fastq  L2G2_S13_L001_R2_001.fastq


## ###
## Input files:
## ### 

# Alexander
# [1] "A2A_S1_L001_R1_001.fastq"
# [2] "A2A_S1_L001_R2_001.fastq"
# [3] "A2G1_S6_L001_R1_001.fastq"
# [4] "A2G1_S6_L001_R2_001.fastq"
# [5] "A2G2_S12_L001_R1_001.fastq"
# [6] "A2G2_S12_L001_R2_001.fastq"

# Lara
# [8] "L2A_S2_L001_R1_001.fastq"
# [9] "L2A_S2_L001_R2_001.fastq"
# [10] "L2G1_S7_L001_R1_001.fastq"
# [11] "L2G1_S7_L001_R2_001.fastq"
# [12] "L2G2_S13_L001_R1_001.fastq"
# [13] "L2G2_S13_L001_R2_001.fastq"

## Note use: Rscript main.R [sample_name.fastq]

# Parse arguments into
args = commandArgs(trailingOnly=TRUE)

# args <- c("A2A_S1_L001_R1_001.fastq", "A2A_S1_L001_R2_001.fastq")

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please provide the directory to the data files (in .fastq format)", call.=FALSE)
} else {
  file_R1 <- args[1]
  file_R2 <- args[2]
}

## Command line argument structure
docker_prefix <- " docker run -v"
directory <- " /Users/emoseley/VHH/data:/data "
immcantation <- " immcantation/suite:4.4.0 "
data_mount <- "/data/"

# Input file metadata
metadata <- strsplit(file_R1, '_')[[1]]

# Create ID for file data
id <- paste(metadata[1], metadata[2], 'Proc', sep = '_')


## ###
## (1) Paired-End Assembly (Takes a long time).
## ###

# Logfile name
logfile <- paste(data_mount, id, '_AP.log', sep = '')

# Initial data for command structure
program <- " AssemblePairs.py "

# Arguments
arguments <- paste(" --coord illumina --rc tail --maxerror 0.1 --outname ",
                   id, " --log ", logfile, sep = '')

# Formulate command from components
print(command <- paste(docker_prefix, 
                 directory, 
                 immcantation, 
                 program, 
                 'align -1 ', 
                 data_mount, 
                 args[2],
                 ' -2 ', 
                 data_mount, 
                 args[1], 
                 arguments, sep = ''))



# Run command
system(command)

#START> AssemblePairs
#COMMAND> align
#FILE1> A2Asample_S1_L001_R2_001.fastq
#FILE2> A2Asample_S1_L001_R1_001.fastq
#COORD_TYPE> illumina
#ALPHA> 1e-05
#MAX_ERROR> 0.1
#MIN_LEN> 8
#MAX_LEN> 1000
#SCAN_REVERSE> False
#NPROC> 6

#PROGRESS> 15:28:38 |####################| 100% (1,250) 0.1 min
  
#OUTPUT> A2Asample_S1_Proc_assemble-pass.fastq
#PAIRS> 1250
#PASS> 660
#FAIL> 590
#END> AssemblePairs

# Parse log

program <- "ParseLog.py -l"

print(command <- paste(docker_prefix, directory, immcantation, program, logfile, 
                        '-f ID LENGTH OVERLAP ERROR PVALUE', sep = ' '))

system(command)

#START> ParseLog
#FILE> A2Asample_S1_Proc_AP.log

#PROGRESS> 15:29:52 (1250) 0.0 min

#OUTPUT> A2Asample_S1_Proc_AP_table.tab
#RECORDS> 1250
#PASS> 1250
#FAIL> 0
#END> ParseLog

## ###
## 2) (Quality Control Step) Filter Sequence
## ###

### Removing Low-Quality Reads
program <- "FilterSeq.py quality -s "
infile <- paste(data_mount, id, '_assemble-pass.fastq', sep = '')
arguments <- " -q 30 --outname "
logfile <- paste(' --log ', data_mount, id, '_FS.log', sep = '')


print(command <- paste(docker_prefix, 
                 directory, 
                 immcantation, 
                 program, 
                 infile,
                 arguments,
                 id,
                 logfile,
                 sep = ''))

system(command)

#START> FilterSeq
#COMMAND> quality
#FILE> A2Asample_S1_Proc_assemble-pass.fastq
#INNER> False
#MIN_QUAL> 30.0
#NPROC> 6

#PROGRESS> 15:52:36 |####################| 100% (660) 0.0 min
  
#OUTPUT> A2Asample_S1_Proc_quality-pass.fastq
#SEQUENCES> 660
#PASS> 541
#FAIL> 119
#END> FilterSeq

# Parse log
program <- "ParseLog.py -l "

logfile <- paste(data_mount, id, '_FS.log -f ID QUALITY', sep = '')

print(command <- paste(docker_prefix, 
                 directory,
                 immcantation,
                 program, logfile, sep = ' '))

system(command)

#START> ParseLog
#FILE> A2Asample_S1_Proc_FS.log
#
#PROGRESS> 16:03:35 (660) 0.0 min
#
#OUTPUT> A2Asample_S1_Proc_FS_table.tab
#RECORDS> 660
#PASS> 660
#FAIL> 0
#END> ParseLog

## ###
## Masking of primer regions
## ###

program <- "MaskPrimers.py score -s "
infile <- paste(data_mount, id, '_quality-pass.fastq', sep = '')
barcode <- " -p /data/AlpVHH.fasta --barcode --start 8 --mode cut --pf VPRIMER --maxerror 0.01"
outname <- paste(' --outname ', id, '-FWD', sep = '')
logfile <- paste(' --log ', data_mount, id, '-FWD-MPV.log ', sep = '')

print(command <- paste(docker_prefix,
                 directory,
                 immcantation,
                 program,
                 infile,
                 barcode,
                 outname,
                 logfile, sep = ''))


system(command)

#START> MaskPrimers
#COMMAND> score
#SEQ_FILE> A2Asample_S1_Proc_quality-pass.fastq
#PRIMER_FILE> AlpVHH.fasta
#MODE> cut
#MAX_ERROR> 0.01
#START_POS> 8
#REV_PRIMER> False
#BARCODE> True
#BARCODE_FIELD> BARCODE
#PRIMER_FIELD> VPRIMER
#NPROC> 6

#PROGRESS> 16:27:24 |####################| 100% (541) 0.0 min
  
#OUTPUT> A2Asample_S1_Proc-FWD_primers-pass.fastq
#SEQUENCES> 541
#PASS> 410
#FAIL> 131
#END> MaskPrimers

## Second one
program <- "MaskPrimers.py score -s "
infile <- paste(data_mount, id, '-FWD_primers-pass.fastq', sep = '')
barcode <- " -p /data/AlpVHH.fasta --barcode --start 8 --mode cut --revpr --pf CPRIMER --maxerror 0.01"
outname <- paste(' --outname ', id, '-REV', sep = '')
logfile <- paste(' --log ', data_mount, id, '-REV-MPC.log ', sep = '')

print(command <- paste(docker_prefix,
                 directory,
                 immcantation,
                 program,
                 infile,
                 barcode,
                 outname,
                 logfile, sep = ''))

system(command)

#START> MaskPrimers
#COMMAND> score
#SEQ_FILE> A2Asample_S1_Proc_quality-pass.fastq
#PRIMER_FILE> AlpVHH.fasta
#MODE> cut
#MAX_ERROR> 0.01
#START_POS> 8
#REV_PRIMER> True
#BARCODE> True
#BARCODE_FIELD> BARCODE
#PRIMER_FIELD> CPRIMER
#NPROC> 6

#PROGRESS> 16:30:13 |####################| 100% (541) 0.0 min
  
#OUTPUT> A2Asample_S1_Proc-REV_primers-pass.fastq
#SEQUENCES> 541
#PASS> 407
#FAIL> 134
#END> MaskPrimers

## Parse Log

program <- "ParseLog.py -l "

print(command <- paste(docker_prefix,
                 directory,
                 immcantation,
                 program,
                 paste(data_mount, id, '-FWD-MPV.log ', sep = ''), 
                 paste(data_mount, id, '-REV-MPC.log -f ID PRIMER ERROR', sep = ''), sep = ''))


system(command)

#START> ParseLog
#FILE> A2Asample_S1_Proc-FWD-MPV.log

#PROGRESS> 16:36:55 (541) 0.0 min

#OUTPUT> A2Asample_S1_Proc-FWD-MPV_table.tab
#RECORDS> 541
#PASS> 541
#FAIL> 0
#END> ParseLog

#START> ParseLog
#FILE> A2Asample_S1_Proc-REV-MPC.log

#PROGRESS> 16:36:55 (541) 0.0 min

#OUTPUT> A2Asample_S1_Proc-REV-MPC_table.tab
#RECORDS> 541
#PASS> 541
#FAIL> 0
#END> ParseLog

## ###
## (4) Parse Headers
## ###


program <- "ParseHeaders.py collapse -s "
# Arguments includes input file data
arguments <- paste(data_mount, id, '-REV_primers-pass.fastq -f BARCODE --act cat', sep = '')


print(command <- paste(docker_prefix, 
                 directory, 
                 immcantation, 
                 program, 
                 arguments, sep = ''))

system(command)


#START> ParseHeaders
#COMMAND> collapse
#FILE> A2Asample_S1_Proc-REV_primers-pass.fastq
#ACTIONS> cat
#FIELDS> BARCODE

#PROGRESS> 16:45:43 |####################| 100% (407) 0.0 min
  
#OUTPUT> A2Asample_S1_Proc-REV_primers-pass_reheader.fastq
#SEQUENCES> 407
#END> ParseHeaders

## ###
## (5) Build Consensus Sequences
## ###

program <- "BuildConsensus.py -s "
infile <- paste(data_mount, id, '-REV_primers-pass_reheader.fastq', sep = '')
arguments <- " --pf VPRIMER --maxerror 0.1 -n 2 --bf BARCODE -o "
outfile <- paste(data_mount, id, '-BC1.fastq', sep = '')

print(command <- paste(docker_prefix, 
                 directory, 
                 immcantation,
                 program,
                 infile,
                 arguments,
                 outfile, sep = ''))

system(command)


#START> BuildConsensus
#FILE> A2Asample_S1_Proc-REV_primers-pass_reheader.fastq
#BARCODE_FIELD> BARCODE
#MIN_COUNT> 2
#MIN_FREQUENCY> 0.6
#MIN_QUALITY> 0
#MAX_GAP> None
#PRIMER_FIELD> VPRIMER
#PRIMER_FREQUENCY> None
#MAX_ERROR> 0.1
#MAX_DIVERSITY> None
#DEPENDENT> False
#COPY_FIELDS> None
#COPY_ACTIONS> None
#NPROC> 6

#ERROR> Processing sequence set with ID: ACGACCTT


## ###
## IMGT Database (Gene Annotation)
## ###


#### Using Database

#It looks like `igblastn` requires `.fasta` files, so we will convert our file into that format.

# load the fastq
# tmp <- ShortRead::readFastq("~/VHH/data/ALPACA1-BC1.fastq")

# write as fasta
# ShortRead::writeFasta(tmp, file = "~/VHH/data/ALPACA1-BC1.fasta")


## Move results to a new directory
current_wd <- gsub(' ', '', strsplit(directory, ':')[[1]][1])

# Set working directory
setwd(current_wd)

# Pull results files by their ID var
# results <- list.files()[grepl(id, list.files())]

# Make directory
system(paste("mkdir ", id, sep = ''))

#file.copy(from=results, to=id, 
#          overwrite = TRUE, recursive = FALSE, 
#          copy.mode = TRUE)

# Cleanup

# next, mkdir id
# cp id* /id
