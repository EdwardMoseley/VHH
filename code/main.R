## ###
## Immcantation Pipeline
## Edward Moseley, Cowen Lab
## Tufts University
## 18Sep23
## ###

## ###
## Input files:

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
# args = commandArgs(trailingOnly=TRUE)

args <- c("A2A_S1_L001_R1_001.fastq", "A2A_S1_L001_R2_001.fastq")

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please provide the directory to the data files (in .fastq format)", call.=FALSE)
} else {
  file_R1 <- args[1]
  file_R2 <- args[2]
}


## Command line argument structure
docker_prefix <- "docker run -v"
directory <- "/Users/emoseley/VHH/data:/data"
data_mount <- "/data/"
immcantation <- "immcantation/suite:4.4.0"

# Input file metadata
metadata <- strsplit(file_R1, '_')[[1]]

# Create ID for file data
id <- paste(metadata[1], metadata[2], sep = '_')

## Initial data for command structure
program <- "AssemblePairs.py"

logfile <- paste(data_mount, id, '_AP.log', sep = '')

arguments <- paste("--coord illumina --rc tail --maxerror 0.1 --outname ",
                   id, " --log ", logfile, sep = '')

## Paired-End Assembly (Takes a long time).
command <- paste(docker_prefix, directory, immcantation, program, 'align -1', data_mount, args[2],
                 '-2 ', data_mount, args[1], arguments, sep = ' ')

#system_command <- "docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 \
#AssemblePairs.py align -1 /data/A2A_S1_L001_R2_001.fastq -2 /data/A2A_S1_L001_R1_001.fastq \
#--coord illumina --rc tail --maxerror 0.1 --outname ALPACA1 --log /data/ALPACA_AP.log"


# Parse log

system_command <- paste(docker_prefix, directory, immcantation, "ParseLog.py -l", logfile, 
                        '-f ID LENGTH OVERLAP ERROR PVALUE', sep = ' ')

system(system_command)
