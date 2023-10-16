library(ShortRead)

# Rscript /Users/emoseley/VHH/code/IMGT.R A2A_S1_Proc-BC1.fastq

# Arguments
args <- c("L2G2_S13_Proc-BC1.fastq")
# args = commandArgs(trailingOnly=TRUE)

# Directory
directory <- "/Users/emoseley/VHH/data/"

# Set current working directory to ncbi-igblast-1.17.1
setwd("~/ncbi-igblast-1.17.1")

file_origin <- paste("~/ncbi-igblast-1.17.1/", args[1], sep = '')
file_output <- gsub(".fastq", ".fasta", file_origin)
results <- gsub(".fastq", '.fmt7', args)

# load the fast q
tmp <- ShortRead::readFastq(file_origin)

# write as fasta
ShortRead::writeFasta(tmp, file = file_output)

outfmt <- '\'7 std qseq sseq btop\''

command <- paste("bin/igblastn -germline_db_V database/alpaca_vhh_clean/alpaca_V -germline_db_D database/alpaca_vhh_clean/alpaca_D -germline_db_J database/alpaca_vhh_clean/alpaca_J -domain_system imgt -ig_seqtype Ig -organism alpaca -auxiliary_data optional_file/alpaca_gl.aux -outfmt ", outfmt, " -query ", file_output, " -out ", paste("~/ncbi-igblast-1.17.1/", results, sep = ''), sep = '')

print(command)

system(command)
