
# Pull in utility functions
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

root_dir <- "/cluster/tufts/cowenlab/emosel01/data/2024_TUCF_7981/"

files <- list.files(root_dir)

# May be a good idea to change this script to account for different file types (.fastq, .tsv, etc.)

for (f in files){
  cat("Loading ", paste(root_dir, f, sep = ''), '\n')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)
  cat("Writing: ", paste(root_dir, gsub(".csv", ".fasta", f), sep = ''), '\n')
  write_fasta(tmp, "sequence", "sequence_id", paste(root_dir, gsub(".csv", ".fasta", f), sep = ''))
}



