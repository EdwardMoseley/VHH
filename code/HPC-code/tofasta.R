
# Pull in utility functions
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

root_dir <- "/cluster/tufts/cowenlab/emosel01/data/hclust/"

files <- list.files(root_dir)

for (f in files){
  cat("Loading ", paste(root_dir, f, sep = ''), '\n')
  tmp <- read.csv(paste0(root_dir, f), sep = '\t', header = T, stringsAsFactors = F)
  cat("Writing: ", paste(root_dir, gsub(".tsv", ".fasta", f), sep = ''), '\n')
  write_fasta(tmp, "sequence", "sequence_id", paste(root_dir, gsub(".tsv", ".fasta", f), sep = ''))
}



