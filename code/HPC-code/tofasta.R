
# Pull in utility functions
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

root_dir <- "/cluster/tufts/cowenlab/emosel01/data/XIV/"

files <- list.files(root_dir)[grepl(".tsv", list.files(root_dir))]
print("Currently Processing Files: ")
print(files)

# May be a good idea to change this script to account for different file types (.fastq, .tsv, etc.)

for (f in files){
  cat("Loading ", paste(root_dir, f, sep = ''), '\n')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F, sep = '\t')
  #tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  #tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
  cat("Writing: ", paste(root_dir, gsub(".tsv", ".fasta", f), sep = ''), '\n')
  write_fasta(tmp, "sequence", "sequence_id", paste(root_dir, gsub(".tsv", ".fasta", f), sep = ''))
}



