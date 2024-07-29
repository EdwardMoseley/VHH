## VHHQueryBack
## Relies on MMseqs2
# https://github.com/soedinglab/MMseqs2

### Strategy ###
## First, we convert .tsv files to .fasta files after translating VHHs
## the key is the sequence_id, and the sequence is the VHH
## These are saved in the target_sequences directory

## Note that data and queries are already VHHs in fastas

# Utilities
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

# Path to MMseqs2
mmseqs_path <- "mmseqs "

# Root directory
root_dir <- "/cluster/tufts/cowenlab/emosel01/data/"
# These directories are just children of the root
query_dir <- "/cluster/tufts/cowenlab/emosel01/data/query_sequences/"
output_dir <- "/cluster/tufts/cowenlab/emosel01/data/tmp/"
data_dir <- "/cluster/tufts/cowenlab/emosel01/data/query_back/"
# Results file
results_file <- "XMY-H5_query_results29Jul24.csv"

files <- list.files(data_dir)
print(files)

# Load Query File
query_file <- list.files(query_dir)[1]
cat("Working on: ", query_file)

# Create data frame for results
results <- data.frame()

# Iterate over query databases
for (f in files){
  cat("Working on: ", f, '\n')
  #Build command
  command <- paste(mmseqs_path, "easy-search",
                   # Path to queries
                   paste0(query_dir, "queries.fasta"),
                   # Path to database
                   paste0(data_dir, f),
                   # Results ard held in
                   paste0(output_dir, "alnResult.m8"),
                   output_dir,
                   ## Options
                   " --search-type 3 --format-mode 4 -v 2", sep = ' ')
  
  cat("Calling command:\n")
  cat(command, '\n')
  system(command)
  
  # Load res
  tmp_frame <- read.csv(paste0(output_dir, "alnResult.m8"), header = T, stringsAsFactors = F, sep = '\t')
  # Change names quickly
  colnames(tmp_frame)[1] <- "sequence_id"
  # Order to pull top hits
  tmp_frame <- tmp_frame[order(tmp_frame$fident, decreasing = T),]
  tmp_frame$file <- rep(f, times = nrow(tmp_frame))
  # row bind and keep it moving
  results <- rbind(results, tmp_frame)
}

cat("Got results and writing...\n")
cat(paste0(output_dir, results_file,"\n"))

write.csv(results, file = paste0(output_dir, results_file), row.names = F)
