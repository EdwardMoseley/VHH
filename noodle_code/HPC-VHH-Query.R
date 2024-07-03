# VHH-Query
## Relies on MMseqs2
# https://github.com/soedinglab/MMseqs2

# Utilities
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

# Path to MMseqs2
mmseqs_path <- "mmseqs "

# Root directory
root_dir <- "/cluster/tufts/cowenlab/emosel01/data/"
# These directories are just children of the root
query_dir <- "/cluster/tufts/cowenlab/emosel01/data/queries/"
output_dir <- "/cluster/tufts/cowenlab/emosel01/data/tmp/"
data_dir <- "/cluster/tufts/cowenlab/emosel01/data/Assembled/"

files <- list.files(data_dir)[grepl("_db-pass.tsv", list.files(data_dir))]
files <- files[!grepl("demand",files)]

# Load Query File
query_file <- list.files(query_dir)[1]

# Queries
queries <- read.csv(paste0(query_dir, query_file), header = T, stringsAsFactors = F)

# Write query to fasta
write_fasta(queries, "VHH", "query_seq", paste0(output_dir, "queries.fasta"))





# Perform Queries on FASTA files


results <- data.frame()

for (f in files){
  #f <- "10_S19_db-pass.tsv"
  # Load data from one spot
  tmp <- read.csv(paste0(data_dir, f), header = T, stringsAsFactors = F, sep = '\t')
  # Translate VHHs
  tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  # Trim the beginning
  tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
  # Write it as a .fasta
  write_fasta(tmp, "VHH", "sequence_id", paste0(output_dir, "database.fasta"))
  # Call mmseqs2
  system(paste(mmseqs_path, "easy-search",
               "/Users/emoseley/VHH/code/VHH-Query/tmp/queries.fasta",
               "/Users/emoseley/VHH/code/VHH-Query/tmp/database.fasta",
               "/Users/emoseley/VHH/code/VHH-Query/tmp/alnResult.m8",
               "/Users/emoseley/VHH/code/VHH-Query/tmp/",
               " --search-type 3 --format-mode 4", sep = ' '))
  
  # Load res
  tmp_frame <- read.csv("~/VHH/code/VHH-Query/tmp/alnResult.m8", header = T, stringsAsFactors = F, sep = '\t')
  # Change names quickly
  colnames(tmp_frame)[2] <- "sequence_id"
  # Order to pull top hits
  tmp_frame <- tmp_frame[order(tmp_frame$fident, decreasing = T),]
  # Pull only top hit
  tmp_frame <- tmp_frame[!duplicated(tmp_frame$query),]
  # Merge to get VHH
  tmp_frame <- merge(tmp_frame, tmp[,c("sequence_id", "VHH", "junction_aa")], by = "sequence_id")
  # Change names again
  colnames(queries)[1] <- "query"
  tmp_frame <- merge(tmp_frame, queries, by = "query", suffixes = c("", "_query"), all.x = T)
  tmp_frame$database <- rep(f, times = nrow(tmp_frame))
  results <- rbind(results, tmp_frame)
}


# Trim junction from database seqs
results$junction_aa <- trim_CDR3(results$junction_aa)

# order
results <- results[order(results$fident, decreasing = T),]


write.csv("~/Downloads/Query_Results.csv", row.names = F)