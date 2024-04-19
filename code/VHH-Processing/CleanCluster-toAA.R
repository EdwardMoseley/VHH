## ###
## CleanCluster-toAA
## First, and final, step in generating 
## clonal clusters for expert analysis.
## ###

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
args <- commandArgs(trailingOnly=TRUE)

file_name <- args[1]

## Amino Acid-only database

# Input data: "20_VecF2_R3_S20-pooled-families.csv"
# Output Data: "translated_dbs/20_VecF2_R3_S20-pooled-families_VHH_database.csv"

# store directory info
root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/"
output_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"

# Elevator music
cat("Loading: ", paste(root_path, file_name, sep = ''), '\n', sep = '')
# Load the file
tmp <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
# Capture source file
tmp$source_file <- rep(gsub("-rule-pooled-families.csv","", file_name), each = nrow(tmp))
# Translate - VHH
tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
# Trim sequences
tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
# Drop unproductive sequences
all_seqs <- nrow(tmp)
# Keep productive only
tmp <- unproductive(tmp)
filtered_seqs <- nrow(tmp)
# Elevator music
cat(file_name, '\n')
cat(round(100 - (filtered_seqs/all_seqs*100), 2), "% of sequences unproductive.\n", sep = '')
# Keep only necessary information
tmp <- tmp[,c("sequence_id","VHH","junction_aa","duplicate_count", "clone_id", "source_file")]
# Hold another frame for collapsing VHHs
tmp_frame <- tmp[ ,c("VHH", "duplicate_count")]
# Drop duplicate count, as we will add it back later
tmp$duplicate_count <- NULL
# Aggregate
tmp_frame <- aggregate(duplicate_count ~ VHH, data = tmp_frame, FUN = sum)
# Merge back to other data
tmp <- merge(tmp_frame, tmp, by = "VHH")
# Drop duplicates
tmp <- tmp[!duplicated(tmp$VHH),]
# Rename file
file_name <- gsub("-pooled-families.csv", "-pooled-families_VHH_database.csv", file_name)
# Write
write.csv(tmp, file = paste(output_dir, file_name, sep = ''), row.names = F)
# Elevator Music
cat("Finished writing: ", file_name, '\n', sep = '')

# Make system call for next script
#system(paste("Rscript CleanCluster-MMseqs2.R ", file_name, sep = ''))

## Route the traffic out
if (grepl("-pooled-families_VHH_database.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-MMseqs2.R", file_name), '\n')
  system(paste("Rscript CleanCluster-MMseqs2.R", file_name))
} else if (grepl("later file name for final pooling...", file_name)){
  cat("ERROR in [CleanCluster-toAA.R]: Wrong file name?")
}