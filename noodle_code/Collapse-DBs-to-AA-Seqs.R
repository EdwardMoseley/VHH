## Amino Acid-only database

# Imports
source("~/VHH/lib/utilities.R")
suppressPackageStartupMessages(library(stringdist))

# store directory info
root_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/"
output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"

# Load up data
files <- list.files(root_dir)
files <- files[grepl(".csv", files)]

#files <- c("18_VecF3_R4_S18_min_threshold_clone-pass.csv", 
#           "19_VecF1_R3_S19_min_threshold_clone-pass.csv")

# Create a data frame for storing results
# results <- data.frame()

# Iterate over files
for (f in files){
  # Load the file
  tmp <- read.csv(paste(root_dir, f, sep = ''), header = T, stringsAsFactors = F)
  # Capture source file
  tmp$source_file <- rep(f, each = nrow(tmp))
  # Trim CDR3
  tmp$junction_aa <- trim_CDR3(tmp$junction_aa)
  # Translate - VHH
  tmp$VHH <- translate_DNA(tmp$sequence_alignment)
  # Trim sequences
  tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
  # Drop unproductive sequences
  all_seqs <- nrow(tmp)
  # Keep productive only
  tmp <- unproductive(tmp)
  filtered_seqs <- nrow(tmp)
  # Elevator music
  cat(f, '\n')
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
  f <- gsub("_min_threshold_clone-pass.csv", "_VHH_database.csv", f)
  # Write
  write.csv(tmp, file = paste(output_dir, f, sep = ''), row.names = F)
  #results <- rbind(tmp, results)
}



