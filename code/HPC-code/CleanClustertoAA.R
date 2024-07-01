## ###
## CleanCluster-toAA
## Converts DNA to Amino Acid Sequences
## Removes unproductive sequences
## Trims vector sequence from VHHs (11 AAs)
## Writes a database of VHHs 
## Calls: CleanCluster-MMseqs2.R
## ###

## Load Utilities
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

## Draw file info from commandArgs
args <- commandArgs(trailingOnly=TRUE)

file_name <- args[1]

## Amino Acid-only database

# Input data: "20_VecF2_R3_S20-pooled.csv"
# Output Data: "translated_dbs/20_VecF2_R3_S20-pooled-families_VHH_database.csv"

# store directory info
root_path <- "/cluster/tufts/cowenlab/emosel01/data/forClustering/tmp/"
output_dir <- "/cluster/tufts/cowenlab/emosel01/data/forClustering/tmp/"

# Elevator music
cat("Loading: ", paste(root_path, file_name, sep = ''), '\n', sep = '')
# Load the file
tmp <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
# Capture source file
tmp$source_file <- rep(gsub("_pooled.csv","", file_name), each = nrow(tmp))
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

## Check if there are sequence id's, if not, it is the Global Set and they need be added
if (!("sequence_id" %in% colnames(tmp))){
  cat("No sequence id, so adding it.\n")
  tmp$sequence_id <- paste("GS-", 1:nrow(tmp), sep = '')
}


# Keep only necessary information
tmp <- tmp[,c("sequence_id", "sequence", "sequence_alignment", "VHH", "junction_aa", "junction", "duplicate_count", "clone_id", "source_file")]

## ###
## Keep an eye out for bugs here
## ###

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
file_name <- gsub("_pooled.csv", "_pooled_VHH_database.csv", file_name)
# Write
write.csv(tmp, file = paste(output_dir, file_name, sep = ''), row.names = F)
# Elevator Music
cat("Finished writing: ", paste(output_dir, file_name, sep = ''), '\n', sep = '')
