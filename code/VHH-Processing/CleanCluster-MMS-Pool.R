## ### 
## Pooling CDR3's by MMseqs2 Clustering
## ###

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
# args <- commandArgs(trailingOnly=TRUE)

output_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/tmp/"

args <- c("GS_sample_min_threshold_clone-pass.csv")

# test if there is at least one argument: if not, return an error
if (grepl("-Pooled-Cleaning-Reclustering-Initial-Results.csv", args[1])) {
  file_name <- args[1]
  cat("Final Processing of: ", file_name, '\n', sep = '')
  # Define root_path
  root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"
  # Load Data
  tmp <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
} else if (grepl("_min_threshold_clone-pass.csv", args[1])){
  file_name <- args[1]
  cat("Initial Processing of: ", file_name, '\n', sep = '')
  # Define root_path
  root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/"
  # Load data
  tmp <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
  # Trim CDR3s
  tmp$junction_aa <- trim_CDR3(tmp$junction_aa)
} else{
  stop("Please provide the file in [prefix]_min_threshold_clone-pass.csv format OR [prefix]-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv format.", call.=FALSE)
}

## Check if there is a sequence_id asssocaited with the file. If not, add an arbitrary one.
## Note: This ID will likely derive from the Global Set
if (!("sequence_id" %in% colnames(tmp))){
  cat("Adding sequence_id\n")
  tmp$sequence_id <- paste0("GS-", 1:nrow(tmp))
}

if (!("VHH" %in% colnames(tmp))){
  cat("Translating VHH's\n")
  cat("Starting with ", nrow(tmp), " rows", '\n', sep = '')
  tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  cat("Starting with ", length(unique(tmp$VHH)), " Unique VHHs", '\n', sep = '')
  cat("Starting with ", sum(tmp$duplicate_count), " overall VHHs (including duplicates)", '\n', sep = '')
  tmp <- unproductive(tmp)
  cat('...\n')
  cat("Working with ", nrow(tmp), " rows after removing nonproductive VHHs ", '\n', sep = '')
  cat("Working with ", length(unique(tmp$VHH)), " Unique VHHs after removing nonproductive VHHs", '\n', sep = '')
  cat("Working with ", sum(tmp$duplicate_count), " overall VHHs (including duplicates) after removing nonproductive VHHs", '\n', sep = '')
}


## ###
## Begin with <6 CDR3s
## ### 60% H

tmp_frame <- tmp[nchar(tmp$junction_aa) < 6,]

cat("There are ", nrow(tmp_frame), " rows with < 6 amino acids in the CDR3.\n", sep = '')
cat("There are ", length(unique(tmp_frame$VHH)), " unique VHHs with < 6 amino acids in the CDR3.\n", sep = '')
cat("There are ", sum(tmp_frame$duplicate_count), " sequences (including duplicates) with < 6 amino acids in the CDR3.\n", sep = '')

homology <- 0.9

## ###
## Note error: 
# No k-mer could be extracted for the database /Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/tmp/tmp/7182412242506947631/clu_tmp/4287131758479141581/input_step_redundancy.
# Maybe the sequences length is less than 14 residues.

## We're going to fix this with Glycine

tmp_frame$junction_aa <- paste0("GGGGGGGGGGAAAAAAAAAAWWWWWWWWWW", tmp_frame$junction_aa, "GGGGGGGGGGAAAAAAAAAAWWWWWWWWWW")
# tmp_frame$junction_aa <- paste0("XXXXXXXXXXX", tmp_frame$junction_aa, "XXXXXXXXXXX")


## ###

# write fasta for mmseqs2
write_fasta(tmp_frame, "junction_aa", "sequence_id", paste(output_dir, "tmp.fasta", sep = ''))

# Prepare mmseqs command
command <- paste("/usr/local/Cellar/mmseqs2/15-6f452/bin/mmseqs ", 
                 "easy-cluster -v 2 --min-seq-id", paste0(homology, " "), 
                 paste(output_dir, "tmp.fasta", sep = ''), 
                 paste(output_dir, "tmp-H95", sep = ''),
                 paste(output_dir, "tmp", sep = ''))


## Call mmseqs to do its job
system(command)

## After clustering, we're going to merge clones by CDR3

dat <- read.csv(paste(output_dir,"tmp-H95_cluster.tsv", sep = ''), header = F, stringsAsFactors = F, sep = '\t')
# Assign Variable Names
colnames(dat) <- c("centroid", "sequence_id")

## dat is the clustered data, tmp is the entire data set
dat <- merge_db_to_centroids(dat, tmp)

## Now, we want to iterate over each centroid and assign the largest most common clone id to the whole crew
cat("There are ", length(unique(dat$clone_id)), " Unique Hierarchical Clone IDs before cleaning\n", sep = '')
cat("There are ", length(unique(dat$centroid_sequence_id)), " Unique Centroid (subcluster) IDs before cleaning\n", sep = '')
cat("There are ", length(unique(dat$junction_aa)), " Unique CDR3s before cleaning\n", sep = '')

















































# Input: "20_VecF2_R3_S20-rule-pooled-families_VHH_database.csv"
# Output: "20_VecF2_R3_S20-rule-pooled-families_subclustered_database.csv"

# store directory info
#root_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"
#output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/tmp/"


## ###
## MMseqs Clustering
## Our goal is to cluster with different homology demends
## depending on the length of the CDR3
## < 6, 6 - 10, 10+
## MMseqs2 will return a centroid, but we will want to pull clonal statistics
## ###

## Strategy:
## For each |CDR3|, perform clustering by homology (100%, 80%, 60%)

## ###
## Begin with >10 CDR3s
## ### 60% H

tmp_frame <- tmp[nchar(tmp$junction_aa) > 10,]

homology <- 0.60

# write fasta for mmseqs2
write_fasta(tmp_frame, "sequence", "sequence_id", paste(output_dir, "tmp.fasta", sep = ''))

# Prepare mmseqs command
command <- paste("/usr/local/Cellar/mmseqs2/15-6f452/bin/mmseqs ", 
                 "easy-cluster -v 2 --min-seq-id ", paste0(homology, " "), 
                 paste(output_dir, "tmp.fasta", sep = ''), 
                 paste(output_dir, "tmp-H95", sep = ''),
                 paste(output_dir, "tmp", sep = ''))


## Call mmseqs to do its job
# system(command)

# Load results
dat <- read.csv(paste(output_dir,"tmp-H95_cluster.tsv", sep = ''), header = F, stringsAsFactors = F, sep = '\t')
# Assign Variable Names
colnames(dat) <- c("centroid", "sequence_id")

dat <- merge_db_to_centroids(dat, tmp)

## ### ## ### ## ###







## ###
## Begin with >10 CDR3s
## ### 60% H

tmp_frame <- tmp[nchar(tmp$junction_aa) < 6,]

homology <- 1.00

# write fasta for mmseqs2
write_fasta(tmp_frame, "sequence", "sequence_id", paste(output_dir, "tmp.fasta", sep = ''))

# Prepare mmseqs command
command <- paste("/usr/local/Cellar/mmseqs2/15-6f452/bin/mmseqs ", 
                 "easy-cluster -v 2 --min-seq-id ", paste0(homology, " "), 
                 paste(output_dir, "tmp.fasta", sep = ''), 
                 paste(output_dir, "tmp-H95", sep = ''),
                 paste(output_dir, "tmp", sep = ''))


## Call mmseqs to do its job
system(command)



# Load results
dat <- read.csv(paste(output_dir,"tmp-H95_cluster.tsv", sep = ''), header = F, stringsAsFactors = F, sep = '\t')
# Assign Variable Names
colnames(dat) <- c("centroid", "sequence_id")

# dat <- merge_db_to_centroids(dat, tmp)

#function(db, tmp){
  db <- merge(dat, tmp, by = "sequence_id")
  # Pull out data and prepare to change column names for ease
  tmp <- dat[,c("sequence_id", "centroid")]
  # Update column names of tmp
  colnames(tmp) <- c("tmp_sequence_id", "sequence_id")
  # Join for Centroid sequences
  tmp <- merge(tmp, dat[,c("sequence_id", "VHH")], by = "sequence_id")
  # Correct column names
  colnames(tmp) <- c("centroid_sequence_id", "sequence_id", "centroid")
  # Remove "centroid" from dat
  dat$centroid <- NULL
  # Final merge
  dat <- merge(dat, tmp, by = "sequence_id")
#  return(dat)
#}








































file_name <- gsub("VHH_database.csv", "subclustered_database.csv", file_name)

cat("Writing: ", paste(root_dir, file_name, sep = ''), '\n', sep = '')

write.csv(results, file = paste(root_dir, file_name, sep = ''), row.names = F)

## Route the traffic out
if (grepl("_subclustered_database.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-Refining.R", file_name), '\n')
  system(paste("Rscript CleanCluster-Refining.R", file_name))
} else if (grepl("later file name for final pooling...", file_name)){
  cat("ERROR in [CleanCluster-toAA.R]: Wrong file name?")
  # change file name
  # write
}






