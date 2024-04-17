## ###
## Pooling with CleanClust
## ###


# Load utilities
source("~/VHH/lib/utilities.R")
suppressPackageStartupMessages(library(stringdist))
suppressPackageStartupMessages(library(shazam))

## ###
## Variables
## dat <- input data (clonal data from ChangeO)
##    - Updated throughout this process to clean clone_ids
## results <- Clone_id, CDR3 statistics for use in cluster pooling
## df <- temporary subset of "results" for pooling

## ###
## Prep Environment
## ###

root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/"

# file_name <- "18_VecF3_R4_S18_min_threshold_clone-pass.csv"
file_name <- "19_VecF1_R3_S19_min_threshold_clone-pass.csv"
# file_name <- "13_VecF1_R2_S13_min_threshold_clone-pass.csv"
# file_name <- "20_VecF2_R3_S20-pooled-families.csv"

dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)

# Trim CDR3s
dat$junction_aa <- trim_CDR3(dat$junction_aa)

## ###

## ###
## Begin Main
## ###

## Aggregate CDR3s and generate statistics

results <- clonal_statistics(dat)

## ###
## Merge IDENTICAL CDR3 familes
## ###

output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/tmp/"

easyclust_for_pooling <- function(results, output_dir, homology, short_seqs = FALSE){
  # write fasta for mmseqs2
  write_fasta(results, "median_CDR3", "clone_id", paste(output_dir, "tmp.fasta", sep = ''))
  
  if (short_seqs){
    # https://github.com/soedinglab/MMseqs2/issues/149
    # Prepare mmseqs command
    command <- paste("/usr/local/Cellar/mmseqs2/15-6f452/bin/mmseqs ", 
                     "easy-cluster --single-step-clustering -s 6 -v 3 --min-seq-id", homology, 
                     paste(output_dir, "tmp.fasta", sep = ''), 
                     paste(output_dir, "tmp-H100", sep = ''),
                     paste(output_dir, "tmp", sep = ''))
  } else {
    # Prepare mmseqs command
    command <- paste("/usr/local/Cellar/mmseqs2/15-6f452/bin/mmseqs ", 
                     "easy-cluster -v 2 --min-seq-id", homology, 
                     paste(output_dir, "tmp.fasta", sep = ''), 
                     paste(output_dir, "tmp-H100", sep = ''),
                     paste(output_dir, "tmp", sep = ''))
  }
  ## Call mmseqs to do its job
  system(command)
  # Load results
  tmp <- read.csv(paste(output_dir,"tmp-H100_cluster.tsv", sep = ''), header = F, stringsAsFactors = F, sep = '\t')
  # Assign Variable Names
  colnames(tmp) <- c("centroid", "clone_id")
  # First merge
  tmp <- merge(results, tmp, by = "clone_id")
  # Prep for second
  colnames(tmp)[which(colnames(tmp) == "clone_id")] <- "sequence_clone_id"
  colnames(tmp)[which(colnames(tmp) == "centroid")] <- "clone_id"
  # Second merge
  tmp <- merge(tmp, results, by = c("clone_id", "median_CDR3"), suffixes = c("", "_centroid"))
  # Clean
  tmp$sequence_clone_id <- as.numeric(tmp$sequence_clone_id)
  tmp$abundance <- as.numeric(tmp$abundance)
  tmp$abundance_centroid <- as.numeric(tmp$abundance_centroid)
  
  return(tmp)
}

# Perform clustering
tmp <- easyclust_for_pooling(results, output_dir, 1, FALSE)

# Order
tmp <- tmp[order(tmp$clone_id, tmp$sequence_clone_id, decreasing = T),]

# Keep only observations with duplicates
tmp <- tmp[duplicated(tmp$clone_id) | duplicated(tmp$clone_id, fromLast=TRUE),]

# Iterate over clone_ids
for (nm in unique(tmp$clone_id)){
  cat(nm, '\n')
  # Subset
  tmp_frame <- tmp[tmp$clone_id == nm,]
  # Order by abundance
  tmp_frame <- tmp_frame[order(tmp_frame$abundance, decreasing = T),]
  # Now we merge groups
  top_clone_id <- tmp_frame$clone_id[1]
  # Iterate over 2..length(tablet), adding other groups to the top group
  for (tmp_clone_id in tmp_frame$sequence_clone_id){
    cat(tmp_clone_id, '\n')
    dat[dat$clone_id == tmp_clone_id,]$clone_id <- rep(nm, 
                                                       each = nrow(dat[dat$clone_id == tmp_clone_id,]))
  }
}



## ###
## Now we pool: >10AA groups, >60% Homology
## ###

## After dupes (which accounts for <= 5 CDRs), we go to >10AA groups

results <- clonal_statistics(dat)

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 10,]

# Perform clustering
tmp <- easyclust_for_pooling(df, output_dir, 0.6)

# Order
tmp <- tmp[order(tmp$clone_id, tmp$sequence_clone_id, decreasing = T),]

# Keep only observations with duplicates
tmp <- tmp[duplicated(tmp$clone_id) | duplicated(tmp$clone_id, fromLast=TRUE),]

if (any(duplicated(tmp$clone_id))){
  # Iterate over clone_ids
  for (nm in unique(tmp$clone_id)){
    cat(nm, '\n')
    # Subset
    tmp_frame <- tmp[tmp$clone_id == nm,]
    # Order by abundance
    tmp_frame <- tmp_frame[order(tmp_frame$abundance, decreasing = T),]
    # Now we merge groups
    top_clone_id <- tmp_frame$clone_id[1]
    # Iterate over 2..length(tablet), adding other groups to the top group
    for (tmp_clone_id in tmp_frame$sequence_clone_id){
      cat(tmp_clone_id, '\n')
      dat[dat$clone_id == tmp_clone_id,]$clone_id <- rep(nm, 
                                                         each = nrow(dat[dat$clone_id == tmp_clone_id,]))
    }
  }
} else {
  print("No pooling required!")
}


## ###
## Now we pool: 6<x<10 groups, >80% Homology
## ###

# Recalculate stats
results <- clonal_statistics(dat)

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 6 & nchar(results$median_CDR3) < 10,]

# Perform clustering
tmp <- easyclust_for_pooling(df, output_dir, 0.7, TRUE)

# Order
tmp <- tmp[order(tmp$clone_id, tmp$sequence_clone_id, decreasing = T),]

# Keep only observations with duplicates
tmp <- tmp[duplicated(tmp$clone_id) | duplicated(tmp$clone_id, fromLast=TRUE),]

if (any(duplicated(tmp$clone_id))){
  # Iterate over clone_ids
  for (nm in unique(tmp$clone_id)){
    cat(nm, '\n')
    # Subset
    tmp_frame <- tmp[tmp$clone_id == nm,]
    # Order by abundance
    tmp_frame <- tmp_frame[order(tmp_frame$abundance, decreasing = T),]
    # Now we merge groups
    top_clone_id <- tmp_frame$clone_id[1]
    # Iterate over 2..length(tablet), adding other groups to the top group
    for (tmp_clone_id in tmp_frame$sequence_clone_id){
      cat(tmp_clone_id, '\n')
      dat[dat$clone_id == tmp_clone_id,]$clone_id <- rep(nm, 
                                                         each = nrow(dat[dat$clone_id == tmp_clone_id,]))
    }
  }
} else {
  print("No pooling required!")
}












