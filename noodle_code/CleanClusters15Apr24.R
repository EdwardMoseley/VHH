## ###
## CleanCluster-Scratchpad
## 15Apr24
## ###

## Trick for merging families:

# Idea: merging families with MAV CDR3’s of a given % Homology cannot change the MAV of the largest family

## Utilities

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
# file_name <- "19_VecF1_R3_S19_min_threshold_clone-pass.csv"
# file_name <- "13_VecF1_R2_S13_min_threshold_clone-pass.csv"
file_name <- "20_VecF2_R3_S20-pooled-families.csv"

dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)

# Trim CDR3s
dat$junction_aa <- trim_CDR3(dat$junction_aa)

## ###

## ###
## Begin Main
## ###

# Hold results in results
results <- data.frame()
# Iterate over clone ids
for (nm in unique(dat$clone_id)){
  tmp_frame <- dat[dat$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  abundance <- sum(tmp_frame$duplicate_count)

  # Keep the integer as the cluster_id
  clone_id <- nm
  
  results <- rbind(results, cbind(clone_id, median_CDR3, abundance))
}

## ###
## Handle exact matches first (this will cover cases where |CDR3| < 6 AA's)
## ###

## Identify groups with the same median CDR3
df <- results[duplicated(results$median_CDR3) | duplicated(results$median_CDR3, fromLast=TRUE),]

for (median_CDR3 in unique(df$median_CDR3)){
  # Subset from the statistics database
  tmp <- df[df$median_CDR3 == median_CDR3,]
  # Subset from the clonal database
  tmp_frame <- dat[(dat$clone_id %in% tmp$clone_id),]
  # Tabulate clone_ids to find which family is larger
  tablet <- table(tmp_frame$clone_id)
  # Order
  tablet <- tablet[order(tablet, decreasing = T)]
  # Iterate over 2..length(tablet), adding other groups to the top group
  for (nm in as.numeric(names(tablet)[2:length(tablet)])){
    # Subset the old clone_ids and replace with the largest group to which they could belong ([1])
    dat[dat$clone_id == nm,]$clone_id <- rep(as.numeric(names(tablet[1])), each = nrow(dat[dat$clone_id == nm,]))
  }
}

## ###
## Double-Check!
## ###

# Hold results in results
results <- data.frame()
# Iterate over clone ids
for (nm in unique(dat$clone_id)){
  tmp_frame <- dat[dat$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  abundance <- sum(tmp_frame$duplicate_count)
  
  # Keep the integer as the cluster_id
  clone_id <- nm
  
  results <- rbind(results, cbind(clone_id, median_CDR3, abundance))
}

## Identify groups with the same median CDR3
df <- results[duplicated(results$median_CDR3) | duplicated(results$median_CDR3, fromLast=TRUE),]
if (nrow(df) == 0) print("No duplicates left!")

## ###
## Now we pool: >10AA groups, >60% Homology
## ###

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 10,]

# Apply function to each median_CDR3 value
tmp <- lapply(df$median_CDR3, find_nearest_neighbor, df)

# Extract nearest neighbor and distance from results
df$nearest_neighbor <- sapply(tmp, function(x) x$nearest_neighbor)
df$distance <- sapply(tmp, function(x) x$distance)

# Before pooling, we demand all CDR3s have identical lengths
df <- df[nchar(df$median_CDR3) == nchar(df$nearest_neighbor),]

# Calculate homology in terms of percents
df$homology <- round(100 - df$distance/nchar(df$nearest_neighbor)*100,2)

# Drop nearest-neighbors that don't make cutoff
df <- df[df$homology >= 60,]

df$abundance <- as.numeric(df$abundance)

# Now let's grab nearest_neighbor CDR3s
tmp <- results[,c("clone_id", "median_CDR3", "abundance")]
colnames(tmp) <- c("nearest_clone_id", "nearest_neighbor", "nearest_abundance")

# Merge together
df <- merge(df, tmp, by = c("nearest_neighbor"))
# Order
df <- df[order(df$abundance, decreasing = T),]

## ###
## Pooling subroutine
## ###

for (median_CDR3 in unique(df$median_CDR3)){
  # Subset from the statistics database
  tmp <- df[df$median_CDR3 == median_CDR3,]
  # Subset from the clonal database
  dat[(dat$clone_id %in% tmp$nearest_clone_id),]$clone_id <- rep(df[df$median_CDR3 == median_CDR3,]$clone_id, 
                                                                 each = nrow(dat[(dat$clone_id %in% tmp$nearest_clone_id),]))
}


## ###
## Double-Check!
## ###

# Hold results in results
results <- data.frame()
# Iterate over clone ids
for (nm in unique(dat$clone_id)){
  tmp_frame <- dat[dat$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  abundance <- sum(tmp_frame$duplicate_count)
  
  # Keep the integer as the cluster_id
  clone_id <- nm
  
  results <- rbind(results, cbind(clone_id, median_CDR3, abundance))
}

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 10,]

# Apply function to each median_CDR3 value
tmp <- lapply(df$median_CDR3, find_nearest_neighbor, df)

# Extract nearest neighbor and distance from results
df$nearest_neighbor <- sapply(tmp, function(x) x$nearest_neighbor)
df$distance <- sapply(tmp, function(x) x$distance)

# Before pooling, we demand all CDR3s have identical lengths
df <- df[nchar(df$median_CDR3) == nchar(df$nearest_neighbor),]

# Calculate homology in terms of percents
df$homology <- round(100 - df$distance/nchar(df$nearest_neighbor)*100,2)

# Drop nearest-neighbors that don't make cutoff
df <- df[df$homology >= 60,]

if (nrow(df) > 0){
  df$abundance <- as.numeric(df$abundance)
  
  # Now let's grab nearest_neighbor CDR3s
  tmp <- results[,c("clone_id", "median_CDR3", "abundance")]
  colnames(tmp) <- c("nearest_clone_id", "nearest_neighbor", "nearest_abundance")
  
  # Merge together
  df <- merge(df, tmp, by = c("nearest_neighbor"))
  # Order
  df <- df[order(df$abundance, decreasing = T),]
} else {
  print("Looks good!")
}

## ###
## Now we pool: 6<x<10 groups, >80% Homology
## ###

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 6 & nchar(results$median_CDR3) < 10,]

# Apply function to each median_CDR3 value
tmp <- lapply(df$median_CDR3, find_nearest_neighbor, df)

# Extract nearest neighbor and distance from results
df$nearest_neighbor <- sapply(tmp, function(x) x$nearest_neighbor)
df$distance <- sapply(tmp, function(x) x$distance)

# Before pooling, we demand all CDR3s have identical lengths
df <- df[nchar(df$median_CDR3) == nchar(df$nearest_neighbor),]

# Calculate homology in terms of percents
df$homology <- round(100 - df$distance/nchar(df$nearest_neighbor)*100,2)

# Drop nearest-neighbors that don't make cutoff
df <- df[df$homology >= 80,]

df$abundance <- as.numeric(df$abundance)

# Now let's grab nearest_neighbor CDR3s
tmp <- results[,c("clone_id", "median_CDR3", "abundance")]
colnames(tmp) <- c("nearest_clone_id", "nearest_neighbor", "nearest_abundance")

# Merge together
df <- merge(df, tmp, by = c("nearest_neighbor"))
# Order
df <- df[order(df$abundance, decreasing = T),]

## ###
## Pooling subroutine
## ###

for (median_CDR3 in unique(df$median_CDR3)){
  # Subset from the statistics database
  tmp <- df[df$median_CDR3 == median_CDR3,]
  # Subset from the clonal database
  dat[(dat$clone_id %in% tmp$nearest_clone_id),]$clone_id <- rep(df[df$median_CDR3 == median_CDR3,]$clone_id, 
                                                                 each = nrow(dat[(dat$clone_id %in% tmp$nearest_clone_id),]))
}
## ###
## Double-Check!
## ###

# Hold results in results
results <- data.frame()
# Iterate over clone ids
for (nm in unique(dat$clone_id)){
  tmp_frame <- dat[dat$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  abundance <- sum(tmp_frame$duplicate_count)
  
  # Keep the integer as the cluster_id
  clone_id <- nm
  
  results <- rbind(results, cbind(clone_id, median_CDR3, abundance))
}

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 6 & nchar(results$median_CDR3) < 10,]

# Apply function to each median_CDR3 value
tmp <- lapply(df$median_CDR3, find_nearest_neighbor, df)

# Extract nearest neighbor and distance from results
df$nearest_neighbor <- sapply(tmp, function(x) x$nearest_neighbor)
df$distance <- sapply(tmp, function(x) x$distance)

# Before pooling, we demand all CDR3s have identical lengths
df <- df[nchar(df$median_CDR3) == nchar(df$nearest_neighbor),]

# Calculate homology in terms of percents
df$homology <- round(100 - df$distance/nchar(df$nearest_neighbor)*100,2)

# Drop nearest-neighbors that don't make cutoff
df <- df[df$homology >= 80,]

if (nrow(df) > 0){
  df$abundance <- as.numeric(df$abundance)
  
  # Now let's grab nearest_neighbor CDR3s
  tmp <- results[,c("clone_id", "median_CDR3", "abundance")]
  colnames(tmp) <- c("nearest_clone_id", "nearest_neighbor", "nearest_abundance")
  
  # Merge together
  df <- merge(df, tmp, by = c("nearest_neighbor"))
  # Order
  df <- df[order(df$abundance, decreasing = T),]
} else {
  print("Looks good!")
}

file_name <- gsub("_min_threshold_clone-pass.csv", "-rule-pooled-families.csv", file_name)

write.csv(dat, file = paste(root_path, file_name, sep = ''), row.names = F)

rm(df)

## ###
## End [VHH]-Clean-and-Outlier-Repooling
## ###


## ###
## BEGIN Collapse-DBs-to-AA-Seqs.R
## ###

## Amino Acid-only database

# Input data: "20_VecF2_R3_S20-pooled-families.csv"
# Output Data: "translated_dbs/20_VecF2_R3_S20-pooled-families_VHH_database.csv"

# store directory info
output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"

# Iterate over files
#for (f in files){
# Load the file
tmp <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
# Capture source file
tmp$source_file <- rep(gsub("-rule-pooled-families.csv","_min_threshold_clone-pass.csv", file_name), each = nrow(tmp))
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
#}

## ###
## END Collapse-DBs-to-AA-Seqs.R
## ###

## ###
## BEGIN AA-DB-mmseqs2-Clustering.R
## ###

# AA-DB-mmseqs2-Clustering

# Input: "20_VecF2_R3_S20-rule-pooled-families_VHH_database.csv"
# Output: "20_VecF2_R3_S20-rule-pooled-families_subclustered_database.csv"

# store directory info
root_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"
output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/tmp/"

# Load the file
tmp <- read.csv(paste(root_dir, file_name, sep = ''), header = T, stringsAsFactors = F)
# Create a data frame for holding results
results <- data.frame()
# Iterate over Hierarchically Clustered clone_ids, re-clustering
for (nm in unique(tmp$clone_id)){
  # Subset observations
  tmp_frame <- tmp[tmp$clone_id == nm,]
  # write fasta for mmseqs2
  write_fasta(tmp_frame, "VHH", "sequence_id", paste(output_dir, "tmp.fasta", sep = ''))
  # Prepare mmseqs command
  command <- paste("/usr/local/Cellar/mmseqs2/15-6f452/bin/mmseqs ", 
                   "easy-cluster -v 2 --min-seq-id 0.90 ", 
                   paste(output_dir, "tmp.fasta", sep = ''), 
                   paste(output_dir, "tmp-H95", sep = ''),
                   paste(output_dir, "tmp", sep = ''))
  ## Call mmseqs to do its job
  system(command)
  # Load results
  dat <- read.csv(paste(output_dir,"tmp-H95_cluster.tsv", sep = ''), header = F, stringsAsFactors = F, sep = '\t')
  # Assign Variable Names
  colnames(dat) <- c("centroid", "sequence_id")
  # Create a centroid id
  dat$centroid_id <- as.numeric(factor(dat$centroid))
  # We're going to use the sequence ids as keys, so we want to clean them
  tmp$sequence_id <- clean_seq_ids(tmp$sequence_id)
  dat$centroid <- clean_seq_ids(dat$centroid)
  dat$sequence_id <- clean_seq_ids(dat$sequence_id)
  # merge results
  dat <- merge_db_to_centroids(dat, tmp)
  # Combine results and move on
  results <- rbind(results, dat)
}

file_name <- gsub("VHH_database.csv", "subclustered_database-H90.csv", file_name)

write.csv(results, file = paste(root_dir, file_name, sep = ''), row.names = F)

rm(tmp, tmp_frame)

## ###
## END AA-DB-mmseqs2-Clustering
## ###

## ###
## BEGIN "[VHH]-Cluster-Refining"
## ###

#Input: 20_VecF2_R3_S20-pooled-families_subclustered_database.csv
#Output:
#(1) -> S19-pooled-clean-database.csv
#(2) -> S19-pooled-outlier-database.csv

results_file <- gsub("-pooled-families_subclustered_database-H90.csv", "-pooled-clean-database.csv", file_name)
outliers_file <- gsub("-pooled-families_subclustered_database-H90.csv", "-pooled-outlier-database.csv", file_name)

## Here, we want to:
#1. Remove outlier subfamilies
# * 0-5AA CDR3 < 100% Homology to Herarchically-Clustered VHH
# * 6-10AA CDR3 < 80% Homology to Herarchically-Clustered VHH
# * >10AA CDR3 < 60% Homology to Herarchically-Clustered VHH
#2. Place them into a "Outlier Database"
#3. Re-define Hiearchical Clones
#4. Sub-cluster each clone with MMseqs2

#dat <- read.csv("~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/19_VecF1_R3_S19_subclustered_database.csv", 
#                header = T, stringsAsFactors = F)

root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"

# dat <- results
dat <- read.csv(paste(root_dir, file_name, sep = ''), header = T, stringsAsFactors = F)


### Family Statistics:
# Add duplicate info here, also?

# Hold results in results
results <- data.frame()
# Iterate over clone ids
for (nm in unique(dat$clone_id)){
  tmp_frame <- dat[dat$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  # Distance from median
  tmp_frame$family_median_VHH_homology <- 100 - (as.vector(adist(median_VHH, tmp_frame$VHH))/nchar(median_VHH))*100
  tmp_frame$family_median_CDR3_homology <- 100 - (as.vector(adist(median_CDR3, tmp_frame$junction_aa))/nchar(median_CDR3))*100
  # Stack
  results <- rbind(results, tmp_frame)
}

### Sub-Family Statistics:

### Note: Cleaning can probably be added here, too?
dat <- results
results <- data.frame()
for (clone_id in unique(dat$clone_id)){
  # Subset clone id
  tmp <- dat[dat$clone_id == clone_id,]
  # For each centroid id
  for (nm in unique(tmp$centroid_id)){
    #cat("Centroid ID:", nm, '\n')
    tmp_frame <- tmp[tmp$centroid_id == nm,]
    # Unroll strings
    unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
    unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
    # Pull statistics
    median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
    median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
    # Distance from median
    tmp_frame$subfamily_median_VHH_homology <- 100 - (as.vector(adist(median_VHH, tmp_frame$VHH))/nchar(median_VHH))*100
    # Subfamiliy Median
    tmp_frame$subfamily_median_CDR3_homology <- 100 - (as.vector(adist(median_CDR3, tmp_frame$junction_aa))/nchar(median_CDR3))*100
    # Stack frames
    results <- rbind(results, tmp_frame)
  }
}

# Ordering
results <- results[order(results$family_median_VHH_homology, 
                         results$subfamily_median_VHH_homology, 
                         results$duplicate_count, 
                         decreasing = T),]

# Order statistic
results$order_index <- 1:nrow(results)

## Iterate Over families and centroids, culling members for further clustering

dat <- results

#1. Remove outlier subfamilies
#  * 0-5AA CDR3 < 100% Homology to Herarchically-Clustered VHH
#  * 6-10AA CDR3 < 80% Homology to Herarchically-Clustered VHH
#  * >10AA CDR3 < 60% Homology to Herarchically-Clustered VHH
#2. Place them into a "Outlier Database"
#3. Re-define Hiearchical Clones
#4. Sub-cluster each clone with MMseqs2

# Store results (keepers)
results <- data.frame()
# Store outliers (for re-clustering)
outliers <- data.frame()

# Iterate over unique clone ids
for (clone_id in unique(dat$clone_id)){
  cat("Clone:", clone_id, '\n')
  # Subset clone id
  tmp <- dat[dat$clone_id == clone_id,]
  
  # For each centroid id
  for (nm in unique(tmp$centroid_id)){
    cat("Subcluster ID:", nm, '\n')
    # Subset centroid id
    tmp_frame <- tmp[tmp$centroid_id == nm,]
    # Calculate quick stats for rule-based cleaning
    CDR3_length <- nchar(tmp_frame$junction_aa[1])
    avg_CDR3_percent_homology <- mean(tmp_frame$family_median_CDR3_homology)
    
    # Rule-based Cleanup here
    if (CDR3_length < 5 & avg_CDR3_percent_homology < 100){
      outliers <- rbind(outliers, tmp_frame)
    } else if (CDR3_length < 10 & avg_CDR3_percent_homology < 80){
      outliers <- rbind(outliers, tmp_frame)
    } else if (CDR3_length >= 10 & avg_CDR3_percent_homology < 60){
      outliers <- rbind(outliers, tmp_frame)
    } else {
      results <- rbind(results, tmp_frame)
    }
  }
  
}


## Check Results
cat("Percent VHH outliers: ", round(nrow(outliers)/nrow(dat), 2), '\n')
cat("Subfamily outliers:\n")
print(table(outliers$clone_id)[order(table(outliers$clone_id), decreasing = T)])


root_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"

results_file <- gsub("-pooled-families_subclustered_database-H90.csv", "-pooled-clean-database-H90.csv", file_name)
outliers_file <- gsub("-pooled-families_subclustered_database-H90.csv", "-pooled-outlier-database-H90.csv", file_name)

cat("Writing: ", paste(root_dir, results_file, sep = ''), '\n', sep = '')
write.csv(results, 
          file = paste(root_dir, results_file, sep = ''), 
          row.names = F)
cat("Writing: ", paste(root_dir, outliers_file, sep = ''), '\n', sep = '')
write.csv(outliers, 
          file = paste(root_dir, outliers_file, sep = ''), 
          row.names = F)

rm(tmp, tmp_frame)


## ###
## END "[VHH]-Cluster-Refining"
## ###

## ###
## BEGIN "[VHH]-Outlier-Re-Clustering"
## ###

#Inputs: 
#(1) -> 20_VecF2_R3_S20-pooled-families_subclustered_database.csv
#(2) -> 20_VecF2_R3_S20_min_threshold_clone-pass.csv"
#Output:
#(1) -> S20-pooled-outlier-reclustered-database.csv

# results_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-clean-database.csv", f)
# outliers_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-outlier-database.csv", f)

## Load Outliers

root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"

# Load outliers
dat <- read.csv(paste(root_dir, outliers_file, sep = ''), header = T, stringsAsFactors = F)

## Load Source Database
tmp <- read.csv(paste("~/VHH/data/Cluster-DBs/phylogenetic_dbs/", gsub("-rule-pooled-families_subclustered_database-H90.csv", "_min_threshold_clone-pass.csv", file_name), sep = ''),
                header = T, stringsAsFactors = F)
# Trim CDR3s
tmp$junction_aa <- trim_CDR3(tmp$junction_aa)

# Translate VHHs
tmp$VHH <- translate_DNA(tmp$sequence_alignment)
# Trim sequences
tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
# We need those VHHs from tmp which are represented in outlier file
dat <- tmp[tmp$VHH %in% dat$VHH,]

## Prepare for clonal inference

dat$clone_id <- NULL
dat$dist_nearest <- NULL

## Clonal Inference

cat("Calculating distToNearest...\n")

db <- distToNearest(dat, model = "ham", normalize = "len", 
                    vCallColumn = "v_call", 
                    jCallColumn = "j_call",
                    nproc = 4)


# Determine threshold (based on median)

thr <- summary(db$dist_nearest)[["Median"]]

new_thr <- thr + (thr*0.25)
print("25% Higher Threshold:")
print(round(new_thr, 2))

cat("Inferring Clones:\n")
system.time(
# Clonal assignment using hierarchical clustering
db <- scoper::hierarchicalClones(db, 
                                 threshold=new_thr, 
                                 v_call="v_call", 
                                 j_call = "j_call")
)

db <- db@db

## Sub-cluster with MMseqs2

output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/tmp/"

results <- data.frame()

# Iterate over Hierarchically Clustered clone_ids, re-clustering
for (nm in unique(db$clone_id)){
  # nm <- 8 (51 obs)
  # Subset observations
  tmp_frame <- db[db$clone_id == nm,]
  # If the family has more than one member, we will re-cluster
  if (nrow(tmp_frame) > 1){
    # write fasta for mmseqs2
    write_fasta(tmp_frame, "VHH", "sequence_id", paste(output_dir, "tmp.fasta", sep = ''))
    # Prepare mmseqs command
    command <- paste("/usr/local/Cellar/mmseqs2/15-6f452/bin/mmseqs ", 
                     "easy-cluster -v 2 --min-seq-id 0.90 ", 
                     paste(output_dir, "tmp.fasta", sep = ''), 
                     paste(output_dir, "tmp-H95", sep = ''),
                     paste(output_dir, "tmp", sep = ''))
    ## Call mmseqs to do its job
    system(command)
    # Load results
    dat <- read.csv(paste(output_dir,"tmp-H95_cluster.tsv", sep = ''), header = F, stringsAsFactors = F, sep = '\t')
    # Assign Variable Names
    colnames(dat) <- c("centroid", "sequence_id")
    # Create a centroid id
    dat$centroid_id <- as.numeric(factor(dat$centroid))
    # We're going to use the sequence ids as keys, so we want to clean them
    db$sequence_id <- clean_seq_ids(db$sequence_id)
    dat$centroid <- clean_seq_ids(dat$centroid)
    dat$sequence_id <- clean_seq_ids(dat$sequence_id)
    # merge results
    dat <- merge_db_to_centroids(dat, db)
  } else {
    # otherwise, we'll make up a centroid id and a sequence id, and
    tmp_frame$centroid_id <- tmp_frame$sequence_id[1]
    tmp_frame$centroid_sequence_id <- tmp_frame$sequence_id[1]
    tmp_frame$centroid <- tmp_frame$VHH[1]
    # Now we assign to dat for rbind
    dat <- tmp_frame
  }
  # Combine results and move on
  results <- rbind(results, dat)
}

results$centroid_id <- as.numeric(factor(results$centroid_id))

### Temporary Write

# Note:
# input_file <- "19_VecF1_R3_S19-pooled-outlier-database-H90.csv"
file_name <- gsub("_database-H90.csv", "-reclustered-database-H90.csv", file_name)

write.csv(results, 
          file = paste(root_dir, file_name, sep = ''), row.names = F)

## ###
## END "[VHH]-Outlier-Re-Clustering"
## ###

## ###
## BEGIN [VHH]-Outlier-DB-Statistics
## ###

## Load Data

root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"
clean_file <- results_file
outlier_file <- file_name

# Load cleaned database
clean_db <- read.csv(paste(root_dir, clean_file, sep = ''), header = T, stringsAsFactors = F)

# Load re-clustered outliers
outlier_db <- read.csv(paste(root_dir, outlier_file, sep = ''), header = T, stringsAsFactors = F)

### Family Statistics:
clean_db <- family_statistics(clean_db)
outlier_db <- family_statistics(outlier_db)

### Sub-Family Statistics:
clean_test_db <- subfamily_statistics(clean_db)
outlier_test_db <- subfamily_statistics(outlier_db)

## Harmonize Outlier DB

# Clean column names
outlier_test_db <- outlier_test_db[colnames(outlier_test_db) %in% colnames(clean_test_db)]

### Add a large number to outlier clone IDs so they are distinguishable
outlier_test_db$clone_id <- outlier_test_db$clone_id + 100000
outlier_test_db$centroid_id <- outlier_test_db$centroid_id + 100000

## Aggregate for correct duplicate counts
tmp <- outlier_test_db
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
# Assign to outlier db
outlier_test_db <- tmp
# Cleanup
rm(tmp, tmp_frame)

### Stack and Write

clean_test_db$group <- rep("clean", each = nrow(clean_test_db))
outlier_test_db$group <- rep("outliers", each = nrow(outlier_test_db))

# Add source file to outlier group
outlier_test_db$source_file <- rep(clean_test_db$source_file[1], each = nrow(outlier_test_db))

results <- rbind(clean_test_db, outlier_test_db)

# Ordering
results <- results[order(results$family_median_VHH_homology, 
                         results$subfamily_median_VHH_homology, 
                         results$duplicate_count, 
                         decreasing = T),]

# Order statistic
results$order_index <- 1:nrow(results)

root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"
output_file <- gsub("-pooled-clean-database-H90.csv", "-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv", clean_file)

write.csv(results, file = paste(root_dir, output_file, sep = ''), row.names = F)

## ###
## END [VHH]-Outlier-DB-Statistics
## ###

## ###
## BEGIN Final Pooling
## ###

#root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"

#clean_file <- "18_VecF3_R4_S18-rule-pooled-clean-database-H90.csv"

input_file <- gsub("-pooled-clean-database-H90.csv", "-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv", clean_file)

dat <- read.csv(paste(root_dir, input_file, sep = ''), header = T, stringsAsFactors = F)

## ###
## Begin Main
## ###

# Hold results in results
results <- data.frame()
# Iterate over clone ids
for (nm in unique(dat$clone_id)){
  tmp_frame <- dat[dat$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  abundance <- sum(tmp_frame$duplicate_count)
  
  # Keep the integer as the cluster_id
  clone_id <- nm
  
  results <- rbind(results, cbind(clone_id, median_CDR3, abundance))
}

## ###
## Handle exact matches first (this will cover cases where |CDR3| < 6 AA's)
## ###

## Identify groups with the same median CDR3
df <- results[duplicated(results$median_CDR3) | duplicated(results$median_CDR3, fromLast=TRUE),]

for (median_CDR3 in unique(df$median_CDR3)){
  # Subset from the statistics database
  tmp <- df[df$median_CDR3 == median_CDR3,]
  # Subset from the clonal database
  tmp_frame <- dat[(dat$clone_id %in% tmp$clone_id),]
  # Tabulate clone_ids to find which family is larger
  tablet <- table(tmp_frame$clone_id)
  # Order
  tablet <- tablet[order(tablet, decreasing = T)]
  # Iterate over 2..length(tablet), adding other groups to the top group
  for (nm in as.numeric(names(tablet)[2:length(tablet)])){
    # Subset the old clone_ids and replace with the largest group to which they could belong ([1])
    dat[dat$clone_id == nm,]$clone_id <- rep(as.numeric(names(tablet[1])), each = nrow(dat[dat$clone_id == nm,]))
  }
}

## ###
## Double-Check!
## ###

# Hold results in results
results <- data.frame()
# Iterate over clone ids
for (nm in unique(dat$clone_id)){
  tmp_frame <- dat[dat$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  abundance <- sum(tmp_frame$duplicate_count)
  
  # Keep the integer as the cluster_id
  clone_id <- nm
  
  results <- rbind(results, cbind(clone_id, median_CDR3, abundance))
}

## Identify groups with the same median CDR3
df <- results[duplicated(results$median_CDR3) | duplicated(results$median_CDR3, fromLast=TRUE),]
if (nrow(df) == 0) print("No duplicates left!")

## ###
## Now we pool: >10AA groups, >60% Homology
## ###

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 10,]

# Apply function to each median_CDR3 value
tmp <- lapply(df$median_CDR3, find_nearest_neighbor, df)

# Extract nearest neighbor and distance from results
df$nearest_neighbor <- sapply(tmp, function(x) x$nearest_neighbor)
df$distance <- sapply(tmp, function(x) x$distance)

# Before pooling, we demand all CDR3s have identical lengths
df <- df[nchar(df$median_CDR3) == nchar(df$nearest_neighbor),]

# Calculate homology in terms of percents
df$homology <- round(100 - df$distance/nchar(df$nearest_neighbor)*100,2)

# Drop nearest-neighbors that don't make cutoff
df <- df[df$homology >= 60,]

df$abundance <- as.numeric(df$abundance)

# Now let's grab nearest_neighbor CDR3s
tmp <- results[,c("clone_id", "median_CDR3", "abundance")]
colnames(tmp) <- c("nearest_clone_id", "nearest_neighbor", "nearest_abundance")

# Merge together
df <- merge(df, tmp, by = c("nearest_neighbor"))
# Order
df <- df[order(df$abundance, df$homology, decreasing = T),]

## ###
## Pooling subroutine
## ###

holder <- dat

dat <- holder

# Store for debugging convenience
# median_CDR3 <- "CAKDPGPGTWYSGSYFRDDDW"

# Goal: Merge smaller group IN TO larger group every time

# Iterate down the list (ordered by abundance)

for (median_CDR3 in unique(df$median_CDR3)){
  # Subset from the statistics database
  tmp <- df[df$median_CDR3 == median_CDR3,]
  # Subset from the clonal database
  dat[(dat$clone_id %in% tmp$nearest_clone_id),]$clone_id <- rep(df[df$median_CDR3 == median_CDR3,]$clone_id, 
                                                                 each = nrow(dat[(dat$clone_id %in% tmp$nearest_clone_id),]))
}

## ###
## Double-Check!
## ###

# Hold results in results
results <- data.frame()
# Iterate over clone ids
for (nm in unique(dat$clone_id)){
  tmp_frame <- dat[dat$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  abundance <- sum(tmp_frame$duplicate_count)
  
  # Keep the integer as the cluster_id
  clone_id <- nm
  
  results <- rbind(results, cbind(clone_id, median_CDR3, abundance))
}

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 10,]

# Apply function to each median_CDR3 value
tmp <- lapply(df$median_CDR3, find_nearest_neighbor, df)

# Extract nearest neighbor and distance from results
df$nearest_neighbor <- sapply(tmp, function(x) x$nearest_neighbor)
df$distance <- sapply(tmp, function(x) x$distance)

# Before pooling, we demand all CDR3s have identical lengths
df <- df[nchar(df$median_CDR3) == nchar(df$nearest_neighbor),]

# Calculate homology in terms of percents
df$homology <- round(100 - df$distance/nchar(df$nearest_neighbor)*100,2)

# Drop nearest-neighbors that don't make cutoff
df <- df[df$homology >= 60,]

if (nrow(df) > 0){
  df$abundance <- as.numeric(df$abundance)
  # Now let's grab nearest_neighbor CDR3s
  tmp <- results[,c("clone_id", "median_CDR3", "abundance")]
  colnames(tmp) <- c("nearest_clone_id", "nearest_neighbor", "nearest_abundance")
  # Merge together
  df <- merge(df, tmp, by = c("nearest_neighbor"))
  # Order
  df <- df[order(df$abundance, decreasing = T),]
} else {
  print("Looks good!")
}

## ###
## Now we pool: 6<x<10 groups, >80% Homology
## ###

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 6 & nchar(results$median_CDR3) < 10,]

# Apply function to each median_CDR3 value
tmp <- lapply(df$median_CDR3, find_nearest_neighbor, df)

# Extract nearest neighbor and distance from results
df$nearest_neighbor <- sapply(tmp, function(x) x$nearest_neighbor)
df$distance <- sapply(tmp, function(x) x$distance)

# Before pooling, we demand all CDR3s have identical lengths
df <- df[nchar(df$median_CDR3) == nchar(df$nearest_neighbor),]

# Calculate homology in terms of percents
df$homology <- round(100 - df$distance/nchar(df$nearest_neighbor)*100,2)

# Drop nearest-neighbors that don't make cutoff
df <- df[df$homology >= 80,]

df$abundance <- as.numeric(df$abundance)

# Now let's grab nearest_neighbor CDR3s
tmp <- results[,c("clone_id", "median_CDR3", "abundance")]
colnames(tmp) <- c("nearest_clone_id", "nearest_neighbor", "nearest_abundance")

# Merge together
df <- merge(df, tmp, by = c("nearest_neighbor"))
# Order
df <- df[order(df$abundance, decreasing = T),]

## ###
## Pooling subroutine
## ###

# Goal: Merge smaller group IN TO larger group every time

for (median_CDR3 in unique(df$median_CDR3)){
  # Subset from the statistics database
  tmp <- df[df$median_CDR3 == median_CDR3,]
  # Subset from the clonal database
  dat[(dat$clone_id %in% tmp$nearest_clone_id),]$clone_id <- rep(df[df$median_CDR3 == median_CDR3,]$clone_id, 
                                                                 each = nrow(dat[(dat$clone_id %in% tmp$nearest_clone_id),]))
}

## ###
## Double-Check!
## ###

# Hold results in results
results <- data.frame()
# Iterate over clone ids
for (nm in unique(dat$clone_id)){
  tmp_frame <- dat[dat$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  abundance <- sum(tmp_frame$duplicate_count)
  
  # Keep the integer as the cluster_id
  clone_id <- nm
  
  results <- rbind(results, cbind(clone_id, median_CDR3, abundance))
}

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 6 & nchar(results$median_CDR3) < 10,]

# Apply function to each median_CDR3 value
tmp <- lapply(df$median_CDR3, find_nearest_neighbor, df)

# Extract nearest neighbor and distance from results
df$nearest_neighbor <- sapply(tmp, function(x) x$nearest_neighbor)
df$distance <- sapply(tmp, function(x) x$distance)

# Before pooling, we demand all CDR3s have identical lengths
df <- df[nchar(df$median_CDR3) == nchar(df$nearest_neighbor),]

# Calculate homology in terms of percents
df$homology <- round(100 - df$distance/nchar(df$nearest_neighbor)*100,2)

# Drop nearest-neighbors that don't make cutoff
df <- df[df$homology >= 80,]

if (nrow(df) > 0){
  df$abundance <- as.numeric(df$abundance)
  
  # Now let's grab nearest_neighbor CDR3s
  tmp <- results[,c("clone_id", "median_CDR3", "abundance")]
  colnames(tmp) <- c("nearest_clone_id", "nearest_neighbor", "nearest_abundance")
  
  # Merge together
  df <- merge(df, tmp, by = c("nearest_neighbor"))
  # Order
  df <- df[order(df$abundance, decreasing = T),]
} else {
  print("Looks good!")
}

file_name <- gsub("-rule-pooled-families_subclustered-reclustered-database-H90.csv", "-H90-Cleaned10Apr24.csv", file_name)

write.csv(dat, file = paste(root_path, file_name, sep = ''), row.names = F)
