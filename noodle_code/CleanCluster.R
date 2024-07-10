# CleanCluster Algorithm

## Quick utility

check_statistics <- function(dat, tmp){
  cat("Sequence ID Count identical: ", 
      length(unique(dat$sequence_id)) == length(unique(tmp$sequence_id)), 
      '\n')
  
  cat("Clone ID Count identical: ", 
      length(unique(dat$clone_id)) == length(unique(tmp$clone_id)), 
      '\n')
  
  cat("DNA Sequence Count identical: ", 
      length(unique(dat$sequence)) == length(unique(tmp$sequence)), 
      '\n')
  
  cat("Sequence Alignment Count identical: ", 
      length(unique(dat$sequence_alignment)) == length(unique(tmp$sequence_alignment)), 
      '\n')
}

print_statistics <- function(dat){
  cat("Sequence ID Count identical: ", 
      length(unique(dat$sequence_id)), 
      '\n')
  cat("Clone ID Count identical: ", 
      length(unique(dat$clone_id)), 
      '\n')
  cat("DNA Sequence Count identical: ", 
      length(unique(dat$sequence)), 
      '\n')
  cat("Sequence Alignment Count identical: ", 
      length(unique(dat$sequence_alignment)), 
      '\n')
}

## ###
## CleanCluster-Pool
## First, and final, step in generating 
## clonal clusters for expert analysis.
## ###

## First pass call:
# Rscript CleanCluster-Pool.R 21_VecF4_R3_S21_min_threshold_clone-pass.csv
# Rscript CleanCluster-Pool.R Assembled_9_S18_min_threshold_clone-pass.csv
# Percent VHH outliers:  3.42 
# Rscript CleanCluster-Pool.R Assembled_D4_S6_min_threshold_clone-pass.csv
## Last pass call:
# Rscript CleanCluster-Pool.R 19_VecF1_R3_S19-rule-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv 

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
# args <- commandArgs(trailingOnly=TRUE)

args <- c("20_VecF2_R3_S20_min_threshold_clone-pass.csv")

# Rscript CleanCluster-Pool.R Global_Clone_Experiment_master_set-D1-09-10-11_min_threshold_clone-pass.csv
# Rscript CleanCluster-Pool.R Global_Clone_Experiment_D1_S3_demand_two_min_threshold_clone-pass.csv
# Rscript CleanCluster-Pool.R Clones_Inferred_Assembled_9_S18_min_threshold_clone-pass.csv
# Rscript CleanCluster-Pool.R Clones_Inferred_Assembled_41_S16_min_threshold_clone-pass.csv

# test if there is at least one argument: if not, return an error
#if (grepl("-Pooled-Cleaning-Reclustering-Initial-Results.csv", args[1])) {
#  file_name <- args[1]
#  cat("Final Processing of: ", file_name, '\n', sep = '')
#  # Define root_path
#  root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"
#  # Load Data
#  dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
#} else 
  
if (grepl("_min_threshold_clone-pass.csv", args[1])){
  file_name <- args[1]
  cat("Initial Processing of: ", file_name, '\n', sep = '')
  # Define root_path
  root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/"
  # Load data
  dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
  # Trim CDR3s
  dat$junction_aa <- trim_CDR3(dat$junction_aa)
} else{
  stop("Please provide the file in [prefix]_min_threshold_clone-pass.csv format OR [prefix]-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv format.", call.=FALSE)
}

# Print clone counts
cat(file_name, " begins with ", length(unique(dat$clone_id)), " clonal families.\n", sep = '')

## Begin Cleaning

# Generate statistics
results <- clonal_statistics(dat)

## Handle Exact Matches first
### This will also handle <6AA CDR3s

# Identify duplicate median_CDR3s
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
  for (nm in names(tablet)[2:length(tablet)]){
    # Subset the old clone_ids and replace with the largest group to which they could belong ([1])
    dat[dat$clone_id == nm,]$clone_id <- rep(names(tablet[1]), each = nrow(dat[dat$clone_id == nm,]))
  }
}

## Double-Check

# Hold results in results
results <- clonal_statistics(dat)

## Identify groups with the same median CDR3
df <- results[duplicated(results$median_CDR3) | duplicated(results$median_CDR3, fromLast=TRUE),]
if (nrow(df) == 0) print("No duplicates left!")

## Update
cat(file_name, " has ", length(unique(dat$clone_id)), " clonal families after pooling families with identical mode CDR3s.\n", sep = '')

## Now we pool: >10AA groups, >60% Homology

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 10,]

tmp <- create_similarity_matrix(df)

# Drop nearest-neighbors that don't make cutoff
tmp <- tmp[tmp$value >= 60,]

# Pull in clone ids for merging
tmp <- merge(tmp, df, by = "median_CDR3")

# Order
tmp <- tmp[order(as.numeric(tmp$abundance), as.numeric(tmp$value), decreasing = T),]

## Pooling with % Homology

done_list <- vector()
# Begin merging routine
for (nm in unique(tmp$median_CDR3)){
  cat("Starting: ", nm, '\n', sep = '')
  # Pull observations we need
  tmp_frame <- tmp[(tmp$median_CDR3 == nm | tmp$variable == nm),]
  # Check if we've done this one yet
  if (!(nm %in% done_list)){
    cat("Merging: ", nm, '\n', sep = '')
    done_list <- c(done_list, tmp_frame$median_CDR3)
    # Begin adding lower-abundance groups to greater-abundance ones
    for (i in 2:(nrow(tmp_frame))){
      if (nrow(dat[dat$clone_id == tmp_frame$clone_id[i],]) >= 1){
        dat[dat$clone_id == tmp_frame$clone_id[i],]$clone_id <- tmp_frame$clone_id[1]
      } else {
        cat("No observations for Clone ID: ", tmp_frame$clone_id[i],'\n', sep = '')
      }
    }
  } else {
    cat("Passing: ", nm, '\n', sep = '')
  }
}

## ###
## Double Check
## ###

results <- clonal_statistics(dat)

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 10,]

tmp <- create_similarity_matrix(df)

# Drop nearest-neighbors that don't make cutoff
tmp <- tmp[tmp$value >= 60,]

if (nrow(tmp) == 0){
  print("No similar >10AA CDR3s left!")
} else {
  print("Did not exhaust >10AA CDR3s.")
}

cat(file_name, " has ", length(unique(dat$clone_id)), " clonal families after pooling 60% homologous families with >10AA CDR3s.\n", sep = '')

## Now we pool: 6<x<10 groups, >80% Homology

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 6 & nchar(results$median_CDR3) < 10,]

# Similarity matrix
tmp <- create_similarity_matrix(df)

# Drop nearest-neighbors that don't make cutoff
tmp <- tmp[tmp$value >= 80,]

# Pull in clone ids for merging
tmp <- merge(tmp, df, by = "median_CDR3")

# Order
tmp <- tmp[order(as.numeric(tmp$abundance), as.numeric(tmp$value), decreasing = T),]

done_list <- vector()
# Begin merging routine
for (nm in unique(tmp$median_CDR3)){
  cat("Starting: ", nm, '\n', sep = '')
  # Pull observations we need
  tmp_frame <- tmp[(tmp$median_CDR3 == nm | tmp$variable == nm),]
  # Check if we've done this one yet
  if (!(nm %in% done_list)){
    cat("Merging: ", nm, '\n', sep = '')
    done_list <- c(done_list, tmp_frame$median_CDR3)
    # Begin adding lower-abundance groups to greater-abundance ones
    for (i in 2:(nrow(tmp_frame))){
      if (nrow(dat[dat$clone_id == tmp_frame$clone_id[i],]) >= 1){
        dat[dat$clone_id == tmp_frame$clone_id[i],]$clone_id <- tmp_frame$clone_id[1]
      } else {
        cat("No observations for Clone ID: ", tmp_frame$clone_id[i],'\n', sep = '')
      }
    }
  } else {
    cat("Passing: ", nm, '\n', sep = '')
  }
}

## ###
## Double Check
## ###

# Hold results in results
results <- clonal_statistics(dat)

## Now we pool: 6<x<10 groups, >80% Homology

## Identify groups of interest (by AA size)
df <- results[nchar(results$median_CDR3) > 6 & nchar(results$median_CDR3) < 10,]

# Similarity matrix
tmp <- create_similarity_matrix(df)

# Drop nearest-neighbors that don't make cutoff
tmp <- tmp[tmp$value >= 80,]

if (nrow(tmp) == 0){
  print("No similar 6AA<x<10AA CDR3s left!")
} else {
  print("Have not exhausted 6AA<x<10AA CDR3s.")
}

cat(file_name, " has ", length(unique(dat$clone_id)), " clonal families after pooling 80% homologous families with 6AA<x<10AA CDR3s.\n", sep = '')

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## Important Variables at the end of Pooling:
# dat <- all the important data we want
#   ** This is CleanCluster-toAA.R 20_VecF2_R3_S20-rule-pooled-families.csv
# file_name <- the input file
# root_path <- folder the input file is in

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##


## Route the traffic out
if (grepl("_min_threshold_clone-pass.csv", file_name)){
  # change file name
  file_name <- gsub("_min_threshold_clone-pass.csv", "-rule-pooled-families.csv", file_name)
  # write.csv(dat, file = paste(root_path, file_name, sep = ''), row.names = F)
  # -rule-pooled-families.csv
  cat(paste("Rscript CleanCluster-toAA.R", file_name), '\n')
  # Pass the file to the next script in the pipeline
}

## Cleanup dfs (results are just summary stats)
rm(df, tmp, tmp_frame, results)

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##









































## *** *** *** *** *** *** *** *** *** *** *** *** *** ##
## *** *** *** ***  POOLING CHECK  *** *** *** *** *** ##
## *** *** *** *** *** *** *** *** *** *** *** *** *** ##

# tmp <- read.csv(paste0(root_path, file_name), header = T, stringsAsFactors = F)

# check_statistics(dat, tmp)

print_statistics(dat)

## *** *** *** *** *** *** *** *** *** *** *** *** *** ##





















## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##


## Begin toAA

## ###
## CleanCluster-toAA
## ###


# holder <- dat


# tmp <- holder

## Amino Acid-only database

# Input data: "20_VecF2_R3_S20-pooled-families.csv"
# Output Data: "translated_dbs/20_VecF2_R3_S20-pooled-families_VHH_database.csv" 

## No need to output data, we have a lot of ram
output_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"

# Assign source file information
tmp$source_file <- rep(gsub("-rule-pooled-families.csv", "",file_name), times = nrow(tmp))
# Translate - VHH
tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)


# Trim sequences (or not yet?)
## tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))

# Drop unproductive sequences
all_seqs <- nrow(tmp)
# Keep productive only
tmp <- unproductive(tmp)
filtered_seqs <- nrow(tmp)

# Elevator music
cat(file_name, '\n')
cat(round(100 - (filtered_seqs/all_seqs*100), 2), "% of sequences unproductive.\n", sep = '')

# Keep only necessary information
# tmp <- tmp[,c("sequence_id","VHH","junction_aa","duplicate_count", "clone_id", "source_file")]

## ###
## May be a good idea to pull "junction_aa" here, as well(?)
## ###

# Hold another frame for collapsing VHHs
tmp_frame <- tmp[ ,c("VHH", "duplicate_count")]
# Drop duplicate count, as we will add it back later
tmp$duplicate_count <- NULL
# Aggregate the VHH's by their duplicate count, so we can get exact numbers on them
tmp_frame <- aggregate(duplicate_count ~ VHH, data = tmp_frame, FUN = sum)
# Merge back to other data on VHH, so that each VHH's whole duplciate count is known
tmp <- merge(tmp_frame, tmp, by = "VHH")
# Drop duplicates
tmp <- tmp[!duplicated(tmp$VHH),]


# Rename file_name variable
file_name <- gsub("-rule-pooled-families.csv", "-pooled-families_VHH_database.csv", file_name)

# Write
# cat("Writing... ")
## write.csv(tmp, file = paste(output_dir, file_name, sep = ''), row.names = F)
# Elevator Music
cat("Finished writing: ", file_name, '\n', sep = '')

## Instead of writing, assign tmp to dat and we will reload the data to check
## dat <- tmp

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## Important Variables at the end of Pooling:
# dat <- all the important data we want
# file_name <- the input file
# root_path <- folder the input file is in
# output_dir

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## Route the traffic out
if (grepl("-pooled-families_VHH_database.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-MMseqs2.R", file_name), '\n')
  # system(paste("Rscript CleanCluster-MMseqs2.R", file_name))
}

## Cleanup dfs (results are just summary stats)
rm(tmp_frame)

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##














## *** *** *** *** *** *** *** *** *** *** *** *** *** ##
## *** *** *** *** toAA  CHECK *** *** *** *** *** *** ##
## *** *** *** *** *** *** *** *** *** *** *** *** *** ##

# output_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"

# tmp <- read.csv(paste0(output_dir, file_name), header = T, stringsAsFactors = F)

# check_statistics(dat, tmp)

#table(tmp$clone_id, dat$clone_id)


## *** *** *** *** *** *** *** *** *** *** *** *** *** ##



























## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##


## CleanCluster-MMseqs2


# AA-DB-mmseqs2-Clustering

# Input: "20_VecF2_R3_S20-rule-pooled-families_VHH_database.csv"
# Output: "20_VecF2_R3_S20-rule-pooled-families_subclustered_database.csv"

# store directory info
root_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"
output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/tmp/"

# Load the file ("20_VecF2_R3_S20-rule-pooled-families_VHH_database.csv")
## tmp <- read.csv(paste(root_dir, file_name, sep = ''), header = T, stringsAsFactors = F)

# Loading the file means assigning dat to tmp
tmp <- dat

# rm(dat)


## ### 
## To AA moved here
## ###


tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)

# Trim sequences (or not yet?)
tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))

# Drop unproductive sequences
all_seqs <- nrow(tmp)
# Keep productive only
tmp <- unproductive(tmp)
filtered_seqs <- nrow(tmp)

# Elevator music
cat(file_name, '\n')
cat(round(100 - (filtered_seqs/all_seqs*100), 2), "% of sequences unproductive.\n", sep = '')

## ###
## END to-AA
## ###




holder <- tmp













# Create a data frame for holding results
results <- data.frame()
# Iterate over Hierarchically Clustered clone_ids, re-clustering
for (nm in unique(tmp$clone_id)[1:2]){
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
  
  ## ###
  ## Here, we need to fix it to make sure that identical centroid_ids dont exist
  ## Actually, is that important given that centroid IDs are nested in Clonal Families?
  ## ###
  
  # Create a centroid id
  # dat$centroid_id <- as.numeric(factor(dat$centroid))
  # Centroid ID will be family ID plus unique centroid identifier
  dat$centroid_id <- paste(nm, '-', as.numeric(factor(dat$centroid)), sep = '')
  
  ## ### Fix probably ends here
  
  # We're going to use the sequence ids as keys, so we want to clean them
  tmp$sequence_id <- clean_seq_ids(tmp$sequence_id)
  dat$centroid <- clean_seq_ids(dat$centroid)
  dat$sequence_id <- clean_seq_ids(dat$sequence_id)
  # merge results
  dat <- merge_db_to_centroids(dat, tmp)
  # Combine results and move on
  results <- rbind(results, dat)
}

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## Important Variables at the end of Pooling:
# dat <- all the important data we want
# file_name <- the input file
# root_path <- folder the input file is in
# output_dir

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##


file_name <- gsub("VHH_database.csv", "subclustered_database.csv", file_name)

cat("Writing: ", paste(root_dir, file_name, sep = ''), '\n', sep = '')

## write.csv(results, file = paste(root_dir, file_name, sep = ''), row.names = F)

## Route the traffic out
if (grepl("_subclustered_database.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-Refining.R", file_name), '\n')
##   system(paste("Rscript CleanCluster-Refining.R", file_name))
}


## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##
































## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## ###
## BEGIN "[VHH]-Cluster-Refining"
## ###

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
args <- commandArgs(trailingOnly=TRUE)

# args <- c("GS_NGO0265_d2VHH_09Jul24-rule-pooled-families_subclustered_database.csv")

# We only need one file name to pull results and outliers
file_name <- args[1]

#Input: 19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv
#Output:
#(1) -> S19-pooled-clean-database.csv
#(2) -> S19-pooled-outlier-database.csv

results_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-clean-database.csv", file_name)
outliers_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-outlier-database.csv", file_name)

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

cat("Loading: ", paste(root_dir, file_name, sep = ''), '\n', sep = '')

# dat <- results
dat <- read.csv(paste(root_dir, file_name, sep = ''), header = T, stringsAsFactors = F)

cat(file_name, " has ", length(unique(dat$clone_id)), " clonal families.\n", sep = '')

### Family Statistics:
# Add duplicate info here, also?

# Hold results in results

### Sub-Family Statistics:

### Note: Cleaning can probably be added here, too?
#dat <- results

results <- data.frame()

for (clone_id in unique(dat$clone_id)){
  # Subset clone id
  tmp <- dat[dat$clone_id == clone_id,]
  
  unrolled_CDR3s <- unroll_dupes_fast(tmp, "junction_aa")
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  
  tmp$family_median_CDR3_homology <- 100 - (as.vector(adist(median_CDR3, tmp$junction_aa))/nchar(median_CDR3))*100
  
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
results <- results[order(results$subfamily_median_VHH_homology, 
                         results$duplicate_count, 
                         decreasing = T),]

# Order statistic
results$order_index <- 1:nrow(results)

## Iterate Over families and centroids, culling members for further clustering





# holder <- results



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
  
  ## family_median_CDR3_homology 
  
  # For each centroid id
  for (nm in unique(tmp$centroid_id)){
    cat("Subcluster ID:", nm, '\n')
    # Subset centroid id
    tmp_frame <- tmp[tmp$centroid_id == nm,]
    # Check if junction_aa is defined
    if (any(is.na(tmp_frame$family_median_CDR3_homology))){
      cat("No CDR3 Homology for Subcluster: ", nm, '\n')
      cat("Skipping...\n")
      next
    }
    # Calculate quick stats for rule-based cleaning
    CDR3_length <- nchar(tmp_frame$junction_aa[1])
    # 
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
cat("Percent VHH outliers: ", round(nrow(outliers)/nrow(dat)*100, 2), '\n')
cat("Subfamily outliers:\n")
print(table(outliers$clone_id)[order(table(outliers$clone_id), decreasing = T)])

root_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"

results_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-clean-database.csv", file_name)
outliers_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-outlier-database.csv", file_name)

cat("Writing: ", paste(root_dir, results_file, sep = ''), '\n', sep = '')
#write.csv(results, 
#          file = paste(root_dir, results_file, sep = ''), 
#          row.names = F)
cat("Writing: ", paste(root_dir, outliers_file, sep = ''), '\n', sep = '')
#write.csv(outliers, 
#          file = paste(root_dir, outliers_file, sep = ''), 
#          row.names = F)


## --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## Important Variables at the end of Pooling:
# dat <- all the important data we want
# file_name <- the input file
# root_path <- folder the input file is in
# output_dir

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## Route the traffic out
if (grepl("-pooled-families_subclustered_database.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-Recluster.R", file_name), '\n')
  system(paste("Rscript CleanCluster-Recluster.R", file_name))
  #Rscript CleanCluster-Recluster.R 19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv 
}

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##







































## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##


## ###
## BEGIN "[VHH]-Outlier-Re-Clustering"
## ###

#Writing: /Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/19_VecF1_R3_S19-rule-pooled-clean-database-H90.csv
#Writing: /Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/19_VecF1_R3_S19-rule-pooled-outlier-database-H90.csv
#Rscript CleanCluster-Recluster.R 19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv 
#Error in file(file, "rt") : cannot open the connection
#Calls: read.csv -> read.table -> file
#In addition: Warning message:
#  In file(file, "rt") :
#  cannot open file '/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv': No such file or directory
#Execution halted

## Load Utilities
source("~/VHH/lib/utilities.R")
# library(shazam)

## Draw file info from commandArgs
args <- commandArgs(trailingOnly=TRUE)
# args <- c("19_VecF1_R3_S19-pooled-families_subclustered_database.csv")
# args <- c("19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv")
# file_name <- "20_VecF2_R3_S20-rule-pooled-families_subclustered_database.csv"# 

#args <- c()
#args[1] <- "Global_Clone_Experiment_master_set-D1-09-10-11-rule-pooled-families_subclustered_database.csv"

file_name <- args[1]

#Inputs: 
#(1) -> 20_VecF2_R3_S20-pooled-families_subclustered_database.csv
#(2) -> 20_VecF2_R3_S20_min_threshold_clone-pass.csv"
#Output:
#(1) -> S20-pooled-outlier-reclustered-database.csv

results_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-clean-database.csv", file_name)
outliers_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-outlier-database.csv", file_name)

## Load Outliers
root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"

cat("Loading: ", paste(root_dir, outliers_file, sep = ''), '\n', sep = '')

# Load outliers
dat <- read.csv(paste(root_dir, outliers_file, sep = ''), header = T, stringsAsFactors = F)

cat("Loading: ", paste("~/VHH/data/Cluster-DBs/phylogenetic_dbs/", 
                       gsub("-rule-pooled-families_subclustered_database.csv", 
                            "_min_threshold_clone-pass.csv", 
                            file_name), sep = ''), '\n', sep = '')

## Load Source Database
tmp <- read.csv(paste("~/VHH/data/Cluster-DBs/phylogenetic_dbs/", 
                      gsub("-rule-pooled-families_subclustered_database.csv", 
                           "_min_threshold_clone-pass.csv", 
                           file_name), sep = ''),
                header = T, stringsAsFactors = F)

# Ensure that, if no sequence_id exists, we create it
if (!("sequence_id" %in% colnames(tmp))){
  # Create the sequence_id
  tmp$sequence_id <- paste0(1:nrow(tmp), '-',tmp$junction_aa)
  tmp <- tmp[,c("sequence_id", colnames(tmp)[1:(ncol(tmp)-1)])]
}

# Trim CDR3s
tmp$junction_aa <- trim_CDR3(tmp$junction_aa)
# Translate VHHs
tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
# Trim sequences
tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
# We need those VHHs from tmp which are represented in outlier file
dat <- tmp[tmp$VHH %in% dat$VHH,]
## Prepare for clonal inference

dat$clone_id <- NULL
dat$dist_nearest <- NULL

## Clonal Inference

cat("Calculating distToNearest...\n")

db <- shazam::distToNearest(dat, model = "ham", 
                            normalize = "len", 
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
    
    # dat$centroid_id <- as.numeric(factor(dat$centroid))
    dat$centroid_id <- paste(nm, '-', as.numeric(factor(dat$centroid)), sep = '')
    
    # We're going to use the sequence ids as keys, so we want to clean them
    db$sequence_id <- clean_seq_ids(db$sequence_id)
    dat$centroid <- clean_seq_ids(dat$centroid)
    dat$sequence_id <- clean_seq_ids(dat$sequence_id)
    # merge results
    dat <- merge_db_to_centroids(dat, db)
  } else {
    # otherwise, we'll make up a centroid id and a sequence id, and
    
    #tmp_frame$centroid_id <- tmp_frame$sequence_id[1]
    tmp_frame$centroid_id <- paste(nm, '-', substr(tmp_frame$sequence_id[1], nchar(tmp_frame$sequence_id[1])-10, nchar(tmp_frame$sequence_id[1])), sep = '')
    
    tmp_frame$centroid_sequence_id <- tmp_frame$sequence_id[1]
    tmp_frame$centroid <- tmp_frame$VHH[1]
    # Now we assign to dat for rbind
    dat <- tmp_frame
  }
  # Combine results and move on
  results <- rbind(results, dat)
}


#results$centroid_id <- as.numeric(factor(results$centroid_id))

### Temporary Write

# Note:
# input_file <- "19_VecF1_R3_S19-pooled-outlier-database-H90.csv"
file_name <- gsub("_database.csv", "-reclustered-database.csv", file_name)

cat("Writing: ", paste(root_dir, file_name, sep = ''), '\n', sep = '')

#write.csv(results, 
#          file = paste(root_dir, file_name, sep = ''), row.names = F)

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## Important Variables at the end of Pooling:
# dat <- all the important data we want
# file_name <- the input file
# root_path <- folder the input file is in
# output_dir

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## Route the traffic out
if (grepl("-reclustered-database.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-Union.R", file_name), '\n')
  system(paste("Rscript CleanCluster-Union.R", file_name))
}

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##



































## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## ###
## BEGIN [VHH]-Outlier-DB-Statistics
## ###

# Rscript CleanCluster-Union.R 19_VecF1_R3_S19-pooled-families_subclustered-reclustered-database.csv

## Clean data file:
# Writing: /Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/19_VecF1_R3_S19-rule-pooled-clean-database-H90.csv
## ORIGINAL outlier file
# Writing: /Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/19_VecF1_R3_S19-rule-pooled-outlier-database-H90.csv
## Re-Clustered outlier file:
# Writing: /Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/19_VecF1_R3_S19-pooled-families_subclustered-reclustered-database.csv

## Load Utilities
source("~/VHH/lib/utilities.R")

# Draw file info from commandArgs
args <- commandArgs(trailingOnly=TRUE)


#args <- c("19_VecF1_R3_S19-rule-pooled-families_subclustered-reclustered-database.csv")

#args <- c("20_VecF2_R3_S20-rule-pooled-families_subclustered-reclustered-database.csv")

file_name <- args[1]

## Load Data
root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"
clean_file <- gsub("-pooled-families_subclustered-reclustered-database.csv", "-pooled-clean-database.csv", file_name)
outlier_file <- file_name

cat("Loading: ", paste(root_dir, clean_file, sep = ''), '\n', sep = '')

# Load cleaned database
clean_db <- read.csv(paste(root_dir, clean_file, sep = ''), header = T, stringsAsFactors = F)

cat("Loading: ", paste(root_dir, outlier_file, sep = ''), '\n', sep = '')

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
outlier_test_db$clone_id <- paste("O-", outlier_test_db$clone_id, sep = '')#outlier_test_db$clone_id + 100000
outlier_test_db$centroid_id <- paste("O-", outlier_test_db$centroid_id, sep = '')#outlier_test_db$centroid_id + 100000

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
output_file <- gsub("-pooled-clean-database.csv", "-Pooled-Cleaning-Reclustering-Initial-Results.csv", clean_file)

cat("Writing: ", paste(root_dir, output_file, sep = ''), "\n", sep = '')

# write.csv(results, file = paste(root_dir, output_file, sep = ''), row.names = F)




## Route the traffic out
if (grepl("-reclustered-database.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-Pool.R", output_file), '\n')
  system(paste("Rscript CleanCluster-Pool.R", output_file))
} 



