# CleanCluster Algorithm

## Change Library path to ensure we have access

#.libPaths("/cluster/tufts/hpc/tools/R/4.3.0")
#library("alakazam")
#library("scoper")
#library("shazam")
#source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

# Set root dir

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
args <- c("GS_XLI_d2VHH_21Aug24_min_threshold_clone-pass.csv")
#args <- c("GS_JPL_d2VHH_21Aug24_min_threshold_clone-pass.csv")

file_name <- args[1]
cat("Initial Processing of: ", file_name, '\n', sep = '')
# Define root_path
root_path <- "~/Downloads/"
# Load data
dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
# Trim CDR3s
dat$junction_aa <- trim_CDR3(dat$junction_aa)

# Translate alignment
dat$VHH <- fixed_translate_DNA(dat$sequence_alignment)

# Trim sequences
dat$VHH <- substr(dat$VHH, 11, nchar(dat$VHH))

# Drop unproductive sequences
all_seqs <- nrow(dat)
# Keep productive only
dat <- unproductive(dat)
filtered_seqs <- nrow(dat)

# Elevator music
cat(file_name, '\n')
cat(round(100 - (filtered_seqs/all_seqs*100), 5), "% of sequences unproductive (due to a premature stop codon).\n", sep = '')

# Print clone counts
cat(file_name, " has ", length(unique(dat$clone_id)), " clonal families after removing unproductive sequences.\n", sep = '')

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




## Cleanup dfs (results are just summary stats)
rm(df, tmp, tmp_frame, results)

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##



# > file_name
# [1] "20_VecF2_R3_S20-rule-pooled-families.csv"

file_name <- gsub("_min_threshold_clone-pass.csv", "-processed-pooled.csv", file_name)
cat("Writing: ", paste0(root_path, file_name), '\n')
write.csv(dat, file = paste0(root_path, file_name), row.names = F)



































## *** *** *** *** *** *** *** *** *** *** *** *** *** ##
## *** *** *** ***  POOLING CHECK  *** *** *** *** *** ##
## *** *** *** *** *** *** *** *** *** *** *** *** *** ##

# tmp <- read.csv(paste0(root_path, file_name), header = T, stringsAsFactors = F)

# check_statistics(dat, tmp)

# print_statistics(dat)

## *** *** *** *** *** *** *** *** *** *** *** *** *** ##



## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##


## Begin toAA

## ###
## CleanCluster-toAA
## ###

# Here, the focus is to translate the data to amino acids
# and remove unproductive sequence information before
# calling mmseqs2 for subclustering by homology



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
# tmp <- dat

## Here, we want to cluster based on the Amino Acids
tmp <- dat

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

file_name <- gsub("-processed-pooled.csv", "-subclustered.csv", file_name)
cat("Writing: ", paste(root_path, file_name, sep = ''), '\n', sep = '')
write.csv(results, file = paste(root_path, file_name, sep = ''), row.names = F)


# results <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)



## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## Reload from MMseqs2
#results <- read.csv("~/Downloads/Example Materials/20_VecF2_R3_S20-pooled-families_subclustered_database05Aug24.csv",
#                    header = T, stringsAsFactors = F)

#length(unique(test$centroid_id))
#length(unique(results$centroid_id))

#> length(unique(test$centroid_id))
#[1] 1733
#> length(unique(results$centroid_id))
#[1] 1686

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##








## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##

## ###
## BEGIN "[VHH]-Cluster-Refining"
## ###


# args <- c("GS_NGO0265_d2VHH_09Jul24-rule-pooled-families_subclustered_database.csv")

#Input: 19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv
#Output:
#(1) -> S19-pooled-clean-database.csv
#(2) -> S19-pooled-outlier-database.csv

#results_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-clean-database.csv", file_name)
#outliers_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-outlier-database.csv", file_name)

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

# root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"

# cat("Loading: ", paste(root_dir, file_name, sep = ''), '\n', sep = '')

# dat <- results
# dat <- read.csv(paste(root_dir, file_name, sep = ''), header = T, stringsAsFactors = F)

# cat(file_name, " has ", length(unique(dat$clone_id)), " clonal families.\n", sep = '')

### Family Statistics:
# Add duplicate info here, also?

# Hold results in results

### Sub-Family Statistics:

### Note: Cleaning can probably be added here, too?
#dat <- results



## ###
## GENERATE STATISTICS
## ###

## Store data for now
dat <- results

cat(file_name, " has ", length(unique(dat$clone_id)), " clonal families.\n", sep = '')

results <- data.frame()

for (clone_id in unique(dat$clone_id)){
  cat("Clone ID:", clone_id, '\n')
  # Subset clone id
  tmp <- dat[dat$clone_id == clone_id,]
  # Unroll duplicate CDR3s
  unrolled_CDR3s <- unroll_dupes_fast(tmp, "junction_aa")
  # Pull the most frequent CDR3 for the family
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  # Calculate CDR3 Homology
  tmp$family_median_CDR3_homology <- 100 - (as.vector(adist(median_CDR3, tmp$junction_aa))/nchar(median_CDR3))*100
  
  # For each centroid id
  for (nm in unique(tmp$centroid_id)){
    cat("Centroid ID:", nm, '\n')
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





## ###
## REMOVE OUTLIERS
## ###

## Iterate Over families and centroids, culling members for reclustering

holder <- results
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
    # Check if junction_aa is defined
    if (any(is.na(tmp_frame$family_median_CDR3_homology))){
      cat("No CDR3 Homology for Subcluster: ", nm, '\n')
      cat("Skipping...\n")
      next
    }
    # Calculate quick stats for rule-based cleaning
    CDR3_length <- nchar(tmp_frame$junction_aa[1])
    # Calculate mean homology
    avg_CDR3_percent_homology <- mean(tmp_frame$family_median_CDR3_homology)
    
    
    # Rule-based Cleanup here
    
    ## OUTLIERS
    
    if (CDR3_length < 5 & avg_CDR3_percent_homology < 100){
      outliers <- rbind(outliers, tmp_frame)
    } else if (CDR3_length < 10 & avg_CDR3_percent_homology < 80){
      outliers <- rbind(outliers, tmp_frame)
    } else if (CDR3_length >= 10 & avg_CDR3_percent_homology < 60){
      outliers <- rbind(outliers, tmp_frame)
    
    ## ELSE PUT THEM IN RESULTS
    } else {
      results <- rbind(results, tmp_frame)
    }
  }
  
}



## Check Results
cat("Percent VHH outliers: ", round(nrow(outliers)/nrow(dat)*100, 2), '\n', sep = '')

cleaned <- results

cat("Outliers plus clean sequences equal input data: ", nrow(cleaned) + nrow(outliers) == nrow(dat), '\n', sep = '')

#write.csv(cleaned, file = paste0(root_path, gsub("-subclustered.csv", "-cleaned-set.csv", file_name)), row.names = F)
#write.csv(outliers, file = paste0(root_path, gsub("-subclustered.csv", "-outlier-set.csv", file_name)), row.names = F)




## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##


## ###
## BEGIN "[VHH]-Outlier-Re-Clustering"
## ###

# args <- c("19_VecF1_R3_S19-pooled-families_subclustered_database.csv")
# args <- c("19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv")
# file_name <- "20_VecF2_R3_S20-rule-pooled-families_subclustered_database.csv"# 

#Inputs: 
#(1) -> 20_VecF2_R3_S20-pooled-families_subclustered_database.csv
#(2) -> 20_VecF2_R3_S20_min_threshold_clone-pass.csv"
#Output:
#(1) -> S20-pooled-outlier-reclustered-database.csv

#results_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-clean-database.csv", file_name)
#outliers_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-outlier-database.csv", file_name)

# Load outliers (they're already here)
dat <- outliers # read.csv(paste(root_dir, outliers_file, sep = ''), header = T, stringsAsFactors = F)

# Ensure that, if no sequence_id exists, we create it
if (!("sequence_id" %in% colnames(dat))){
  cat("Adding sequence IDs:\n")
  # Create the sequence_id
  dat$sequence_id <- paste0(1:nrow(dat), '-',dat$junction_aa)
  # Reorder column names to place sequence_id first
  dat <- dat[,c("sequence_id", colnames(dat)[1:(ncol(dat)-1)])]
}


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

## Make the thres
new_thr <- thr + (thr*0.25)
cat("25% Higher Threshold:\n")
cat(round(new_thr, 2), '\n')

cat("Inferring Clones:\n")
system.time(
  # Clonal assignment using hierarchical clustering
  db <- scoper::hierarchicalClones(db, 
                                   threshold=new_thr, 
                                   v_call="v_call", 
                                   j_call = "j_call",
                                   verbose = T)
)

db <- as.data.frame(db@db)

# We no longer need centroid info from the outliers data either
db$centroid_id <- NULL
db$centroid_sequence_id <- NULL

## Sub-cluster with MMseqs2
output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/tmp/"

hclusted_outliers <- db

## Since DB is outliers, we need to clean it
results <- data.frame()

# Iterate over Hierarchically Clustered clone_ids, re-clustering
for (nm in unique(db$clone_id)){
  cat('Clone_ID: ', nm, '\n', sep = '')
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
    dat$centroid_id <- paste(nm, '-', as.numeric(factor(dat$centroid)), sep = '')
    
    # We're going to use the sequence ids as keys, so we want to clean them
    db$sequence_id <- clean_seq_ids(db$sequence_id)
    dat$centroid <- clean_seq_ids(dat$centroid)
    dat$sequence_id <- clean_seq_ids(dat$sequence_id)

    # merge results
    #dat <- merge_db_to_centroids(dat, db)
    
    colnames(dat)[1] <- "centroid_sequence_id"
    
    dat <- merge(dat, db, by = c("sequence_id"))#, "centroid_sequence_id"))
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


file_name <- gsub("-subclustered.csv", "-outliers-subclustered.csv", file_name)
cat("Writing: ", paste(root_path,file_name, sep = ''), '\n', sep = '')
write.csv(results, file = paste(root_path, file_name, sep = ''), row.names = F)




#outliers <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)





## -- --- ##

## Union of Outliers and Cleaned

## Prepend "O-" to indicate outliers
outliers$clone_id <- paste("O-", outliers$clone_id, sep = '')
outliers$centroid_id <- paste("O-", outliers$centroid_id, sep = '')
outliers$centroid_sequence_id <- paste("O-", outliers$centroid_sequence_id, sep = '')

cat(file_name, " has ", length(unique(outliers$clone_id)), " clonal families.\n", sep = '')
cat(file_name, " has ", length(unique(outliers$centroid_sequence_id)), " clonal subfamilies.\n", sep = '')


### Stack and Write
#cleaned$group <- rep("clean", each = nrow(cleaned))
#outliers$group <- rep("outliers", each = nrow(outliers))

dat <- rbind(cleaned, outliers)

file_name <- gsub("-outliers-subclustered.csv", "-pre-final-pool.csv", file_name)
cat("Writing: ", paste(root_path,file_name, sep = ''), '\n', sep = '')
write.csv(dat, file = paste(root_path, file_name, sep = ''), row.names = F)


# dat <- read.csv(paste0(root_path, "GS_Schisto_d2VHH_14Aug24-pre-final-pool.csv"), header = T, stringsAsFactors = F)

#dat <- read.csv("~/Downloa

## --- --- ##

## Final Pool


# Print clone counts
cat(file_name, " has ", length(unique(dat$clone_id)), " clonal families after processing outliers.\n", sep = '')

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

file_name <- gsub("-pre-final-pool.csv", "-cleaned-GS.csv", file_name)
cat("Writing: ", paste(root_path,file_name, sep = ''), '\n', sep = '')
write.csv(dat, file = paste(root_path, file_name, sep = ''), row.names = F)


#tab <- top_groups(dat, "clone_id")

#for (nm in head(names(tab), 5)){
#  tmp_frame <- dat[dat$clone_id == nm,]
#  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
#  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
#  # Pull statistics
#  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
#  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
#  
#  if (nrow(tmp_frame) > 1){
#    cat(paste('Family: ',unique(tmp_frame$clone_id), '\n', sep = ''))
#    # cat(paste("Group: ", unique(tmp_frame$group)), '\n', sep = '')
#    cat("Overall Members (Including Duplicate DNA Reads): ", sum(tmp_frame$duplicate_count), '\n', sep = '')
#    #cat("Number of Sub-families: ", length(unique(tmp_frame$centroid_id)), '\n', sep = '')
#    cat("Unique VHHs: ", length(unique(tmp_frame$VHH)), '\n', sep = '')
#    cat("Unique CDR3s: ", length(unique(tmp_frame$junction_aa)), '\n', sep = '')
#    cat("Most-frequent VHH: ", median_VHH, '\n', sep = '')
#    print(ggseqlogo(as.character(tmp_frame$junction_aa), seq_type='aa', method = 'p') +
#            theme_classic() +
#            ggtitle(paste(" CDR3 Logo Plot for cluster with the most common CDR3 sequence:\n               ",
#                          median_CDR3, sep = '')))
#  }
#}