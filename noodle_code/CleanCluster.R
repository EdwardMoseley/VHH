# CleanCluster Algorithm
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
if (grepl("-Pooled-Cleaning-Reclustering-Initial-Results.csv", args[1])) {
  file_name <- args[1]
  cat("Final Processing of: ", file_name, '\n', sep = '')
  # Define root_path
  root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"
  # Load Data
  dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
} else if (grepl("_min_threshold_clone-pass.csv", args[1])){
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
# file_name <- the input file
# root_path <- folder the input file is in

## --- --- --- --- --- --- --- --- --- --- --- --- --- ##


## Route the traffic out
if (grepl("_min_threshold_clone-pass.csv", file_name)){
  # change file name
  file_name <- gsub("_min_threshold_clone-pass.csv", "-rule-pooled-families.csv", file_name)
  cat("Writing: ")
  # write.csv(dat, file = paste(root_path, file_name, sep = ''), row.names = F)
  # -rule-pooled-families.csv
  cat(paste("Rscript CleanCluster-toAA.R", file_name), '\n')
  # Pass the file to the next script in the pipeline
}

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ##





















































## Begin toAA

## ###
## CleanCluster-toAA
## ###

## Load Utilities
# source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
# args <- commandArgs(trailingOnly=TRUE)

#args <- c("GS_d2_18Jun24-rule-pooled-families.csv")
# file_name <- args[1]

## Amino Acid-only database

# Input data: "20_VecF2_R3_S20-pooled-families.csv"
# Output Data: "translated_dbs/20_VecF2_R3_S20-pooled-families_VHH_database.csv" 

# store directory info
# root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/"
## No need to output data, we have a lot of ram
# output_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"





# LET'S WAIT TO Trim sequences
## ==tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))


# Elevator music
cat("Loading: ", paste(root_path, file_name, sep = ''), '\n', sep = '')
# This is just loading the original file
tmp <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
# Capture source file
tmp$source_file <- rep(file_name, times = nrow(tmp))
# Translate - VHH
tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
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
file_name <- gsub("_min_threshold_clone-pass.csv", "-pooled-families_VHH_database.csv", file_name)

# Write
# cat("Writing... ")
## write.csv(tmp, file = paste(output_dir, file_name, sep = ''), row.names = F)
# Elevator Music
cat("Finished writing: ", file_name, '\n', sep = '')

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
} else if (grepl("later file name for final pooling...", file_name)){
  cat("ERROR in [CleanCluster-toAA.R]: Wrong file name?")
}