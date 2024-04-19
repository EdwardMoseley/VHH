## ###
## CleanCluster-Pool
## First, and final, step in generating 
## clonal clusters for expert analysis.
## ###

## First pass call:
# Rscript CleanCluster-Pool.R 19_VecF1_R3_S19_min_threshold_clone-pass.csv
## Last pass call:
# Rscript CleanCluster-Pool.R 19_VecF1_R3_S19-rule-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv 

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
#args <- commandArgs(trailingOnly=TRUE)

#args <- "19_VecF1_R3_S19_min_threshold_clone-pass.csv"
args <- "19_VecF1_R3_S19-rule-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv"

# test if there is at least one argument: if not, return an error
if (grepl("-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv", args[1])) {
  file_name <- args[1]
  cat("Processing: ", file_name, '\n', sep = '')
  # Define root_path
  root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"
  # Load Data
  dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
} else if (grepl("_min_threshold_clone-pass.csv", args[1])){
  file_name <- args[1]
  cat("Processing: ", file_name, '\n', sep = '')
  # Define root_path
  root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/"
  # Load data
  dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)
  # Trim CDR3s
  dat$junction_aa <- trim_CDR3(dat$junction_aa)
} else{
  stop("Please provide the directory to the data file (in [prefix]_min_threshold_clone-pass.csv format) OR...", call.=FALSE)
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
  for (nm in as.numeric(names(tablet)[2:length(tablet)])){
    # Subset the old clone_ids and replace with the largest group to which they could belong ([1])
    dat[dat$clone_id == nm,]$clone_id <- rep(as.numeric(names(tablet[1])), each = nrow(dat[dat$clone_id == nm,]))
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
tmp <- tmp[order(as.numeric(tmp$abundance), decreasing = T),]


## Pooling with % Homology

#done_list <- vector()
# Begin merging routine
# For each unique CDR3
#for (nm in unique(tmp$median_CDR3)){
#  cat("Starting: ", nm, '\n', sep = '')
#  # Pull ALL observations involving that CDR3
#  tmp_frame <- tmp[(tmp$median_CDR3 == nm | tmp$variable == nm),]
#  # Check if we've done this one yet
#  if (!(nm %in% done_list)){
#    cat("Merging: ", nm, '\n', sep = '')
#    # If not, add it to the list
#    done_list <- c(done_list, tmp_frame$median_CDR3)
#    # Begin adding lower-abundance groups to greater-abundance ones
#    for (i in 2:(nrow(tmp_frame))){
#      dat[dat$clone_id == tmp_frame$clone_id[i],]$clone_id <- tmp_frame$clone_id[1]
#    }
#  } else {
#    # If we have done it, we pass over it
#    cat("Passing: ", nm, '\n', sep = '')
#  }
#}

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
tmp <- tmp[order(as.numeric(tmp$abundance), decreasing = T),]

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

## Route the traffic out
if (grepl("_min_threshold_clone-pass.csv", file_name)){
  file_name <- gsub("_min_threshold_clone-pass.csv", "-rule-pooled-families.csv", file_name)
  write.csv(dat, file = paste(root_path, file_name, sep = ''), row.names = F)
  cat(paste("Rscript CleanCluster-toAA.R", file_name), '\n')
  # Pass the file to the next script in the pipeline
  system(paste("Rscript CleanCluster-toAA.R", file_name))
} else if (grepl("-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv", file_name)){
  # change file name
  file_name <- gsub("-Pooled-Cleaning-Reclustering-Initial-Results-H90.csv", "-H90-Cleaned19Apr24.csv", file_name)
  # write
  write.csv(dat, file = paste(root_path, file_name, sep = ''), row.names = F)
}


