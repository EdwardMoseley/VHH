## ###
## BEGIN "[VHH]-Cluster-Refining"
## ###

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
args <- commandArgs(trailingOnly=TRUE)

# args <- c("19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv")

# We only need one file name to pull results and outliers
file_name <- args[1]

#Input: 19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv
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


## Route the traffic out
if (grepl("-pooled-families_subclustered_database-H90.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-Recluster.R", file_name), '\n')
  system(paste("Rscript CleanCluster-Recluster.R", file_name))
  #Rscript CleanCluster-Recluster.R 19_VecF1_R3_S19-rule-pooled-families_subclustered_database-H90.csv 
} else if (grepl("later file name for final pooling...", file_name)){
  cat("ERROR in [CleanCluster-toAA.R]: Wrong file name?")
  # change file name
  # write
}