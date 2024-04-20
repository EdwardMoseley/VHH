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
# args <- commandArgs(trailingOnly=TRUE)


#args <- c("19_VecF1_R3_S19-rule-pooled-families_subclustered-reclustered-database.csv")

args <- c("20_VecF2_R3_S20-rule-pooled-families_subclustered-reclustered-database.csv")

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

write.csv(results, file = paste(root_dir, output_file, sep = ''), row.names = F)

## ###
## END [VHH]-Outlier-DB-Statistics
## ###


## ###
## Here, we want to return the data to the pooling script for a final pooling run
## ###

## Route the traffic out
if (grepl("-reclustered-database.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-Pool.R", output_file), '\n')
  system(paste("Rscript CleanCluster-Pool.R", output_file))
} else if (grepl("later file name for final pooling...", file_name)){
  cat("ERROR in [CleanCluster-toAA.R]: Wrong file name?")
  # change file name
  # write
}