## CleanClusterMMseqs2.R

# AA-DB-mmseqs2-Clustering

## Load Utilities
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

## Draw file info from commandArgs
args <- commandArgs(trailingOnly=TRUE)

file_name <- args[1]

# Input: "20_VecF2_R3_S20-rule-pooled-families_VHH_database.csv"
# Output: "20_VecF2_R3_S20-rule-pooled-families_subclustered_database.csv"

# store directory info
root_dir <- "/cluster/tufts/cowenlab/emosel01/data/forClustering/tmp/"
output_dir <- "/cluster/tufts/cowenlab/emosel01/data/forClustering/tmp/"

# Load the file
tmp <- read.csv(paste(root_dir, file_name, sep = ''), header = T, stringsAsFactors = F)
# Create a data frame for holding results
results <- data.frame()
# Iterate over Hierarchically Clustered clone_ids, re-clustering
for (nm in unique(tmp$clone_id)[1:10]){
  # Subset observations
  tmp_frame <- tmp[tmp$clone_id == nm,]
  
  ## ###
  ## In this for-loop, it would be a good idea to check if there is more than one clone in the family first!
  ## ###
  
  # write fasta for mmseqs2
  write_fasta(tmp_frame, "VHH", "sequence_id", paste(output_dir, "tmp.fasta", sep = ''))
  # Prepare mmseqs command
  command <- paste("mmseqs ", 
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

file_name <- gsub("VHH_database.csv", "subclustered_database.csv", file_name)

cat("Writing: ", paste(root_dir, file_name, sep = ''), '\n', sep = '')

#Save Results
write.csv(results, file = paste(root_dir, file_name, sep = ''), row.names = F)

