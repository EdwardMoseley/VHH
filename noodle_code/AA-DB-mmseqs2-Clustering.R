# AA-DB-mmseqs2-Clustering

# Imports
source("~/VHH/lib/utilities.R")
suppressPackageStartupMessages(library(stringdist))

# store directory info
root_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/"
output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/tmp/"

# Load up data
files <- list.files(root_dir)
files <- files[grepl("VHH_database.csv", files)]

for (f in files){
  # Load the file
  tmp <- read.csv(paste(root_dir, f, sep = ''), header = T, stringsAsFactors = F)
  # Create a data frame for holding results
  results <- data.frame()
  # Iterate over Hierarchically Clustered clone_ids, re-clustering
  for (nm in unique(tmp$clone_id)){
    # nm <- 29
    # Subset observations
    tmp_frame <- tmp[tmp$clone_id == nm,]
    # write fasta for mmseqs2
    write_fasta(tmp_frame, "VHH", "sequence_id", paste(output_dir, "tmp.fasta", sep = ''))
    # Prepare mmseqs command
    command <- paste("/usr/local/Cellar/mmseqs2/15-6f452/bin/mmseqs ", 
                     "easy-cluster -v 3 --min-seq-id 0.95 ", 
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
  write.csv(results, file = paste(root_dir, gsub("VHH_database.csv", "subclustered_database.csv", f), sep = ''), row.names = F)
}

