## ### 
## Pooling CDR3's by MMseqs2 Clustering
## ###

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
# args <- commandArgs(trailingOnly=TRUE)

args <- c("GS_d2_18Jun24_min_threshold_clone-pass.csv")

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













































results <- data.frame()

# Iterate over Hierarchically Clustered clone_ids, re-clustering
for (nm in unique(db$clone_id)){
  # Subset observations
  tmp_frame <- db[db$clone_id == nm,]
  # If the family has more than one member, we will re-cluster
  if (nrow(tmp_frame) > 1){
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