## ###
## Build Clonal Databases
## Contains information and statistics at the clonal family level
## 15Mar24
## ###

# Load utilities
source("~/VHH/lib/utilities.R")

data_dir <- "~/VHH/data/Cluster-DBs"
res_dir <- "~/VHH/data/Cluster-DBs/results"

files <- list.files(data_dir)[grepl(".tsv", list.files(data_dir))]

## Data frame for storing all results
consensus_clusters <- data.frame()

# Temporary data frame
tmp_res <- data.frame()

# For each database file
for (f in files){
  
  dat <- read.csv(paste(data_dir, f, sep = '/'), header = T, stringsAsFactors = F, sep = '\t')
  
  cat("Loading file: ", f, '\n')
  # for each family in the database
  tab <- top_groups(dat, "clone_id")
  
  for (nm in names(tab)){
    # Subset the clones of interest
    tmp_frame <- dat[dat$clone_id == as.numeric(nm),]
    tmp_frame$junction_aa <- trim_CDR3(tmp_frame$junction_aa)
    cluster_id <- nm
    
    ## Note: here we need to ensure that we account for the 
    unrolled_CDR3s <- unroll_dupes(tmp_frame)
    
    consensus_sequence <- strict_consensus(unrolled_CDR3s)
    
    Dup_count <- sum(tmp_frame$duplicate_count)
    DNA_seqs_in_cluster <- length(unique(tmp_frame$sequence))
    CDR3s_in_cluster <- length(unique(tmp_frame$junction_aa))
    
    # Bind results
    tmp_res <- rbind(tmp_res, cbind(f,
                                    cluster_id,
                                    consensus_sequence,
                                    Dup_count,
                                    DNA_seqs_in_cluster,
                                    CDR3s_in_cluster))
    
  }
  cat("Completed file: ", f, '\n')
  file_out <- paste(res_dir, gsub('.tsv', '-consensuses26Feb24.csv', f), sep = '/')
  write.csv(tmp_res, file_out, row.names = F)
}