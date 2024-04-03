## ###
## Build Clonal Databases
## Contains information and statistics at the clonal family level
## 15Mar24
## ###

# Load utilities
source("~/VHH/lib/utilities.R")

data_dir <- "~/VHH/data/Cluster-DBs"
res_dir <- "~/VHH/data/Cluster-DBs/cluster_database"

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
  tab <- table(dat$clone_id) # top_groups(dat, "clone_id")
  # Instantiate a counter
  counter <- 0
  # Reinstantiate tmp_res
  tmp_res <- data.frame()
  
  for (nm in names(tab)){
    # Remind where we are
    if (counter %% 100 == 0){
      cat("Progress: ", round(counter/length(tab), 2), "%", '\n', sep = '')
    }
    # Increment Counter
    counter <- counter + 1
    # Subset the clones of interest
    tmp_frame <- dat[dat$clone_id == as.numeric(nm),]
    tmp_frame$junction_aa <- trim_CDR3(tmp_frame$junction_aa)
    # Keep the integer as the cluster_id
    cluster_id <- nm
    
    ## Note: here we need to ensure that we account for the 
    unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
    
    # Make consensus sequence
    consensus_sequence <- strict_consensus(unrolled_CDR3s)
    
    # Calculate entropy
    CDR3_entropy <- string_entropy(unrolled_CDR3s)
    
    # Sum duplicate count for the number of TOTAL DNA seqs
    Dup_count <- sum(tmp_frame$duplicate_count)
    # Unique DNA sequences
    DNA_seqs_in_cluster <- length(unique(tmp_frame$sequence))
    # Unique CDR3s
    CDR3s_in_cluster <- length(unique(tmp_frame$junction_aa))
    
    # Bind results
    tmp_res <- rbind(tmp_res, cbind(f,
                                    cluster_id,
                                    consensus_sequence,
                                    Dup_count,
                                    DNA_seqs_in_cluster,
                                    CDR3s_in_cluster,
                                    CDR3_entropy))
    
  }
  
  cat("Completed file:", f, '\n')
  file_out <- paste(res_dir, gsub('.tsv', '-consensuses15Mar24.csv', f), sep = '/')
  cat("Writing:", file_out, '\n')
  #write.csv(tmp_res, file_out, row.names = F)
}