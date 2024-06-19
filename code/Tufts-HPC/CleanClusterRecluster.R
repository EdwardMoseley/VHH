## ###
## BEGIN "[VHH]-Outlier-Re-Clustering"
## ###

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
args <- commandArgs(trailingOnly=TRUE)

file_name <- args[1]

#Inputs: 
#(1) -> GS_sample_pooled_clean_database.csv
#(2) -> GS_sample_pooled_outlier_database.csv
#Output:
#(1) -> S20-pooled-outlier-reclustered-database.csv

results_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-clean-database.csv", file_name)
outliers_file <- gsub("-pooled-families_subclustered_database.csv", "-pooled-outlier-database.csv", file_name)

## Load Outliers
root_dir <- "/cluster/tufts/cowenlab/emosel01/data/forClustering/tmp/"

cat("Loading: ", paste(root_dir, outliers_file, sep = ''), '\n', sep = '')

# Load outliers
dat <- read.csv(paste(root_dir, outliers_file, sep = ''), header = T, stringsAsFactors = F)

cat("Loading: ", paste("/cluster/tufts/cowenlab/emosel01/data/forClustering/tmp/", 
                       gsub("-rule-pooled-families_subclustered_database.csv", 
                            "_min_threshold_clone-pass.csv", 
                            file_name), sep = ''), '\n', sep = '')

## Load Source Database
tmp <- read.csv(paste("/cluster/tufts/cowenlab/emosel01/data/forClustering/tmp/", 
                      gsub("-rule-pooled-families_subclustered_database.csv", 
                           "_min_threshold_clone-pass.csv", 
                           file_name), sep = ''),
                header = T, stringsAsFactors = F)

# Ensure that, if no sequence_id exists, we create it
if (!("sequence_id" %in% colnames(tmp))){
  # Create the sequence_id
  tmp$sequence_id <- paste0(1:nrow(tmp), '-',tmp$junction_aa)
  tmp <- tmp[,c("sequence_id", colnames(tmp)[1:(ncol(tmp)-1)])]
}

# Trim CDR3s
tmp$junction_aa <- trim_CDR3(tmp$junction_aa)
# Translate VHHs
tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
# Trim sequences
tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
# We need those VHHs from tmp which are represented in outlier file
dat <- tmp[tmp$VHH %in% dat$VHH,]
## Prepare for clonal inference

dat$clone_id <- NULL
dat$dist_nearest <- NULL

## Clonal Inference

cat("Calculating distToNearest...\n")

db <- shazam::distToNearest(dat, model = "ham", 
                    normalize = "len", 
                    vCallColumn = "v_call", 
                    jCallColumn = "j_call",
                    nproc = 1)

# Determine threshold (based on median)
thr <- summary(db$dist_nearest)[["Median"]]

new_thr <- thr + (thr*0.25)
print("25% Higher Threshold:")
print(round(new_thr, 2))

cat("Inferring Clones:\n")
system.time(
  # Clonal assignment using hierarchical clustering
  db <- scoper::hierarchicalClones(db, 
                                   threshold=new_thr, 
                                   v_call="v_call", 
                                   j_call = "j_call")
)

db <- db@db

## Sub-cluster with MMseqs2
output_dir <- "/Users/emoseley/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/tmp/"

results <- data.frame()

# Iterate over Hierarchically Clustered clone_ids, re-clustering
for (nm in unique(db$clone_id)){
  # nm <- 8 (51 obs)
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


#results$centroid_id <- as.numeric(factor(results$centroid_id))

### Temporary Write

# Note:
# input_file <- "19_VecF1_R3_S19-pooled-outlier-database-H90.csv"
file_name <- gsub("_database.csv", "-reclustered-database.csv", file_name)

cat("Writing: ", paste(root_dir, file_name, sep = ''), '\n', sep = '')

write.csv(results, 
          file = paste(root_dir, file_name, sep = ''), row.names = F)


## Route the traffic out
if (grepl("-reclustered-database.csv", file_name)){
  # Pass the file to the next script in the pipeline
  cat(paste("Rscript CleanCluster-Union.R", file_name), '\n')
  system(paste("Rscript CleanCluster-Union.R", file_name))
} else if (grepl("later file name for final pooling...", file_name)){
  cat("ERROR in [CleanCluster-toAA.R]: Wrong file name?")
  # change file name
  # write
}