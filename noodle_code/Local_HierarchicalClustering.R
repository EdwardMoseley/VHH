## Local_HierarchicalClustering.R

# CleanCluster Algorithm

## Change Library path to ensure we have access
#HPC .libPaths("/cluster/tufts/hpc/tools/R/4.3.0")

## Load Utilities
source("~/VHH/lib/utilities.R")
#HPC source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

library("alakazam")
library("scoper")
library("shazam")
library("hash")

# Set root dir
# root_dir <- "/cluster/tufts/cowenlab/emosel01/data/rp_vs_ch2/"
root_dir <- "~/Downloads/2024_TUCF_7981/"

# output_file <- "/cluster/tufts/cowenlab/emosel01/data/GS_Schisto_d2VHH_13Aug24_min_threshold_clone-pass.csv"
output_file <- "~/Downloads/2024_TUCF_7981/GS_Schisto_d2VHH_14Aug24_min_threshold_clone-pass.csv"

files <- list.files(root_dir)[grepl(".csv", list.files(root_dir))]

GS_name <- list.files(root_dir)[grepl("GS", list.files(root_dir))]

files <- files[files != GS_name]

cat("Processing files:\n")

for (f in files){
  cat(f, '\n')
}

dat <- data.frame()
for (f in files){
  cat("Loading:", f,'\n')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)#, sep = '\t')
  tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
  tmp <- tmp[,c("sequence", "locus", "sequence_alignment", "junction", "junction_aa", "v_call", "d_call", "j_call", "VHH", "duplicate_count")]
  dat <- rbind(dat, tmp)
}

# count rows
cat("[Dat] Row count after loading:", nrow(dat), "\n")

# Remove unproductive VHHs
dat <- unproductive(dat)

# count rows
cat("[Dat] Row count after removing unproductive sequences:", nrow(dat), "\n")

## Aggregate duplicate count by VHH
tmp <- aggregate(duplicate_count ~ VHH, FUN = sum, data = dat)

cat("[tmp] Row count after aggregating to the VHH level:", nrow(tmp), "\n")

## Demand duplicate VHHs
tmp <- tmp[tmp$duplicate_count >= 2,]

cat("[tmp] Row count after demanding duplicates:", nrow(tmp), "\n")

## Use hash table for quick lookups
# Create a hash table from tmp
hash_table <- hash(tmp$VHH, TRUE)

system.time(
  # Filter dat using the hash table
  dat <- dat[sapply(dat$VHH, function(x) has.key(x, hash_table)), ]
)

# Cleanup
rm(hash_table)

cat("[Dat] Row count after hash lookups in tmp for 2+ VHHs:", nrow(dat), "\n")

dat <- dat[!duplicated(dat$sequence),]

cat("[Dat] Row count after removing redundant (duplicate) sequence observations:", nrow(dat), "\n")


system.time(
  # Now that we've saved some time, we can merge back for VHH dupe counts
  dat <- merge(dat, tmp, by = "VHH")
)

# count rows
cat("[Dat] Row count after Merge:", nrow(dat), "\n")

colnames(dat)[which(colnames(dat) == "duplicate_count.x")] <- "duplicate_count"
colnames(dat)[which(colnames(dat) == "duplicate_count.y")] <- "VHH_duplicate_count"

# Sequences are the unique identifier, 
# so we drop duplicated entries after the VHH-level Merge

# Assign Global Set IDs
dat$sequence_id <- paste0("GS-", 1:nrow(dat))

cat("Calculating distToNearest...\n")

db <- shazam::distToNearest(dat, model = "ham", 
                            normalize = "len", 
                            vCallColumn = "v_call", 
                            jCallColumn = "j_call",
                            locusColumn = "locus",
                            nproc = 1)

### Find Threshold

# Find threshold using density method
output <- shazam::findThreshold(db$dist_nearest, method="density")
threshold <- output@threshold

# Decrease threshold
new_thr <- threshold - (threshold*0.25)
cat("25% Lower Threshold:\n")
cat(round(new_thr, 2), '\n')
cat('\n')

cat("Inferring Clones:\n")

# Clonal assignment using hierarchical clustering
db <- scoper::hierarchicalClones(db,
                                 threshold=new_thr,
                                 v_call="v_call",
                                 j_call = "j_call",
                                 verbose = TRUE,
                                 summarize_clones = FALSE)


cat("Writing: ", output_file, '\n')
write.csv(db, file = output_file, row.names = F)

## ###
## Here Clean Clusters
## ###













## ###
## Using Global Set for Lookups
## ###


root_dir <- "~/Downloads/2024_TUCF_7981/"
files <- list.files(root_dir)
GS_name <- list.files(root_dir)[grepl("GS", list.files(root_dir))]
# Drop global set
files <- files[files != GS_name]



# Iterate over files (not including global data set)
data_files <- list()
for (f in files){
  cat("Loading: ", f, '\n')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)
  # Drop clone info if it is there (it should not be)
  tmp$dist_nearest <- NULL
  tmp$clone_id <- NULL
  file_id <- gsub("_min_threshold_clone-pass.csv", "", f)
  tmp$file <- rep(file_id, times = nrow(tmp))
  tmp <- merge(tmp, db[,c("sequence", "clone_id", "VHH")], by = "sequence")
  data_files[[file_id]] <- tmp
}


## ###
## Abundance Tabulation
## ###

## IDEA
## Merge all data sets to the global set, then drop the global set






















































## ###
## All dup function

all_dup <- function(vec){
  return(vec[duplicated(vec) | duplicated(vec, fromLast=TRUE)])
}



###





## 
















root_dir <- "~/Downloads/2024_TUCF_7981/"

files <- list.files(root_dir)

GS_name <- list.files(root_dir)[grepl("GS", list.files(root_dir))]

# Drop global set
files <- files[files != GS_name]

## Load Global Set
GS <- read.csv(paste0(root_dir,"GS_Schisto_d2VHH_13Aug24_min_threshold_clone-pass.csv"),
               header = T, stringsAsFactors = F)


##


test <- GS[GS$sequence %in% duped_seqs,]











# Iterate over files (not including data set)
for (f in files[1]){
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)
  # Drop clone info if it is there
  tmp$dist_nearest <- NULL
  tmp$clone_id <- NULL
  
  test <- merge(tmp, GS, by = intersect(colnames(GS), colnames(tmp)), all.x = T)
  
  
}


## Check to see if the sequences that are not in the global set have low counts
dat_in <- tmp[(tmp$sequence %in% GS$sequence),]
dat_out <- tmp[!(tmp$sequence %in% GS$sequence),]