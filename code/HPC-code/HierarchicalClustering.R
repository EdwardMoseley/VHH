## Change Library path to ensure we have access

.libPaths("/cluster/tufts/hpc/tools/R/4.3.0")

library("alakazam")
library("scoper")
library("shazam")
library("hash")
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

# Set root dir
root_dir <- "/cluster/tufts/cowenlab/emosel01/data/background/"

output_file <- "/cluster/tufts/cowenlab/emosel01/data/GS_background_d2VHH_16Jul24_min_threshold_clone-pass.csv"

files <- list.files(root_dir)[grepl(".tsv", list.files(root_dir))]

# Ignore certain files
# files <- files[!(files %in% c("C1_S1_db-pass.tsv", "C2_S2_db-pass.tsv", "D1_S3_db-pass.tsv", "D2_S4_db-pass.tsv", "D3_S5_db-pass.tsv", "D4_S6_db-pass.tsv"))]

cat("Processing files:\n")

for (f in files){
  cat(f, '\n')
}

dat <- data.frame()
for (f in files){
  cat("Loading:", f,'\n')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F, sep = '\t')
  tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
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

cat("[tmp] Row count after aggregating:", nrow(tmp), "\n")

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

cat("[Dat] Row count after hash lookups in tmp:", nrow(dat), "\n")

system.time(
  # Now that we've saved some time, we can merge back for VHH dupe counts
  dat <- merge(dat, tmp, by = "VHH")
)

# count rows
cat("[Dat] Row count after Merge:", nrow(dat), "\n")

colnames(dat)[which(colnames(dat) == "duplicate_count.x")] <- "duplicate_count"
colnames(dat)[which(colnames(dat) == "duplicate_count.y")] <- "VHH_duplicate_count"

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
