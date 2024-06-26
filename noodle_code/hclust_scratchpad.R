## Change Library path to ensure we have access

# .libPaths("/cluster/tufts/hpc/tools/R/4.3.0")

#if (!require("airr")) install.packages("airr")

library("alakazam")
library("scoper")
library("shazam")
library("hash")


#root_dir <- "/cluster/tufts/cowenlab/emosel01/data/hclust/"
root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs"

#output_file <- "/cluster/tufts/cowenlab/emosel01/data/GS_ImpulseLulu_20Jun24_min_threshold_clone-pass.csv"
files <- list.files(root_dir)
# We only want files with underscores
files <- files[grepl("_min_threshold_clone-pass.csv", files)]
## We only want to pull a couple
files <- files[files %in% c("19_VecF1_R3_S19_min_threshold_clone-pass.csv", "20_VecF2_R3_S20_min_threshold_clone-pass.csv", "21_VecF4_R3_S21_min_threshold_clone-pass.csv")]


dat <- data.frame()
for (f in files){
  cat("Loading:", f,'\n')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)#, sep = '\t')
  tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  tmp <- tmp[,c("sequence", "locus", "sequence_alignment", "junction", "junction_aa", "v_call", "d_call", "j_call", "VHH", "duplicate_count")]
  dat <- rbind(dat, tmp)
}


## Aggregate duplicate count by VHH
tmp <- aggregate(duplicate_count ~ VHH, FUN = sum, data = dat)

## Demand duplicate VHHs
tmp <- tmp[tmp$duplicate_count >= 2,]

## Use hash table for quick lookups

# Create a hash table from tmp
hash_table <- hash(tmp$VHH, TRUE)

# Filter dat using the hash table
dat <- dat[sapply(dat$VHH, function(x) has.key(x, hash_table)), ]

# Now that we've saved some time, we can merge back for VHH dupe counts
dat <- merge(dat, tmp, by = "VHH")

# count rows
cat("Row count:", nrow(dat), "\n")

rm(tmp)

## Clonal Inference

### Dist_to_nearest calculation

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