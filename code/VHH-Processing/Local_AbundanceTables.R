## Local_GS_Lookup

source("~/VHH/lib/utilities.R")

root_dir <- "~/Downloads/2024_TUCF_7981/"
file_name <- "GS_Schisto_d2VHH_14Aug24-cleaned-GS.csv"

files <- list.files(root_dir)

GS_name <- list.files(root_dir)[grepl("GS", list.files(root_dir))]

# Drop global set
files <- files[!(files %in% GS_name)]

## Load Global Set
GS <- read.csv(paste0(root_dir,file_name),
               header = T, stringsAsFactors = F)


## Strategy:

# 
# Merge experiment_data to global_data on sequence_id for clone_ids


## Note:
# If they're in the Global Set, then they MUST
# Have >= 2 VHH observations across multiple sets (no singleton sequences)
# Not have a premature stop codon
# As a result, we can expect to lose a lot of samples from each component
# data set.


results <- list()
for (f in files){
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)
  # Drop clone info if it is there
  tmp$dist_nearest <- NULL
  tmp$clone_id <- NULL
  
  tmp <- merge(GS[,c("sequence", "clone_id")], tmp, by = "sequence", all.x = T)
  tmp <- tmp[!is.na(tmp$duplicate_count),]
  tmp$file <- rep(f, times = nrow(tmp))
  results[[f]] <- tmp
}


for (f in names(results)){
  cat(f, '\n')
  cat(length(unique(results[[f]]$clone_id)), '\n')
}


# Aggregate each experiment's duplicate_count to the clone_id level

dat <- results

for (nm in names(dat)){
  dat[[nm]] <- aggregate(duplicate_count ~ clone_id, FUN = sum, data = dat[[nm]])
}

# Merge duplicate counts on clone id to GS_agg, GS will have CDR3 mode identifier

GS_agg <- data.frame()
# Iterate over clone ids
for (nm in unique(GS$clone_id)){
  tmp_frame <- GS[GS$clone_id == nm,]
  # Unroll
  unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
  unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
  # Pull statistics
  median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
  median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
  abundance <- sum(tmp_frame$duplicate_count)
  
  # Keep the integer as the cluster_id
  clone_id <- nm
  
  GS_agg <- rbind(GS_agg, cbind(clone_id, median_CDR3, median_VHH, abundance))
}


holder <- GS_agg

for (nm in names(dat)){
  colnames(dat[[nm]])[2] <- gsub("_min_threshold_clone-pass.csv", "", nm)
  GS_agg <- merge(GS_agg, dat[[nm]], by = "clone_id", all.x = T)
}

GS_agg[is.na(GS_agg)] <- 0

tmp <- GS_agg

strings <- tmp$median_VHH

CDR3s <- tmp$median_CDR3

# Compute the distance matrix
dist_matrix <- stringdist::stringdistmatrix(strings, method = "lv")

# Perform hierarchical clustering
cluster_tree <- hclust(as.dist(dist_matrix), method = "complete")

tmp <- tmp[cluster_tree$order,]

main_title <- paste("Cluster Dendrogram of CDR3 families Abundance Data\nDemand >= ", 2," VHHs, Distance Metric: Levenshtein", sep = '')

par(family = "Courier", cex = 0.75)#0.55)

# Plot the dendrogram with labels
plot(cluster_tree, 
     main = main_title, 
     xlab = "Overall VHH Count-Mode CDR3", sub = NULL, 
     labels = tmp$median_CDR3,#[cluster_tree$order]), 
     ylab = "Distance (Levenshtein)", 
     cex = 0.75)#0.55)

write.csv(tmp, file = "~/Downloads/Schistosome_Clonal_Family_Abundances.csv", row.names = F)
