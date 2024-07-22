## Load Utilities
source("~/VHH/lib/utilities.R")

# res <- read.csv("~/Downloads/TbpA_abundance_table.csv", header = T, stringsAsFactors = F)
# res <- read.csv("~/Downloads/TdjJ_abundance_table.csv", header = T, stringsAsFactors = F)

res <-  read.csv("~/Downloads/uncoated_abundance_table.csv", header = T, stringsAsFactors = F)


files <- unique(res$f)

## Aggregate to clone level

#dfs <- list()

#for (nm in files){
#  cat(nm, '\n')
  #dfs[[nm]] <- clonal_statistics(res[res$f == nm,])
#}


dfs <- list()
# results_out <- data.frame()
for (nm in files){
  results <- data.frame()
  cat("File: ", nm, '\n', sep = '')
  tmp_frame <- res[res$f == nm,]
  for (clone_id in unique(tmp_frame$clone_id)){
    if (!is.na(clone_id)){
      tmp_frame_two <- tmp_frame[tmp_frame$clone_id == clone_id,]
      cat("Clone ID: ", clone_id, '\n', sep = '')
      tmp_tab <- table(tmp_frame_two$junction_aa)[order(table(tmp_frame_two$junction_aa), decreasing = T)]
      median_CDR3 <- names(tmp_tab)[1]
      abundance <- sum(tmp_tab)
      results <- rbind(results, cbind(clone_id, median_CDR3, abundance))
    }
  }
  dfs[[nm]] <- results
}


# Remove clone ids
for (nm in files){
  dfs[[nm]]$clone_id <- NULL
  # We don't need clone_id in this case
}

# Drop global set
dfs[["master"]] <- NULL

# Merge order table
merge_order <- data.frame()

for (nm in names(dfs)){
  cat(nm, '\n')
  cat(nrow(dfs[[nm]]), '\n')
  cat(length(unique(dfs[[nm]]$median_CDR3)), '\n')
  # Look at counts
  # print(as.numeric(sum(dfs[[nm]]$abundance)))
  
  merge_order <- rbind(merge_order, cbind(nm, length(unique(dfs[[nm]]$median_CDR3))))
}

merge_order <- as.data.frame(merge_order)
colnames(merge_order) <- c("run", "counts")
merge_order$counts <- as.numeric(as.character(merge_order$counts))
merge_order <- merge_order[order(merge_order$counts, decreasing = T),]

#for (i in 1:nrow(merge_order)){
  # Rename abundance
#  colnames(dfs[[merge_order$run[i]]])[2] <- paste0("abundance-", merge_order$run[i])
#}

tmp <- merge(dfs[[merge_order$run[1]]],
             dfs[[merge_order$run[2]]],
             by = "median_CDR3", all.x = T, suffixes = c(merge_order$run[1], merge_order$run[2]))


for (i in 3:nrow(merge_order)){
  tmp <- merge(tmp, dfs[[merge_order$run[i]]],
               all.x = T,
               by = "median_CDR3", suffixes = c("", merge_order$run[i]))
}

colnames(tmp) <- gsub("abundance", "", colnames(tmp))

colnames(tmp)[which(colnames(tmp) == "")] <- merge_order[which(colnames(tmp) == "")-1,]$run



holder <- tmp

for (nm in colnames(tmp)[2:ncol(tmp)]){
  tmp[[nm]] <- as.numeric(tmp[[nm]])
}

# Replace NA's with zero
tmp[is.na(tmp)] <- 0



cat("Writing: ", "/cluster/tufts/cowenlab/emosel01/data/GS_uncoated_17Jul24-Cleaned.csv", '\n')
# write.csv(tmp, file = "/cluster/tufts/cowenlab/emosel01/data/GS_TbpA_15Jul24-Cleaned.csv", row.names = F)

write.csv(tmp, file = "~/Downloads/GS_uncoated_17Jul24-Cleaned.csv", row.names = F)