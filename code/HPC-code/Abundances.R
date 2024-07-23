## ###
## Abundance_Tables.R
## ###

## Load Utilities
# source("~/VHH/lib/utilities.R")
#
## Load Utilities
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

## Pull Root Directory
root_dir <- "/cluster/tufts/cowenlab/emosel01/data/"

output_table <- "background_abundance_table.csv"

cat("Creating: ", paste0(root_dir, output_table),"\n", sep = '')

## Draw file info from commandArgs (for use later)
# args <- commandArgs(trailingOnly=TRUE)
# file_name <- args[1] 

## /cluster/tufts/cowenlab/emosel01/data/hclust


# Create a list to use as a dictionary
files <- list()

cat("Locating Files:\n")
## Can change this later to be file_name <- "globals/GS_d2_18Jun24-Cleaned19Apr24.csv"
cat(paste0(root_dir, "globals/GS_background_d2VHH_16Jul24-Cleaned19Apr24.csv\n"))

# Populate key/values with files and their keys
files[["master"]] <- paste0(root_dir, "globals/GS_background_d2VHH_16Jul24-Cleaned19Apr24.csv")

file_list <- list.files(paste0(root_dir, "background/"))

for (nm in gsub("_db-pass.tsv", "", file_list)){
  # Load file
  cat(paste0("background/",nm, "_db-pass.tsv\n"))
  files[[nm]] <- paste0(root_dir, paste0("background/",nm, "_db-pass.tsv"))
}

max_clone_count <- 0
max_id <- "S0"

cat("Locating Files:\n")
for (nm in names(files)){
  cat(nm, '\n')
  if (nm != "master"){
    cat("Loading: ", nm, '\n')
    files[[nm]] <- read.csv(files[[nm]], header = T, stringsAsFactors = F, sep = '\t')
  } else {
    cat("Loading: ", nm, '\n')
    files[[nm]] <- read.csv(files[[nm]], header = T, stringsAsFactors = F)
  }
  
  cat("Pre-Cleaning:\n")
#  print(head(files[[nm]]))
  clone_count <- length(unique(files[[nm]]$clone_id))
  cat(clone_count, '\n')
  cat(sum(files[[nm]]$duplicate_count), '\n')
  
  if (nm != "master"){
    # Trim CDR3s
    files[[nm]]$junction_aa <- trim_CDR3(files[[nm]]$junction_aa)
    # Translate alignments to VHHs
    files[[nm]]$VHH <- fixed_translate_DNA(files[[nm]]$sequence_alignment)
    # Trim sequence, as we do elsewhere
    files[[nm]]$VHH <- substr(files[[nm]]$VHH, 11, nchar(files[[nm]]$VHH))
    # Drop unproductive sequences
    files[[nm]] <- unproductive(files[[nm]])
  }
  
  cat("After-Cleaning:\n")
  cat(length(unique(files[[nm]]$clone_id)), '\n')
  cat(sum(files[[nm]]$duplicate_count), '\n')
  cat('\n')
  
  if (max_clone_count < clone_count){
    max_clone_count <- clone_count
    max_id <- nm
  }
}

# Subset the global set
master <- files[["master"]]

# Iterate over the files to ensure that we have VHH's accounted for
for (f in names(files)){
  cat(f, '\n')
  cat(nrow(files[[f]]), '\n')
  files[[f]] <- files[[f]][which(files[[f]]$VHH %in% master$VHH),]
  cat(nrow(files[[f]]), '\n')
}

## Assign Families from Global Set

## ###
## Note: Is there a bug here?
## We need to perform lookups by full-VHH
## As the CDR3 could have ambiguous/incourrect matches
## ###

res <- data.frame()
for (f in names(files)[names(files) != "master"]){
  # Subset data frame
  tmp <- files[[f]]
  print(f)
  #print(table(tmp$clone_id))
  # Iterate over each VHH
  count <- 0
  max_count <- length(unique(tmp$junction_aa))
  for (CDR3 in unique(tmp$junction_aa)){
    count <- count + 1
    if (round(count/max_count*100, 2) %% 5 == 0) cat("Progress:", round(count/max_count*100, 2), '\n')
    tmp_frame <- tmp[tmp$junction_aa == CDR3,]
    tmp_frame$clone_id <- rep(master[which(master$junction_aa == CDR3),]$clone_id[1], nrow(tmp_frame))
    tmp_frame$f <- rep(f, nrow(tmp_frame))
    # Rbind those
    res <- rbind(res, tmp_frame)
  }
}

cat("Writing: ", paste0(root_dir, output_table), '\n', sep = '')
write.csv(res, file = paste0(root_dir, output_table), row.names = F)

## Aggregate to clone level

dfs <- list()

for (nm in names(files)){
  cat(nm, '\n')
  dfs[[nm]] <- clonal_statistics(res[res$f == nm,])
}

# Remove clone ids
for (nm in names(files)){
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

for (i in 1:nrow(merge_order)){
  # Rename abundance
  colnames(dfs[[merge_order$run[i]]])[2] <- paste0("abundance-", merge_order$run[i])
}

tmp <- merge(dfs[[merge_order$run[1]]], 
             dfs[[merge_order$run[2]]], 
             by = "median_CDR3", all.x = T)


for (i in 3:nrow(merge_order)){
  tmp <- merge(tmp, dfs[[merge_order$run[i]]], 
               all.x = T, 
               by = "median_CDR3")
}

holder <- tmp

for (nm in colnames(tmp)[2:ncol(tmp)]){
  tmp[[nm]] <- as.numeric(tmp[[nm]])
}

# Replace NA's with zero
tmp[is.na(tmp)] <- 0

cat("Writing: ", "/cluster/tufts/cowenlab/emosel01/data/GS_background_18Jul24-Cleaned.csv", '\n')
write.csv(tmp, file = "/cluster/tufts/cowenlab/emosel01/data/GS_background_18Jul24-Cleaned.csv", row.names = F)

