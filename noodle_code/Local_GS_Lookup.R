## Local_GS_Lookup

root_dir <- "~/Downloads/2024_TUCF_7981/"

files <- list.files(root_dir)

GS_name <- list.files(root_dir)[grepl("GS", list.files(root_dir))]

# Drop global set
files <- files[files != GS_name]

## Load Global Set
GS <- read.csv(paste0(root_dir,"GS_Schisto_d2VHH_13Aug24_min_threshold_clone-pass.csv"),
               header = T, stringsAsFactors = F)

results <- data.frame()
for (f in files[1]){
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)
  # Drop clone info if it is there
  tmp$dist_nearest <- NULL
  tmp$clone_id <- NULL
  
  tmp <- merge(tmp, GS, by = intersect(colnames(GS), colnames(tmp)), all.x = T)
  
  
}


## Check to see if the sequences that are not in the global set have low counts
dat_in <- tmp[(tmp$sequence %in% GS$sequence),]
dat_out <- tmp[!(tmp$sequence %in% GS$sequence),]