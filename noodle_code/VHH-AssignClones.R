# Assign Clones

## Strategy
# 1. Load Master Data set
# 2. Load Individual Sample Data Sets
# 3. For each VHH in the sample sets
#    * Locate in master set, and assign clone id


root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"


global_set <- read.csv(paste0(root_dir, "GS18May24-Cleaned19Apr24.csv"), header = T, stringsAsFactors = F)

files <- list.files(root_dir)[grepl("Cleaned15May24", list.files(root_dir))]

## Note: Better strategy to merge by VHH?


results <- data.frame()
for (f in files){
  cat("Loading: ", f, '\n')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)
  tmp$database <- rep(f, times = nrow(tmp))
  for (VHH in unique(tmp$VHH)){
    if (VHH %in% global_set$VHH){
      tmp_frame <- tmp[tmp$VHH == VHH,]
      tmp_frame$clone_id <- rep(global_set[global_set$VHH == VHH,]$clone_id[1], times = nrow(tmp_frame))
      tmp_frame$centroid <- rep(global_set[global_set$VHH == VHH,]$centroid[1], times = nrow(tmp_frame))
      tmp_frame$centroid_id <- rep(global_set[global_set$VHH == VHH,]$centroid_id[1], times = nrow(tmp_frame))
    }
  results <- rbind(results, tmp_frame)
  }
}


dat <- data.frame()
for (f in files){
  cat("Loading: ", f, '\n')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)
  tmp$database <- rep(f, times = nrow(tmp))
  dat <- rbind(dat, tmp)
}


counter <- 1
results <- data.frame()
for (VHH in unique(global_set$VHH)){
  if (counter %% 1000 == 0) cat("Progress: ", round(counter/nrow(global_set), 2)*100)
  tmp_frame <- dat[dat$VHH == VHH,]
  tmp_frame$clone_id <- rep(global_set[global_set$VHH == VHH,]$clone_id[1], times = nrow(tmp_frame))
  results <- rbind(results, tmp_frame)
  counter <- counter + 1
}














test <- merge(dat, global_set[,c("VHH", "clone_id")], by = "VHH", all.x = T, suffixes = c("", "_GS"))