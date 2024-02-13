suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tigger))
suppressPackageStartupMessages(library(scoper))
# get the distance to nearest neighbors
suppressPackageStartupMessages(library(shazam))

directory <- "~/VHH/data"
files <- list.files(directory)[grepl("db-pass.tsv", list.files(directory))]

done_list <- c("~/VHH/data/17_VecF3_R3_S17_db-pass.tsv",
               "~/VHH/data/18_VecF3_R4_S18_db-pass.tsv",
               "~/VHH/data/20_VecF2_R3_S20_db-pass.tsv",
               "~/VHH/data/22_VecF4_R4_S22_db-pass.tsv",
               "~/VHH/data/Selected-XLN-XME-XMF-SmPDE-VHHs_9-23_db-pass.tsv",
               "~/VHH/data/13_VecF1_R2_S13_atleast-1_db-pass.tsv",
               #
               "~/VHH/data/1_L1_sh2_S1_db-pass.tsv", 
               "~/VHH/data/10_L3_lh2_S10_db-pass.tsv", 
               "~/VHH/data/11_L3_sh3_S11_db-pass.tsv", 
               "~/VHH/data/13_VecF1_R2_S13_db-pass.tsv",
               "~/VHH/data/14_VecF1_R2_S14_db-pass.tsv")


file_list <- paste(directory, files, sep = '/')
file_list <- file_list[!(file_list %in% done_list)]


# Drop unproductive sequences
for (f in file_list){
  # Read the dataset
  db <- read_rearrangement(file.path(f))
  # Subset productive clones
  cat("Input:", nrow(db),'\n')
  db <- db[db$productive == T,]
  cat("Productive:", nrow(db), '\n')
  
  file_out <- gsub("_db-", "_productive-", f)
  cat('\n')
  cat("Writing: ", f, '\n')
  write.table(db, file = file_out, sep = '\t')
  
  cat("Calculating distToNearest...\n")
  
  db <- distToNearest(db, model = "ham", normalize = "len", 
                      vCallColumn = "v_call", 
                      jCallColumn = "j_call",
                      nproc = 4)
  
  # determine the threshold
  threshold <- findThreshold(db$dist_nearest, method = "density")
  thr <- threshold@threshold
  
  print("Threshold Data:")
  print(thr)
  new_thr <- thr - (thr*0.25)
  print("25% Lower Threshold:")
  print(round(new_thr, 2))
  
  cat('\n')
  
  # Outfile prefix
  prefix <- gsub("_db-pass.tsv", '', gsub("~/VHH/data/", "", f))
  prefix <- paste(prefix, "_min_threshold", sep = '')
  
  
  cat("Inferring Clones:\n")
  system.time(
  # Clonal assignment using hierarchical clustering
  db <- hierarchicalClones(db, threshold=thr, v_call="v_call", j_call = "j_call")
  )
  
  cat("Writing: ", gsub("_productive-pass.tsv", "_min_threshold_clone-pass.tsv", file_out), '\n')
  write.table(db@db, file = gsub("_productive-pass.tsv", "_min_threshold_clone-pass.tsv", file_out), sep = '\t')
  
}