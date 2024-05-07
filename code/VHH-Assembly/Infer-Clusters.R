## Productivity and Infer Clusters

## Load Utilities
source("~/VHH/lib/utilities.R")

## Draw file info from commandArgs
args <- commandArgs(trailingOnly=TRUE)

file_name <- args[1]

cat("Processing: ", file_name, '\n', sep = '')
# Define root_path
root_path <- "~/VHH/data/"
# Load Data
dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)



# files <- list.files(data_dir)[grepl(".tsv", list.files(data_dir))]

## Clonal Inference

cat("Calculating distToNearest...\n")

db <- shazam::distToNearest(dat, model = "ham", 
                            normalize = "len", 
                            vCallColumn = "v_call", 
                            jCallColumn = "j_call",
                            nproc = 4)

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

write.csv(db, file = paste(root_path, file_name, sep = ''), row.names = F)