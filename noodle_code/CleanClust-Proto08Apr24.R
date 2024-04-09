## ###
## CleanClust-Prototype
## 08Apr24
## ###

# Load utilities
source("~/VHH/lib/utilities.R")

root_path <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/"
file_name <- "18_VecF3_R4_S18_min_threshold_clone-pass.csv"

dat <- read.csv(paste(root_path, file_name, sep = ''), header = T, stringsAsFactors = F)

  
# 1. `~/VHH/data/Cluster-DBs/phylogenetic_dbs/19_VecF1_R3_S19_min_threshold_clone-pass.csv`
# 2. `~/VHH/data/Cluster-DBs/phylogenetic_dbs/20_VecF2_R3_S20_min_threshold_clone-pass.csv`
# 3. `~/VHH/data/Cluster-DBs/phylogenetic_dbs/18_VecF3_R4_S18_min_threshold_clone-pass.csv`
# 4. `~/VHH/data/Cluster-DBs/phylogenetic_dbs/21_VecF4_R3_S21_min_threshold_clone-pass.csv`

