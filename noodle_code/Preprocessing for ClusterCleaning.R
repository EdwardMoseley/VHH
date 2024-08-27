
source("~/VHH/lib/utilities.R")

tmp <- read.csv("~/VHH/data/Cluster-DBs/phylogenetic_dbs/Global_Clone_Experiment_master_set-D1-09-10-11_min_threshold_clone-pass.csv", header = T, stringsAsFactors = F)

tmp$junction_aa <- translate_DNA(tmp$junction)

write.csv(tmp, file = "~/VHH/data/Cluster-DBs/phylogenetic_dbs/Global_Clone_Experiment_master_set-D1-09-10-11_min_threshold_clone-pass.csv", row.names = F)

for (i in 1:nrow(tmp)){
  if (tmp$junction_aa[i] != tmp$cdr3_aa[i]){
    print(tmp[i,])
  }
}

tmp <- read.csv("~/VHH/data/Cluster-DBs/phylogenetic_dbs/Global_Clone_Experiment_master_set-D1-09-10-11-rule-pooled-families.csv", header = T, stringsAsFactors = F)

tmp$sequence_id <- paste0(1:nrow(tmp), '-',tmp$junction_aa)

tmp <- tmp[,c("sequence_id", colnames(tmp)[1:(ncol(tmp)-1)])]

write.csv(tmp, file = "~/VHH/data/Cluster-DBs/phylogenetic_dbs/Global_Clone_Experiment_master_set-D1-09-10-11-rule-pooled-families.csv", row.names = F)


# Rscript CleanCluster-Pool.R Global_Clone_Experiment_9_S18_demand_two_min_threshold_clone-pass.csv