## We will load global set
## We will also load the abundance table

## We will pull the most abundant VHH for each clonal family
# Then we will calculate their distance matrix
# Then we will order them accordingly

library(stringdist)

dat <- read.csv("~/Downloads/GS_rp_vs_ch2_d2VHH_24Jul24-Cleaned19Apr24.csv", header = T, stringsAsFactors = F)
# dat <- read.csv("~/Downloads/GS_TbpA_d2VHH_14Jul24-Cleaned19Apr24.csv", header = T, stringsAsFactors = F)
# tmp <- read.csv("~/Downloads/Lst_Abundance_Table.csv", header = T, stringsAsFactors = F)

# tmp <- read.csv("~/Downloads/GS_uncoated_17Jul24-Cleaned.csv", header = T, stringsAsFactors = F)
#tmp <- read.csv("~/Downloads/GS_TbpA_18Jul24-Cleaned.csv", header = T, stringsAsFactors = F)
#tmp <- read.csv("~/Downloads/GS_background_18Jul24-Cleaned.csv", header = T, stringsAsFactors = F)

tmp <- read.csv("~/Downloads/GS_rp_vs_ch2_24Jul24-Cleaned.csv", header = T, stringsAsFactors = F)


## Rename abundance table according to AIRR-Seq standards
colnames(tmp)[1] <- "junction_aa"

## Merge
tmp <- merge(tmp, dat[,c("junction_aa","VHH")], by = "junction_aa", all.x = T)
## Remove deplicates resulting from the cross product
tmp <- tmp[!duplicated(tmp$junction_aa),]
## Drop the * obs
tmp <- tmp[!grepl("\\*",tmp$junction_aa),]
# Now create dendrogram ordering

strings <- tmp$VHH

CDR3s <- tmp$junction_aa

# Compute the distance matrix
dist_matrix <- stringdistmatrix(strings, method = "lv")

# Perform hierarchical clustering
cluster_tree <- hclust(as.dist(dist_matrix), method = "complete")

# test <- tmp[cluster_tree$order,]

main_title <- paste("Cluster Dendrogram of CDR3 families from RP vs. CH2 (Primer) Abundance Data\nDemand >= ", 2," VHHs, Distance Metric: Levenshtein", sep = '')

par(family = "Courier", cex = 0.40)#0.55)

# Plot the dendrogram with labels
plot(cluster_tree, 
     main = main_title, 
     xlab = "Overall VHH Count-Mode CDR3", sub = NULL, 
          labels = tmp$junction_aa,#[cluster_tree$order]), 
     ylab = "Distance (Levenshtein)", 
     cex = 0.40)#0.55)

# write.csv(tmp[cluster_tree$order,], file = "~/Downloads/rp_vs_ch2_Dendrogram_Ordered_24Jul24.csv", row.names = F)

