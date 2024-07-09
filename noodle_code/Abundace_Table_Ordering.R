## We will load global set
## We will also load the abundance table

## We will pull the most abundant VHH for each clonal family
# Then we will calculate their distance matrix
# Then we will order them accordingly

library(stringdist)

dat <- read.csv("~/Downloads/GS_NGO1549_d2VHH_09Jul24-Cleaned19Apr24.csv")
tmp <- read.csv("~/Downloads/NGO1549_Abundance_Table.csv", header = T, stringsAsFactors = F)

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

main_title <- paste("Cluster Dendrogram of CDR3 families from NGO1549 Abundance Data\nDemand >= ", 1," VHHs, Distance Metric: Levenshtein", sep = '')

par(family = "Courier", cex = 1)

# Plot the dendrogram with labels
plot(cluster_tree, 
     main = main_title, 
     xlab = "Overall VHH Count-Mode CDR3", sub = NULL, 
          labels = tmp$junction_aa,#[cluster_tree$order]), 
     ylab = "Distance (Levenshtein)", 
     cex = 0.70)

write.csv(tmp[cluster_tree$order,], file = "~/Downloads/NGO1549_Dendrogram_Ordered_09Jul24.csv", row.names = F)

