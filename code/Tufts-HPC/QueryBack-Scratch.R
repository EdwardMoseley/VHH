#QueryBack-Scratch

root_dir <- "~/Downloads/"
input_file <- "query_results23Jul24.csv"
database_file <- "query_back/A_S57_db-pass.tsv"
output_file <-"A_S57_min_threshold_clone-pass.csv"

tmp <- read.csv(paste0(root_dir, input_file), header = T, stringsAsFactors = F)

tmp <- tmp[grepl(".fasta", tmp$file),]
tmp$file <- gsub("_db-pass.fasta", "", tmp$file)

boxplot(tmp$fident ~ tmp$file, las = 3, ylim = c(0, 1))

db <- read.csv(paste0(root_dir, database_file), header = T, stringsAsFactors = F, sep = '\t')









library("alakazam")
library("scoper")
library("shazam")





## To do
# Hclust one of the data sets so we can pull nearest families

cat("Calculating distToNearest...\n")

db <- shazam::distToNearest(db, model = "ham",
                            normalize = "len",
                            vCallColumn = "v_call",
                            jCallColumn = "j_call",
                            locusColumn = "locus",
                            nproc = 1)

### Find Threshold

# Find threshold using density method
output <- shazam::findThreshold(db$dist_nearest, method="density")
threshold <- output@threshold

# Decrease threshold
new_thr <- threshold - (threshold*0.25)
cat("25% Lower Threshold:\n")
cat(round(new_thr, 2), '\n')
cat('\n')

cat("Inferring Clones:\n")

# Clonal assignment using hierarchical clustering
db <- scoper::hierarchicalClones(db,
                                 threshold=new_thr,
                                 v_call="v_call",
                                 j_call = "j_call",
                                 verbose = TRUE,
                                 summarize_clones = FALSE)



cat("Writing: ", output_file, '\n')
write.csv(db, file = paste0(root_dir, output_file), row.names = F)