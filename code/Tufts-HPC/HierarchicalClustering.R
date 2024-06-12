#if (!require(airr)) install.packages(airr)
#if (!require(alakazam)) install.packages(alakaza)
#if (!require(tigger)) install.packages(tigger)
#if (!require(scoper)) install.packages(scoper)
#if (!require(shazam)) install.packages(shazam)

# Change Library path to ensure we have access
.libPaths("/cluster/tufts/hpc/tools/R/4.3.0")

if (!require("airr")) install.packages("airr")

library("alakazam")
#library("tigger")
library("scoper")
library("shazam")


input_file <- "/cluster/tufts/cowenlab/emosel01/data/GS_sample.csv"
output_file <- "/cluster/tufts/cowenlab/emosel01/data/GS_sample_min_threshold_clone-pass.csv"

db <- read.csv(input_file, header = T, stringsAsFactors = F)

db$locus <- rep("IGH", times = nrow(db))

cat("Calculating distToNearest...\n")

db <- distToNearest(db, model = "ham", normalize = "len",
                    vCallColumn = "v_call",
                    jCallColumn = "j_call",
                    locusColumn = "locus",
                    nproc = 1)

# determine the threshold
threshold <- findThreshold(db$dist_nearest, method = "density")
thr <- threshold@threshold

cat("Threshold Data:\n")
cat(thr,'\n')
new_thr <- thr - (thr*0.25)
cat("25% Lower Threshold:\n")
cat(round(new_thr, 2), '\n')


cat("Inferring Clones:\n")
# Clonal assignment using hierarchical clustering
db <- hierarchicalClones(db, threshold=thr, v_call="v_call", j_call = "j_call")


cat("Writing: ", output_file, '\n')
write.csv(db@db, file = output_file, row.names = F)