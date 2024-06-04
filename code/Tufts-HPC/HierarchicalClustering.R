if (!require(airr)) install.packages(airr)
if (!require(alakazam)) install.packages(alakaza)
if (!require(tigger)) install.packages(tigger)
if (!require(scoper)) install.packages(scoper)
if (!require(shazam)) install.packages(shazam)

input_file <- "/cluster/tufts/cowenlab/emosel01/data/GS_sample.csv"
output_file <- "/cluster/tufts/cowenlab/emosel01/data/GS_sample_min_threshold_clone-pass.csv"

db <- read.csv(input_file, header = T, stringsAsFactors = F)

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


cat("Inferring Clones:\n")
system.time(
  # Clonal assignment using hierarchical clustering
  db <- hierarchicalClones(db, threshold=thr, v_call="v_call", j_call = "j_call")
)

cat("Writing: ", output_file, '\n')
write.csv(db@db, output_file)