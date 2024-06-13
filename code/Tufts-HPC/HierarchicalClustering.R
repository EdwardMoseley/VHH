# Change Library path to ensure we have access
.libPaths("/cluster/tufts/hpc/tools/R/4.3.0")

if (!require("airr")) install.packages("airr")

library("alakazam")
#library("tigger")
library("scoper")
library("shazam")


root_dir <- "/cluster/tufts/cowenlab/emosel01/data/hclust"
output_file <- "/cluster/tufts/cowenlab/emosel01/data/GS_full_12Jun24_min_threshold_clone-pass.csv"
files <- list.files(root_dir)

dat <- data.frame()
for (f in files){
  cat("Loading:", f,'\n')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F, sep = '\t')
  # tmp$junction_aa <- trim_CDR3(tmp$junction_aa)
  # tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  # tmp$VHH <- substr(tmp$VHH, 17, nchar(tmp$VHH))
  # tmp <- tmp[,c("sequence_id", "junction_aa", "VHH", "duplicate_count", "file")]
  tmp <- tmp[,c("sequence", "locus", "sequence_alignment", "junction", "junction_aa", "v_call", "d_call", "j_call", "duplicate_count")]
  #tmp$file <- rep(f, times = nrow(tmp))
  dat <- rbind(dat, tmp)
}

rm(tmp)

tmp <- aggregate(duplicate_count ~ ., FUN = sum, data = dat)

cat("Writing Global Set (Clones Not Yet Inferred)\n")

write.csv(tmp, file = "/cluster/tufts/cowenlab/emosel01/data/hclust/ImpulseLuluGlobalSet12Jun24.csv", row.names = F)

dat <- tmp

rm(tmp)

## Clonal Inference

### Dist_to_nearest calculation

cat("Calculating distToNearest...\n")

db <- shazam::distToNearest(dat, model = "ham", 
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
print("25% Lower Threshold:")
print(round(new_thr, 2))

cat("Inferring Clones:\n")

  # Clonal assignment using hierarchical clustering
db <- scoper::hierarchicalClones(db,
                                   threshold=new_thr,
                                   v_call="v_call",
                                   j_call = "j_call",
                                   verbose = TRUE)



cat("Writing: ", output_file, '\n')
write.csv(db@db, file = output_file, row.names = F)