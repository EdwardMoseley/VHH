# Utilities
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

# Define path to query data
query_file <- "/cluster/tufts/cowenlab/emosel01/data/tmp/query_results04Jul24.csv"

# Data Dir
data_dir <- "/cluster/tufts/cowenlab/emosel01/data/tmp_holder/"


# Note that these files are .tsvs
files <- list.files(data_dir)

# Load query Results
dat <- read.csv(query_file, header = T, stringsAsFactors = F)
# Adjust Column names
colnames(dat)[1] <- "query_sequence"
colnames(dat)[3] <- "sequence_id"

print(head(dat))

cat("Iterating over database files...")
results <- data.frame()
for (f in files){
	cat("Loading: ", paste0(data_dir, f), '\n')
	tmp <- read.csv(paste0(data_dir, f), header = T, stringsAsFactors = F, sep = '\t')
	tmp$file <- rep(gsub(".tsv", ".fasta", f), times = nrow(tmp))
	print(head(tmp))
	# Now merge to our query data
	tmp_frame <- merge(dat, tmp, by = c("sequence_id", "file"))
	# Add to results
	results <- rbind(results, tmp_frame)
}


cat("Writing: ", "/cluster/tufts/cowenlab/emosel01/data/joined_query_results04Jul24.csv\n", sep = '')
write.csv(results, file = "/cluster/tufts/cowenlab/emosel01/data/joined_query_results04Jul24.csv", row.names = F)

cat("Done Writing.\n")

