# VHH-Query
## Relies on MMseqs2
# https://github.com/soedinglab/MMseqs2

# Utilities
source("/cluster/tufts/cowenlab/emosel01/code/lib/utilities.R")

# Path to MMseqs2
mmseqs_path <- "mmseqs "

# Root directory
root_dir <- "/cluster/tufts/cowenlab/emosel01/data/"
# These directories are just children of the root
query_dir <- "/cluster/tufts/cowenlab/emosel01/data/queries/"
# query_sequences04July24.csv
output_dir <- "/cluster/tufts/cowenlab/emosel01/data/tmp/"
data_dir <- "/cluster/tufts/cowenlab/emosel01/data/downloads/"

files <- list.files(data_dir)# [grepl("_db-pass.tsv", list.files(data_dir))]
# files <- files[!grepl("demand",files)]
print(files)

# Load Query File
query_file <- list.files(query_dir)[1]
cat("Working on: ", query_file)

# Queries
queries <- read.csv(paste0(query_dir, query_file), header = T, stringsAsFactors = F)
print(head(queries))

# Write query to fasta
write_fasta(queries, "sequence", "sequence_id", paste0(query_dir, "rev_queries.fasta"))

# Create data frame for results
results <- data.frame()

# Iterate over query databases
for (f in files){
cat("Working on: ", f, '\n')
#Build command
  	command <- paste(mmseqs_path, "easy-search",
                   # Path to queries
                   paste0(query_dir, "rev_queries.fasta"),
                   # Path to database
                   paste0(data_dir, f),
                   # Results ard held in
                   paste0(output_dir, "alnResult.m8"),
                   output_dir,
                   ## Options
                   " --search-type 3 --format-mode 4 -v 2", sep = ' ')
  
  	cat("Calling command:\n")
  	cat(command, '\n')
  	system(command)

	# Load res
	tmp_frame <- read.csv(paste0(output_dir, "alnResult.m8"), header = T, stringsAsFactors = F, sep = '\t')
	# Change names quickly
	colnames(tmp_frame)[1] <- "sequence_id"
	# Order to pull top hits
	tmp_frame <- tmp_frame[order(tmp_frame$fident, decreasing = T),]
	tmp_frame$file <- rep(f, times = nrow(tmp_frame))
	tmp <- merge(queries, tmp_frame, by = "sequence_id")
	print(f)
	print(head(tmp))

	results <- rbind(results, tmp)
}

cat("Got results and writing...\n")
cat(paste0(output_dir, "query_results.csv\n"))

write.csv(results, file = paste0(output_dir, "query_results.csv"), row.names = F)

