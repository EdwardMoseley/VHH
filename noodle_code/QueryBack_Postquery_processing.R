source("~/VHH/lib/utilities.R")
library("writexl")

root_dir <- "~/Downloads/"
input_file <- "query_results04Aug24.csv"
query_sequences <- "query_seqs01Aug24.csv"
database_folder <- "query_back/"
#output_file <-"A_S57_min_threshold_clone-pass.csv"

tmp <- read.csv(paste0(root_dir, input_file), header = T, stringsAsFactors = F)

colnames(tmp)[1] <- "query_identifier"
colnames(tmp)[2] <- "sequence_id"


tmp <- tmp[grepl(".fasta", tmp$file),]
tmp$file <- gsub("_db-pass.fasta", "", tmp$file)

# boxplot(tmp$fident ~ tmp$file, las = 3, ylim = c(0, 1))

# db <- read.csv(paste0(root_dir, database_file), header = T, stringsAsFactors = F, sep = '\t')

suffix <- "_db-pass.tsv"

# Build results data frame
results <- data.frame()

for (f in unique(tmp$file)){
  cat(f, '\n')
  db <- read.csv(paste0(root_dir, database_folder, f, suffix), header = T, stringsAsFactors = F, sep = '\t')
  for (query_seq in unique(tmp$query_identifier)){
    cat(query_seq, '\n')
    tmp_frame <- merge(tmp[tmp$file == f & tmp$query_identifier == query_seq,], db, by = "sequence_id")
    results <- rbind(results, tmp_frame)
  }
}

results$VHH <- fixed_translate_DNA(results$sequence_alignment)
results$VHH <- substr(results$VHH, 11, nchar(results$VHH))

query_db <- read.csv(paste0(root_dir, query_sequences), header = T, stringsAsFactors = F)

dat <- merge(results, query_db, by = "query_identifier")

dat$junction_aa <- trim_CDR3(dat$junction_aa)

keep_cols <- c("query_identifier",
               "sequence_id",
               "junction_aa",
               "duplicate_count",
               "file",
               "fident",
               "evalue",
               "VHH",
               "Query_VHH")


for (nm in unique(dat$query_identifier)){
  cat(nm,'\n')
  tmp_frame <- dat[dat$query_identifier == nm, ]
  tmp_frame <- tmp_frame[order(tmp_frame$fident, decreasing = T),]
  tmp_frame <- tmp_frame[!duplicated(tmp_frame),]
  tmp_frame <- tmp_frame[tmp_frame$fident >= 0.7,]
  #write.csv(tmp_frame[,keep_cols], paste0(root_dir, nm, "_query_results01Aug24.csv"), row.names = F)
  
  
  # Write the data frame to an Excel file
  write_xlsx(tmp_frame[,keep_cols], paste0(root_dir, nm, "_query_results04Aug24.xlsx"))
}


