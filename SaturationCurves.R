# Saturation Curves
source("~/VHH/lib/utilities.R")

#root_dir <- "/Users/emoseley/VHH/data/Assembled/"
root_dir <- "/Users/emoseley/Downloads/VHH Assembly 24Jul24 Materials/2024_TUCF_7981"

file_list <- list.files(root_dir)[grepl("*_db-pass.tsv", list.files(root_dir))]

results <- data.frame()
for (f in file_list){
  cat("Working on ", f, "\n", sep = '')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F, sep = '\t')
  
  tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  
  df <- tmp[,c("VHH", "duplicate_count")]
  
  # Expand the data frame based on duplicate_count
  expanded_df <- df[rep(1:nrow(df), df$duplicate_count), "VHH", drop = FALSE]
  
  # Initialize result vector
  sample_sizes <- seq(1000, nrow(expanded_df), by = 1000)
  unique_counts <- numeric(length(sample_sizes))
  
  
  for (i in seq_along(sample_sizes)) {
    set.seed(i)
    sample_size <- sample_sizes[i]
    sampled_junctions <- sample(expanded_df$VHH, sample_size, replace = TRUE)
    unique_counts[i] <- length(unique(sampled_junctions))
  }
  
  # Combine results into a data frame
  saturation_results <- data.frame(
    sample_size = sample_sizes,
    unique_count = unique_counts
  )
  
  saturation_results$file <- rep(f, times = nrow(saturation_results))
  results <- rbind(results, saturation_results)
}


# results$database <- ifelse(results$file %in% c("C2_S2_db-pass.tsv", "D2_S4_db-pass.tsv", "D3_S5_db-pass.tsv", "D4_S6_db-pass.tsv"), "Source", "Enriched")

library(ggrepel)




selection <- c("C2_S2_db-pass.tsv", 
               "D2_S4_db-pass.tsv", 
               "D3_S5_db-pass.tsv", 
               "D4_S6_db-pass.tsv",
               "37_S11_db-pass.tsv",
               "38_S12_db-pass.tsv",
               "40_S15_db-pass.tsv",
               "42_S17_db-pass.tsv")



tmp <- results[results$file %in% selection,]

ggplot(tmp, aes(x=sample_size, y=unique_count, group=file, color=file)) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("SARS-CoV-2 Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Unique Count") 















































# Saturation Curves
source("~/VHH/lib/utilities.R")

#root_dir <- "/Users/emoseley/VHH/data/Assembled/"
root_dir <- "/Users/emoseley/Downloads/VHH Assembly 24Jul24 Materials/2024_TUCF_7981/"

file_list <- list.files(root_dir)[grepl(".csv", list.files(root_dir))]


results <- data.frame()
for (f in file_list[1:2]){
  cat("Working on ", f, "\n", sep = '')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)#, sep = '\t')
  
  tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  
  df <- tmp[,c("VHH", "duplicate_count")]
  
  # Expand the data frame based on duplicate_count
  expanded_df <- df[rep(1:nrow(df), df$duplicate_count), "VHH", drop = FALSE]
  
  # Initialize result vector
  sample_sizes <- seq(1000, nrow(expanded_df), by = 1000)
  unique_counts <- numeric(length(sample_sizes))
  
  
  for (i in seq_along(sample_sizes)) {
    set.seed(i)
    sample_size <- sample_sizes[i]
    sampled_junctions <- sample(expanded_df$VHH, sample_size, replace = TRUE)
    unique_counts[i] <- length(unique(sampled_junctions))
  }
  
  # Combine results into a data frame
  saturation_results <- data.frame(
    sample_size = sample_sizes,
    unique_count = unique_counts
  )
  
  saturation_results$file <- rep(f, times = nrow(saturation_results))
  results <- rbind(results, saturation_results)
}


### #####


results <- data.frame()
for (f in file_list){
  cat("Working on ", f, "\n", sep = '')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)#, sep = '\t')
  
  tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
  
  df <- tmp[,c("VHH", "duplicate_count")]
  
  # Expand the data frame based on duplicate_count
  expanded_df <- df[rep(1:nrow(df), df$duplicate_count), "VHH", drop = FALSE]
  
  # Determine the sample sizes as 1% increments
  total_rows <- nrow(expanded_df)
  sample_sizes <- seq(0.01, 1, by = 0.01) * total_rows
  sample_sizes <- ceiling(sample_sizes)  # Round up to ensure integer values
  
  unique_counts <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    set.seed(i)
    sample_size <- sample_sizes[i]
    sampled_junctions <- sample(expanded_df$VHH, sample_size, replace = TRUE)
    unique_counts[i] <- length(unique(sampled_junctions))
  }
  
  # Combine results into a data frame
  saturation_results <- data.frame(
    sample_size = sample_sizes,
    unique_count = unique_counts
  )
  
  saturation_results$file <- rep(f, times = nrow(saturation_results))
  results <- rbind(results, saturation_results)
}

tmp <- results


tmp$database <- ifelse(grepl("VecF", tmp$file), "Vectored", "Source")

tmp$file <- gsub("_min_threshold_clone-pass.csv", "", tmp$file)

ggplot(tmp, aes(x=sample_size, y=unique_count, group=file, color=file)) +
  #geom_line(aes(linetype = file), size = 1.2) +
  #geom_line(size = 1.2) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("Schistosome Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Unique Count") 

















































# Saturation Curves
source("~/VHH/lib/utilities.R")

#root_dir <- "/Users/emoseley/VHH/data/Assembled/"
root_dir <- "/Users/emoseley/Downloads/VHH Assembly 24Jul24 Materials/2024_TUCF_7981/"

file_list <- list.files(root_dir)[grepl(".csv", list.files(root_dir))]


dat <- data.frame()
results <- data.frame()
for (f in file_list){
  cat("Working on ", f, "\n", sep = '')
  tmp <- read.csv(paste0(root_dir, f), header = T, stringsAsFactors = F)#, sep = '\t')
  tmp$VHH <- fixed_translate_DNA(tmp$sequence_alignment)
  tmp$VHH <- substr(tmp$VHH, 11, nchar(tmp$VHH))
  ## Keep only productive
  tmp <- unproductive(tmp)
  ## Making Saturation Curve
  df <- tmp[,c("VHH", "duplicate_count")]
  
  # Expand the data frame based on duplicate_count
  expanded_df <- df[rep(1:nrow(df), df$duplicate_count), "VHH", drop = FALSE]
  
  # Determine the sample sizes as 1% increments
  total_rows <- nrow(expanded_df)
  sample_sizes <- seq(0.01, 1, by = 0.01) * total_rows
  sample_sizes <- ceiling(sample_sizes)  # Round up to ensure integer values
  
  unique_counts <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    set.seed(i)
    sample_size <- sample_sizes[i]
    sampled_junctions <- sample(expanded_df$VHH, sample_size, replace = TRUE)
    unique_counts[i] <- length(unique(sampled_junctions))
  }
  
  # Combine results into a data frame
  saturation_results <- data.frame(
    sample_size = sample_sizes,
    unique_count = unique_counts
  )
  
  saturation_results$file <- rep(f, times = nrow(saturation_results))
  
  results <- rbind(results, saturation_results)
  
  Seqs <- length(unique(tmp$sequence_id))
  SeqAligns <- length(unique(tmp$sequence_alignment))
  VHHs <- length(unique(tmp$VHH))
  ClonalFamilies <- length(unique(tmp$clone_id))
  PercUniqueVHHs <- round(VHHs/Seqs, 2)
  
  dat <- rbind(dat, cbind(f, Seqs, SeqAligns, VHHs, ClonalFamilies, PercUniqueVHHs))
}

























## --- ##

source("~/VHH/lib/utilities.R")

dat <- read.csv("~/Downloads/Schisto_run_counts14Aug24.csv", header = T, stringsAsFactors = F)
dat$X <- NULL




options(scipen=999)



dat <- results

dat$database <- ifelse(grepl("VecF", dat$file), "Vectored", "Source")

dat$file <- gsub("_min_threshold_clone-pass.csv", "", dat$file)

tmp <- dat

ggplot(tmp, aes(x=sample_size, y=unique_count, group=file, color=file)) +
  #geom_line(aes(linetype = file), size = 1.2) +
  #geom_line(size = 1.2) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("Schistosome Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Unique Count") 






tmp <- dat

# This gives us the proportion of duplicates
tmp$dup_prop <- 100 - round(tmp$unique_count/tmp$sample_size, 4)*100


ggplot(tmp, aes(x=sample_size, y=dup_prop, group=file, color=file)) +
  #geom_line(aes(linetype = file), size = 1.2) +
  #geom_line(size = 1.2) +
  ylim(c(0, 100)) + xlim(c(0, 500000)) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("Schistosome Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Sequencing Saturation") 














# Could you please create just eight saturation plots of trimmed amino acid 
# sequences - four plots with these source dataset pools: 1-4, 5-8, 9-12, 13-16, 
# and four with these enriched dataset pools: 17, 18, 19-21, 22. 
# This will make identifying each dataset color simpler, 
# and also provide saturation curves containing much larger numbers, 
# some hopefully approaching saturation!

selection <- c("1_L1_sh2_S1", "2_L1_lh2_S2", "3_L1_sh3_S3", "4_L1_lh3_S4")

tmp <- dat[dat$file %in% selection,]

# This gives us the proportion of duplicates
tmp$dup_prop <- 100 - round(tmp$unique_count/tmp$sample_size, 4)*100


ggplot(tmp, aes(x=sample_size, y=dup_prop, group=file, color=file)) +
  #geom_line(aes(linetype = file), size = 1.2) +
  #geom_line(size = 1.2) +
  ylim(c(0, 100)) + xlim(c(0, 500000)) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("Schistosome Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Sequencing Saturation") 



# Could you please create just eight saturation plots of trimmed amino acid 
# sequences - four plots with these source dataset pools: 1-4, 5-8, 9-12, 13-16, 
# and four with these enriched dataset pools: 17, 18, 19-21, 22. 
# This will make identifying each dataset color simpler, 
# and also provide saturation curves containing much larger numbers, 
# some hopefully approaching saturation!


selection <- c("5_L2_sh2_S5",
               "6_L2_lh2_S6",
               "7_L2_sh3_S7",
               "8_L2_lh3_S8")
               
tmp <- dat[dat$file %in% selection,]

# This gives us the proportion of duplicates
tmp$dup_prop <- 100 - round(tmp$unique_count/tmp$sample_size, 4)*100


ggplot(tmp, aes(x=sample_size, y=dup_prop, group=file, color=file)) +
  #geom_line(aes(linetype = file), size = 1.2) +
  #geom_line(size = 1.2) +
  ylim(c(0, 100)) + xlim(c(0, 500000)) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("Schistosome Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Sequencing Saturation") 





# Could you please create just eight saturation plots of trimmed amino acid 
# sequences - four plots with these source dataset pools: 1-4, 5-8, 9-12, 13-16, 
# and four with these enriched dataset pools: 17, 18, 19-21, 22. 
# This will make identifying each dataset color simpler, 
# and also provide saturation curves containing much larger numbers, 
# some hopefully approaching saturation!


selection <- c("9_L3_sh2_S9",
               "10_L3_lh2_S10",
               "11_L3_sh3_S11",
               "12_L3_lh3_S12")

tmp <- dat[dat$file %in% selection,]

# This gives us the proportion of duplicates
tmp$dup_prop <- 100 - round(tmp$unique_count/tmp$sample_size, 4)*100


ggplot(tmp, aes(x=sample_size, y=dup_prop, group=file, color=file)) +
  #geom_line(aes(linetype = file), size = 1.2) +
  #geom_line(size = 1.2) +
  ylim(c(0, 100)) + xlim(c(0, 500000)) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("Schistosome Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Sequencing Saturation") 

















# Could you please create just eight saturation plots of trimmed amino acid 
# sequences - four plots with these source dataset pools: 1-4, 5-8, 9-12, 13-16, 
# and four with these enriched dataset pools: 17, 18, 19-21, 22. 
# This will make identifying each dataset color simpler, 
# and also provide saturation curves containing much larger numbers, 
# some hopefully approaching saturation!


selection <- c("9_L3_sh2_S9",
               "10_L3_lh2_S10",
               "11_L3_sh3_S11",
               "12_L3_lh3_S12")

tmp <- dat[dat$file %in% selection,]

# This gives us the proportion of duplicates
tmp$dup_prop <- 100 - round(tmp$unique_count/tmp$sample_size, 4)*100


ggplot(tmp, aes(x=sample_size, y=dup_prop, group=file, color=file)) +
  #geom_line(aes(linetype = file), size = 1.2) +
  #geom_line(size = 1.2) +
  ylim(c(0, 100)) + xlim(c(0, 500000)) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("Schistosome Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Sequencing Saturation") 



# Could you please create just eight saturation plots of trimmed amino acid 
# sequences - four plots with these source dataset pools: 1-4, 5-8, 9-12, 13-16, 
# and four with these enriched dataset pools: 17, 18, 19-21, 22. 
# This will make identifying each dataset color simpler, 
# and also provide saturation curves containing much larger numbers, 
# some hopefully approaching saturation!


selection <- c("13_VecF1_R2_S13",
               "14_VecF1_R2_S14",
               "15_VecF3_R2_S15",
               "16_VecF4_R2_S16")

tmp <- dat[dat$file %in% selection,]

# This gives us the proportion of duplicates
tmp$dup_prop <- 100 - round(tmp$unique_count/tmp$sample_size, 4)*100


ggplot(tmp, aes(x=sample_size, y=dup_prop, group=file, color=file)) +
  #geom_line(aes(linetype = file), size = 1.2) +
  #geom_line(size = 1.2) +
  ylim(c(0, 100)) + xlim(c(0, 500000)) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("Schistosome Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Sequencing Saturation") 











# Could you please create just eight saturation plots of trimmed amino acid 
# sequences - four plots with these source dataset pools: 1-4, 5-8, 9-12, 13-16, 
# and four with these enriched dataset pools: 17, 18, 19-21, 22. 
# This will make identifying each dataset color simpler, 
# and also provide saturation curves containing much larger numbers, 
# some hopefully approaching saturation!


selection <- c("17_VecF3_R3_S17",
               "18_VecF3_R4_S18",
               "19_VecF1_R3_S19",
               "20_VecF2_R3_S20",
               "21_VecF4_R3_S21",
               "22_VecF4_R4_S22")

tmp <- dat[dat$file %in% selection,]

# This gives us the proportion of duplicates
tmp$dup_prop <- 100 - round(tmp$unique_count/tmp$sample_size, 4)*100


ggplot(tmp, aes(x=sample_size, y=dup_prop, group=file, color=file)) +
  #geom_line(aes(linetype = file), size = 1.2) +
  #geom_line(size = 1.2) +
  ylim(c(0, 100)) + xlim(c(0, 500000)) +
  geom_line(aes(linetype=database), size = 1.2) + 
  theme_minimal() + theme(legend.position="right") + 
  ggtitle("Schistosome Experiment Saturation Curves") +
  xlab("Samples (n)") + ylab("Sequencing Saturation") 
