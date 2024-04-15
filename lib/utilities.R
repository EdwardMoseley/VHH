## Imports
suppressMessages(library(ggplot2))
suppressMessages(library(msa))
suppressMessages(library(ggseqlogo))
suppressMessages(library(gridExtra))

## Translate DNA
translate_DNA <- function(dna_data){
  clean_dna <- gsub("[^ACTG]", "", dna_data)
  DNA_strings <- Biostrings::DNAStringSet(clean_dna)
  translations <- as.character(Biostrings::translate(DNA_strings))
  return(translations)
}

## Trim CDR3s
trim_CDR3 <- function(seq_vec){
  return(substr(seq_vec, 2, nchar(seq_vec)-1))
}

## Top Groups
top_groups <- function(results, family_name){
  return(table(results[[family_name]])[order(table(results[[family_name]]), decreasing = T)])
}

## Drop sequences with premature stop codons
unproductive <- function(results){
  return(results[which(!grepl("\\*", results$VHH, fixed = F)),])
}

## Calculate the consensus sequence
strict_consensus <- function(char_vec) {
  # Check if all character vectors have the same length
  if (!all(sapply(char_vec, function(x) nchar(x) == nchar(char_vec[1])))) {
    stop("All character vectors must have the same length.")
  }
  # Initialize an empty character vector to store the result
  result_sequence <- character(nchar(char_vec[1]))
  # Iterate over each position in the character vectors
  for (i in seq_len(nchar(char_vec[1]))) {
    # Get the characters at the current index for all vectors
    characters_at_position <- sapply(char_vec, function(x) substr(x, i, i))
    # Find the most common character at the current index
    most_common_char <- names(which.max(table(characters_at_position)))
    # Store it
    result_sequence[i] <- most_common_char
  }
  
  # Return
  return(paste(result_sequence, collapse = ""))
}

# Function to calculate Shannon entropy at the string level
string_entropy <- function(x) {
  # Tabulate to count each unique string
  string_counts <- table(x)
  # Convert to probabilities
  probabilities <- string_counts/length(x)
  # Shannon entropy calculation
  entropy <- -sum(probabilities * log2(probabilities))
  # Return
  return(entropy)
}

## "Unrolling" duplicate values
unroll_dupes <- function(dat, unroll_column){
  # Create to store results
  result_vec <- vector()
  # Iterate over the duplicate count column
  for (i in 1:length(dat$duplicate_count)){
    result_vec <- c(result_vec, 
                    rep(dat[[unroll_column]],
                        each = dat$duplicate_count[i]))
  }
  return(result_vec)
}


## This is more effective
unroll_dupes_fast <- function(dat, unroll_column) {
  # Repeat each element of the unroll column based on duplicate count
  result_vec <- rep(dat[[unroll_column]], times = dat$duplicate_count)
  return(result_vec)
}

## Calculate levenshtein distance
calculate_levenshtein <- function(query_seq, strings) {
  distances <- adist(query_seq, strings)
  return(distances)
}

## Write to fasta for clustering in MMseqs2
write_fasta <- function(df, seq_col, id_col, file_out){
  seqinr::write.fasta(sequences = as.list(df[[seq_col]]), names = df[[id_col]], file.out = file_out)
}

## Cleaning colons out of the sequence_ids from immcantation
clean_seq_ids <- function(seq_col){
  seq_col <- gsub(":", '', seq_col)
  return(seq_col)
}

## mmseqs merging data to clusters
merge_db_to_centroids <- function(db, tmp){
  dat <- merge(db, tmp, by = "sequence_id")
  # Pull out data and prepare to change column names for ease
  tmp <- dat[,c("sequence_id", "centroid")]
  # Update column names of tmp
  colnames(tmp) <- c("tmp_sequence_id", "sequence_id")
  # Join for Centroid sequences
  tmp <- merge(tmp, dat[,c("sequence_id", "VHH")], by = "sequence_id")
  # Correct column names
  colnames(tmp) <- c("centroid_sequence_id", "sequence_id", "centroid")
  # Remove "centroid" from dat
  dat$centroid <- NULL
  # Final merge
  dat <- merge(dat, tmp, by = "sequence_id")
  return(dat)
}


# Create VHH and CDR3 statistics at the family level
family_statistics <- function(dat){
  # Hold results in results
  results <- data.frame()
  # Iterate over clone ids
  for (nm in unique(dat$clone_id)){
    tmp_frame <- dat[dat$clone_id == nm,]
    
    # Unroll
    unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
    unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
    
    # Pull statistics
    median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
    median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
    consensus_CDR3 <- strict_consensus(unrolled_CDR3s)
    
    # Distance from median
    tmp_frame$family_median_VHH_homology <- 100 - (as.vector(adist(median_VHH, tmp_frame$VHH))/nchar(median_VHH))*100
    tmp_frame$family_median_CDR3_homology <- 100 - (as.vector(adist(median_CDR3, tmp_frame$junction_aa))/nchar(median_CDR3))*100
    tmp_frame$family_consensus_CDR3_homology <- 100 - (as.vector(adist(consensus_CDR3, tmp_frame$junction_aa))/nchar(consensus_CDR3))*100
    
    results <- rbind(results, tmp_frame)
  }
  # Ordering
  results <- results[order(results$family_median_VHH_homology, 
                           results$family_median_CDR3_homology, 
                           results$duplicate_count, 
                           decreasing = T),]
  
  # Order statistic
  results$order_index <- 1:nrow(results)
  return(results)
}


# Create VHH and CDR3 statistics at the subfamily level
subfamily_statistics <- function(dat){
  results <- data.frame()
  for (clone_id in unique(dat$clone_id)){
    # Subset clone id
    tmp <- dat[dat$clone_id == as.character(clone_id),]
    # For each centroid id
    for (nm in unique(tmp$centroid_id)){
      #cat("Centroid ID:", nm, '\n')
      tmp_frame <- tmp[tmp$centroid_id == nm,]
      # Unroll strings
      unrolled_VHHs <- unroll_dupes_fast(tmp_frame, "VHH")
      unrolled_CDR3s <- unroll_dupes_fast(tmp_frame, "junction_aa")
      # Pull statistics
      median_VHH <- names(table(unrolled_VHHs)[order(table(unrolled_VHHs), decreasing = T)][1])
      median_CDR3 <- names(table(unrolled_CDR3s)[order(table(unrolled_CDR3s), decreasing = T)][1])
      consensus_CDR3 <- strict_consensus(unrolled_CDR3s)
      # Distance from median
      tmp_frame$subfamily_median_VHH_homology <- 100 - (as.vector(adist(median_VHH, tmp_frame$VHH))/nchar(median_VHH))*100
      # Subfamily Median
      tmp_frame$subfamily_median_CDR3_homology <- 100 - (as.vector(adist(median_CDR3, tmp_frame$junction_aa))/nchar(median_CDR3))*100
      # Subfamily Consensus
      tmp_frame$subfamily_consensus_CDR3_homology <- 100 - (as.vector(adist(consensus_CDR3, tmp_frame$junction_aa))/nchar(consensus_CDR3))*100
      # Stack frames
      results <- rbind(results, tmp_frame)
    }
  }
  # Ordering
  results <- results[order(results$family_median_VHH_homology, 
                           results$subfamily_median_VHH_homology, 
                           results$duplicate_count, 
                           decreasing = T),]
  
  # Order statistic
  results$order_index <- 1:nrow(results)
  
  return(results)
}

# Function to find nearest neighbor and calculate edit distance
find_nearest_neighbor <- function(x, df) {
  # Calculate distance from median_CDR3
  distances <- adist(x, df$median_CDR3)
  # Assign Inf to diagnal to ignore self-matches
  diag(distances) <- Inf
  # Nearest is the minimum
  nearest_index <- which.min(distances)
  # Subset the minimum
  nearest_neighbor <- df$median_CDR3[nearest_index]
  
  # If nearest neighbor is same as original value, 
  if (nearest_neighbor == x) {
    # Assign Inf
    distances[nearest_index] <- Inf
    # Find next minimum
    nearest_index <- which.min(distances)
    # pull it
    nearest_neighbor <- df$median_CDR3[nearest_index]
  }
  
  # Calculate edit distance
  distance <- distances[nearest_index]
  
  # Return in list format
  return(list(nearest_neighbor = nearest_neighbor, 
              distance = distance))
}