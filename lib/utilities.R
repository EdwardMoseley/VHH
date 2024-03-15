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