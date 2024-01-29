## Imports
suppressMessages(library(ggplot2))
suppressMessages(library(ggseqlogo))
suppressMessages(library(gridExtra))

## Translate DNA
translate_DNA <- function(dna_data){
  clean_dna <- gsub("[^ACTG]", "", dna_data)
  DNA_strings <- Biostrings::DNAStringSet(clean_dna)
  translations <- as.character(Biostrings::translate(DNA_strings))
  return(translations)
}

## ###
## Trim CDR3
## Trims the first and last Amino Acid of the CDR3
## These are usually C and W
## ###

trim_CDR3 <- function(seq_vec){
  return(substr(seq_vec, 2, nchar(seq_vec)-1))
}

## Top Groups
top_groups <- function(results){
  return(table(results$family)[order(table(results$family), decreasing = T)])
}

## Drop sequences with premature stop codons
unproductive <- function(results){
  return(results[which(!grepl("\\*", results$VHH, fixed = F)),])
}

## Plots sequence logo with a bunch of other info, will need to clean this
plot_logo <- function(tmp_frame){
  common_CDR3 <- names(table(tmp_frame$junction_aa)[order(table(tmp_frame$junction_aa), decreasing = T)][1])
  
  cat("Illumina Run: ", tmp_frame$experiment[1], ' \n', sep = '')
  cat("Group ID: ", unique(tmp_frame$family), '\n', sep = '')
  cat("Unique VHHs: ", length(unique(tmp_frame$VHH)), '\n', sep = '')
  cat("Unique CDR3s: ", length(unique(tmp_frame$junction_aa)), '\n', sep = '')
  cat("Most-frequent CDR3: ", common_CDR3, ' (with ',nrow(tmp_frame[tmp_frame$junction_aa == common_CDR3,]), ' CDR3 observations)\n', sep = '')
  print(ggseqlogo(as.character(tmp_frame$junction_aa), seq_type='aa', method = 'p') + 
          theme_classic() + ggtitle(paste(" CDR3 Logo Plot for cluster with the most common CDR3 sequence:\n               ", common_CDR3, sep = '')))
  
}