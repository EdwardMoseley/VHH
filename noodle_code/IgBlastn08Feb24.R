### IMGT

# List the files 
## Pull the directory
directory <- "~/ncbi-igblast-1.17.1"
## Pull files we've already assembled
files <- list.files(directory)[grepl("_atleast-1", list.files(directory))]
## Pull their prefixes
assembled_files <- unique(gsub("_atleast-1.*", "", files))
## Pull the ones we need
file_set <- assembled_files

## Set the working directory
setwd("~/ncbi-igblast-1.17.1")

## Iterate over each file
for (f in file_set[2:length(file_set)]){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(
    paste("bin/igblastn -germline_db_V database/alpaca_vhh_clean/alpaca_V -germline_db_D ",
          "database/alpaca_vhh_clean/alpaca_D -germline_db_J database/alpaca_vhh_clean/alpaca_J ",
          "-domain_system imgt -ig_seqtype Ig -organism alpaca -auxiliary_data optional_file/alpaca_gl.aux ",
          "-outfmt '7 std qseq sseq btop' ",
          paste("-query ", directory,'/', f, "_atleast-1.fasta ", sep = ''),
          paste("-out ", directory, '/', f, ".fmt7", sep = ''), sep = '')
  )
  
  # Tell the system
  system(
    paste("bin/igblastn -germline_db_V database/alpaca_vhh_clean/alpaca_V -germline_db_D ",
          "database/alpaca_vhh_clean/alpaca_D -germline_db_J database/alpaca_vhh_clean/alpaca_J ",
          "-domain_system imgt -ig_seqtype Ig -organism alpaca -auxiliary_data optional_file/alpaca_gl.aux ",
          "-outfmt '7 std qseq sseq btop' ",
          paste("-query ", directory,'/', f, "_atleast-1.fasta ", sep = ''),
          paste("-out ", directory, '/', f, ".fmt7", sep = ''), sep = '')
  )
}



#```
## Copy back to data
#cp S1-M02-TRIM_atleast-1.fmt7 ~/VHH/data
#cp S1-M02-TRIM_atleast-1.fasta ~/VHH/data
#```

#### Make Database

#```{bash, eval = F}
#docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 \
#MakeDb.py igblast -i /data/S1-M02-TRIM_atleast-1.fmt7 -s /data/S1-M02-TRIM_atleast-1.fasta \
#-r /data/imgt-data/imgt_gapped_alpaca_IGHV.fasta \
#/data/imgt-data/imgt_gapped_alpaca_IGHD.fasta  \
#/data/imgt-data/imgt_gapped_alpaca_IGHJ.fasta \
#--outdir /data/ --failed --format airr --extended
#```
