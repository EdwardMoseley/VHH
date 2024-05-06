## Make Database

### IMGT

# List the files 
## Pull the directory
ncbi_directory <- "~/ncbi-igblast-1.17.1"
## Pull files we've already assembled
files <- list.files(tmp_dir)[grepl("_atleast-1", list.files(tmp_dir))]
## Pull their prefixes
assembled_files <- unique(gsub("_atleast-1.*", "", files))
## Pull the ones we need
file_set <- assembled_files

# First, we'll read the fastq from ~/VHH/data/tmp
# Then we'll move it to ncbi_directory as a fasta

tmp <- readFastq(paste(tmp_dir, paste(assembled_files, "_atleast-1.fastq", sep = ''), sep = '/'))

# We need it in the ncbi directory AND the data directory (for making database)
writeFasta(tmp, file = paste(ncbi_directory, paste(assembled_files, "_atleast-1.fasta", sep = ''), sep = '/'))
writeFasta(tmp, file = paste(tmp_dir, paste(assembled_files, "_atleast-1.fasta", sep = ''), sep = '/'))

# Cleanup
rm(tmp)

## Set the working directory
setwd(ncbi_directory)

## Iterate over each file
for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(
    paste("bin/igblastn -germline_db_V database/alpaca_vhh_clean/alpaca_V -germline_db_D ",
          "database/alpaca_vhh_clean/alpaca_D -germline_db_J database/alpaca_vhh_clean/alpaca_J ",
          "-domain_system imgt -ig_seqtype Ig -organism alpaca -auxiliary_data optional_file/alpaca_gl.aux ",
          "-outfmt '7 std qseq sseq btop' ",
          paste("-query ", ncbi_directory,'/', f, "_atleast-1.fasta ", sep = ''),
          paste("-out ", tmp_dir, '/', f, ".fmt7", sep = ''), sep = '')
  )
  
  # Tell the system
  system(
    paste("bin/igblastn -germline_db_V database/alpaca_vhh_clean/alpaca_V -germline_db_D ",
          "database/alpaca_vhh_clean/alpaca_D -germline_db_J database/alpaca_vhh_clean/alpaca_J ",
          "-domain_system imgt -ig_seqtype Ig -organism alpaca -auxiliary_data optional_file/alpaca_gl.aux ",
          "-outfmt '7 std qseq sseq btop' ",
          paste("-query ", ncbi_directory,'/', f, "_atleast-1.fasta ", sep = ''),
          paste("-out ", tmp_dir, '/', f, ".fmt7", sep = ''), sep = '')
  )
}

## Make Database

# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files(tmp_dir)[grepl(".fmt7", list.files(tmp_dir))]
## Pull their prefixes
assembled_files <- gsub(".fmt7", "", files)
## Pull the ones we need
file_set <- assembled_files

## Iterate over each file
for (f in file_set){
  # Some elevator musi
  print(f)
  
  # So we can see what we're doing
  print(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 ",
          paste("MakeDb.py igblast -i ", sep = ''),
          paste("/data/tmp/", f, ".fmt7", sep = ''),
          paste(" -s ", sep = ''),
          paste("/data/tmp/", f, "_atleast-1.fasta", sep = ''), 
          paste(" -r /data/imgt-data/imgt_gapped_alpaca_IGHV.fasta ", 
                "/data/imgt-data/imgt_gapped_alpaca_IGHD.fasta ",
                "/data/imgt-data/imgt_gapped_alpaca_IGHJ.fasta ",
                " --outdir /data/Assembled/ --failed --format airr --extended ", sep = ''), sep = '')
    
  )
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 ",
          paste("MakeDb.py igblast -i ", sep = ''),
          paste("/data/tmp/", f, ".fmt7", sep = ''),
          paste(" -s ", sep = ''),
          paste("/data/tmp/", f, "_atleast-1.fasta", sep = ''), 
          paste(" -r /data/imgt-data/imgt_gapped_alpaca_IGHV.fasta ", 
                "/data/imgt-data/imgt_gapped_alpaca_IGHD.fasta ",
                "/data/imgt-data/imgt_gapped_alpaca_IGHJ.fasta ",
                " --outdir /data/Assembled/ --failed --format airr --extended ", sep = ''), sep = '')
  )
}


