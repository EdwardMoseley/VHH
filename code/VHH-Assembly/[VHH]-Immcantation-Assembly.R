## ###
## VHH Assembly Script
## ###

## Assemble Pairs

# List the files we need to assemble pairs for
root_dir <- "~/VHH/data/Assembly"
tmp_dir <- "~/VHH/data/tmp"
final_dir <- "~/VHH/data/Assembled"
files <- list.files(root_dir)

# Reads one and two
#R1 <- "_L001_R1_001.fastq"
#R2 <- "_L001_R2_001.fastq"

## Iterate over each file
#for (f in new_files){
  # Some elevator music
#  print(f)
#  print(paste(f, "_L001_R1_001.fastq", sep = ''))# %in% files)
#  print(paste(f, "_L001_R2_001.fastq", sep = ''))
  # So we can see what we're doing
#  print(    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 AssemblePairs.py align -1 ", 
#                  paste("/data/", f, R2, sep = ''),  
#                  " -2 ", paste("/data/", f, R1, sep = ''),
#                  " --coord illumina --rc tail --maxerror 0.2 --outname ", 
#                  paste(f, sep = ''), sep = ''))
  
  # Tell the system
#  system(
#    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 AssemblePairs.py align -1 ", 
#          paste("/data/", f, R2, sep = ''),  
#          " -2 ", paste("/data/", f, R1, sep = ''),
#          " --coord illumina --rc tail --maxerror 0.2 --outname ", 
#          paste(f, sep = ''), sep = '')
#  )
#}



## Quality Check



# Pull files we've already assembled
assembled_files <- list.files(tmp_dir)
# Pull their prefixes
assembled_files <- gsub("_assemble-pass.fastq", "", assembled_files)
# Pull the ones we need
assembled_files <- assembled_files[!grepl("_quality-pass.fastq", assembled_files)]
# unassembled_files <- assembled_files[which(!(assembled_files %in% new_files))]

suffix <- "_assemble-pass.fastq"

## Iterate over each file
for (f in assembled_files){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 FilterSeq.py quality -s ",
              paste("/data/tmp/", f, suffix, sep = '' ), 
              " -q 5 --outname ",
              f, sep = ''))
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 FilterSeq.py quality -s ",
          paste("/data/tmp/", f, suffix, sep = '' ), 
          " -q 5 --outname ",
          f, sep = '')
  )
}



## Mask Primers


# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files(tmp_dir)[grepl("quality-pass", list.files(tmp_dir))]
## Pull their prefixes
assembled_files <- gsub("_quality-pass.fastq", "", files)
## Pull the ones we need
file_set <- assembled_files


## ###
## First is VPRIMER
## ###

## Iterate over each file
for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
              paste("/data/tmp/", f, "_quality-pass.fastq", sep = '' ), 
              " -p /data/AlpVHH_R2_primers_05Feb24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
              paste(f, "_MPV", sep = ''), " --log ",
              paste("/data/tmp/", f, "_MPV_TRIM.log", sep = ''), sep = ''))
  
  # Tell the system
  system(
    print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
                paste("/data/tmp/", f, "_quality-pass.fastq", sep = '' ), 
                " -p /data/AlpVHH_R2_primers_05Feb24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
                paste(f, "_MPV", sep = ''), " --log ",
                paste("/data/tmp/", f, "_MPV_TRIM.log", sep = ''), sep = ''))
  )
}

## ###
## Second is CPRIMER
## ###


# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files(tmp_dir)[grepl("_MPV_primers-pass.fastq", list.files(tmp_dir))]
## Pull their prefixes
assembled_files <- gsub("_MPV_primers-pass.fastq", "", files)
## Pull the ones we need
file_set <- assembled_files

## Iterate over each file
for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
          paste("/data/tmp/", f, "_MPV_primers-pass.fastq", sep = '' ), 
          " -p /data/AlpVHH_R1_primers_05Feb24.fasta --start 0 --mode cut --revpr --pf CPRIMER --outname ",
          paste(f, "_MPC", sep = ''), " --log ",
          paste("/data/tmp/", f, "_MPV_TRIM.log", sep = ''), sep = '')
  )
  
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
          paste("/data/tmp/", f, "_MPV_primers-pass.fastq", sep = '' ), 
          " -p /data/AlpVHH_R1_primers_05Feb24.fasta --start 0 --mode cut --revpr --pf CPRIMER --outname ",
          paste(f, "_MPC", sep = ''), " --log ",
          paste("/data/tmp/", f, "_MPC_TRIM.log", sep = ''), sep = '')
  )
}


## Collapse and Deduplicate




# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files(tmp_dir)[grepl("_MPC_primers-pass.fastq", list.files(tmp_dir))]
## Pull their prefixes
assembled_files <- gsub("_MPC_primers-pass.fastq", "", files)
## Pull the ones we need
file_set <- assembled_files

## Iterate over each file
for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 CollapseSeq.py -s ",
          paste("/data/tmp/", f, "_MPC_primers-pass.fastq", sep = '' ), 
          " -n 20 --inner --uf CPRIMER --cf VPRIMER --act set --outname ",
          paste(f, sep = ''), sep = '')
    
  )
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 CollapseSeq.py -s ",
          paste("/data/tmp/", f, "_MPC_primers-pass.fastq", sep = '' ), 
          " -n 20 --inner --uf CPRIMER --cf VPRIMER --act set --outname ",
          paste(f, sep = ''), sep = '')
  )
}

# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files(tmp_dir)[grepl("_collapse-unique.fastq", list.files(tmp_dir))]
## Pull their prefixes
assembled_files <- gsub("_collapse-unique.fastq", "", files)
## Pull the ones we need
file_set <- assembled_files

## Iterate over each file
for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 SplitSeq.py group -s ",
          paste("/data/tmp/", f, "_collapse-unique.fastq", sep = '' ), 
          " -f DUPCOUNT --num 1 --outname ",
          paste(f, sep = ''), sep = '')
    
  )
  
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 SplitSeq.py group -s ",
          paste("/data/tmp/", f, "_collapse-unique.fastq", sep = '' ), 
          " -f DUPCOUNT --num 1 --outname ",
          paste(f, sep = ''), sep = '')
  )
}

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


