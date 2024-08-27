## ###
## VHH Assembly Script
## ###

if(!require(ShortRead)){
  install.packages("ShortRead")
}

## Libraries
library(ShortRead)

## Assemble Pairs

# List the files we need to assemble pairs for
root_dir <- "~/VHH/data/Assembly"
tmp_dir <- "~/VHH/data/tmp"
final_dir <- "~/VHH/data/Assembled"
files <- list.files(root_dir)


ncbi_directory <- "~/ncbi-igblast-1.17.1"
## Set the working directory
setwd(ncbi_directory)


# files <- c("5_S3_L001_R1_001.fastq")

# Unique because there's two for each
files <- unique(gsub("*_L001_R[1|2]_001.fastq", "", files))
# Drop extensions
files <- files[!grepl("\\.", files)]
files <- files[!(files %in% c("master_set","C1_S1","subsets_db-pass","subsets_min_threshold","Undetermined_S0"))]

# Find prefixes of files again to ensure we don't run through the data twice
completed_files <- gsub(".fmt7", "", list.files(final_dir))

# We want to do 9_S18, so:
#completed_files <- completed_files[!(completed_files %in% "")]#c("9_S18", "10_S19", "11_S20", "12_S21"))] #

# Reads one and two
R1 <- "_L001_R1_001.fastq"
R2 <- "_L001_R2_001.fastq"


for (f in files){
  if (!(f %in% completed_files)){

## Iterate over each file
#for (f in files){
  # Some elevator music
  print(f)
  print(paste(f, "_L001_R1_001.fastq", sep = ''))# %in% files)
  print(paste(f, "_L001_R2_001.fastq", sep = ''))
  # So we can see what we're doing
  print(    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 AssemblePairs.py align -1 ", 
                  paste("/data/Assembly/", f, R2, sep = ''),  
                  " -2 ", paste("/data/Assembly/", f, R1, sep = ''),
                  " --coord illumina --rc tail --maxerror 0.2 --outname ", 
                  paste(f, sep = ''), sep = ''))
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 AssemblePairs.py align -1 ", 
          paste("/data/Assembly/", f, R2, sep = ''),  
          " -2 ", paste("/data/Assembly/", f, R1, sep = ''),
          " --coord illumina --rc tail --maxerror 0.2 --outname ", 
          paste(f, sep = ''), sep = '')
  )
#}



## Quality Check



# Pull files we've already assembled
#assembled_files <- list.files(tmp_dir)
# Pull their prefixes
#assembled_files <- gsub("_assemble-pass.fastq", "", assembled_files)
# Pull the ones we need
#assembled_files <- assembled_files[!grepl("_quality-pass.fastq", assembled_files)]
# unassembled_files <- assembled_files[which(!(assembled_files %in% new_files))]

suffix <- "_assemble-pass.fastq"

## Iterate over each file
#for (f in assembled_files){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 FilterSeq.py quality -s ",
              paste("/data/Assembly/", f, suffix, sep = '' ), 
              " -q 20 --outname ",
              f, sep = ''))
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 FilterSeq.py quality -s ",
          paste("/data/Assembly/", f, suffix, sep = '' ), 
          " -q 20 --outname ",
          f, sep = '')
  )
#}



## Mask Primers


# List the files we need to assemble pairs for
## Pull files we've already assembled
#files <- list.files(tmp_dir)[grepl("quality-pass", list.files(tmp_dir))]
## Pull their prefixes
#assembled_files <- gsub("_quality-pass.fastq", "", files)
## Pull the ones we need
#file_set <- assembled_files
  
suffix <- "_quality-pass.fastq"


## ###
## First is VPRIMER
## ###

## Iterate over each file
#for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
              paste("/data/Assembly/", f, suffix, sep = '' ), 
              " -p /data/AlpVHH_R2_primers_05Feb24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
              paste(f, "_MPV", sep = ''), " --log ",
              paste("/data/Assembly/", f, "_MPV_TRIM.log", sep = ''), sep = ''))
  
  # Tell the system
  system(
    print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
                paste("/data/Assembly/", f, suffix, sep = '' ), 
                " -p /data/AlpVHH_R2_primers_05Feb24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
                paste(f, "_MPV", sep = ''), " --log ",
                paste("/data/Assembly/", f, "_MPV_TRIM.log", sep = ''), sep = ''))
  )
#}

## ###
## Second is CPRIMER
## ###
  
suffix <- "_MPV_primers-pass.fastq"


# List the files we need to assemble pairs for
## Pull files we've already assembled
#files <- list.files(tmp_dir)[grepl("_MPV_primers-pass.fastq", list.files(tmp_dir))]
## Pull their prefixes
#assembled_files <- gsub("_MPV_primers-pass.fastq", "", files)
## Pull the ones we need
#file_set <- assembled_files

## Iterate over each file
#for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
          paste("/data/Assembly/", f, suffix, sep = '' ), 
          " -p /data/AlpVHH_R1_primers_05Feb24.fasta --start 0 --mode cut --revpr --pf CPRIMER --outname ",
          paste(f, "_MPC", sep = ''), " --log ",
          paste("/data/Assembly/", f, "_MPV_TRIM.log", sep = ''), sep = '')
  )
  
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
          paste("/data/Assembly/", f, suffix, sep = '' ), 
          " -p /data/AlpVHH_R1_primers_05Feb24.fasta --start 0 --mode cut --revpr --pf CPRIMER --outname ",
          paste(f, "_MPC", sep = ''), " --log ",
          paste("/data/Assembly/", f, "_MPC_TRIM.log", sep = ''), sep = '')
  )
#}

## ###
## EM QUICK EDIT
## ###
  
##

## Collapse and Deduplicate

suffix <- "_MPC_primers-pass.fastq"


# List the files we need to assemble pairs for
## Pull files we've already assembled
#files <- list.files(tmp_dir)[grepl("_MPC_primers-pass.fastq", list.files(tmp_dir))]
## Pull their prefixes
#assembled_files <- gsub("_MPC_primers-pass.fastq", "", files)
## Pull the ones we need
#file_set <- assembled_files

## Iterate over each file
#for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 CollapseSeq.py -s ",
          paste("/data/Assembly/", f, suffix, sep = '' ), 
          " -n 20 --inner --uf CPRIMER --cf VPRIMER --act set --outname ",
          paste(f, sep = ''), sep = '')
    
  )
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 CollapseSeq.py -s ",
          paste("/data/Assembly/", f, suffix, sep = '' ), 
          " -n 20 --inner --uf CPRIMER --cf VPRIMER --act set --outname ",
          paste(f, sep = ''), sep = '')
  )
#}

# List the files we need to assemble pairs for
## Pull files we've already assembled
#files <- list.files(tmp_dir)[grepl("_collapse-unique.fastq", list.files(tmp_dir))]
## Pull their prefixes
#assembled_files <- gsub("_collapse-unique.fastq", "", files)
## Pull the ones we need
#file_set <- assembled_files
  
suffix <- "_collapse-unique.fastq"

## Iterate over each file
#for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 SplitSeq.py group -s ",
          paste("/data/Assembly/", f, suffix, sep = '' ), 
          " -f DUPCOUNT --num 1 --outname ",
          paste(f, sep = ''), sep = '')
    
  )
  
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 SplitSeq.py group -s ",
          paste("/data/Assembly/", f, suffix, sep = '' ), 
          " -f DUPCOUNT --num 1 --outname ",
          paste(f, sep = ''), sep = '')
  )
  
  
  
  
  
## ###
## Make Database
## ###

# List the files 
## Pull the directory
#ncbi_directory <- "~/ncbi-igblast-1.17.1"
## Pull files we've already assembled
#files <- list.files(tmp_dir)[grepl("_atleast-1", list.files(tmp_dir))]
## Pull their prefixes
#assembled_files <- unique(gsub("_atleast-1.*", "", files))
## Pull the ones we need
#file_set <- assembled_files


# First, we'll read the fastq from ~/VHH/data/tmp
# Then we'll move it to ncbi_directory as a fasta

tmp <- readFastq(paste(root_dir, paste(f, "_atleast-1.fastq", sep = ''), sep = '/'))

# We need it in the ncbi directory AND the data directory (for making database)
writeFasta(tmp, file = paste(ncbi_directory, paste(f, "_atleast-1.fasta", sep = ''), sep = '/'))
writeFasta(tmp, file = paste(final_dir, paste(f, "_atleast-1.fasta", sep = ''), sep = '/'))

# Cleanup
rm(tmp)

# Some elevator music
print(f)

# So we can see what we're doing
print(
  paste("bin/igblastn -germline_db_V database/alpaca_vhh_clean/alpaca_V -germline_db_D ",
        "database/alpaca_vhh_clean/alpaca_D -germline_db_J database/alpaca_vhh_clean/alpaca_J ",
        "-domain_system imgt -ig_seqtype Ig -organism alpaca -auxiliary_data optional_file/alpaca_gl.aux ",
        "-outfmt '7 std qseq sseq btop' ",
        paste("-query ", ncbi_directory,'/', f, "_atleast-1.fasta ", sep = ''),
        paste("-out ", final_dir, '/', f, ".fmt7", sep = ''), sep = '')
)

# Tell the system
system(
  paste("bin/igblastn -germline_db_V database/alpaca_vhh_clean/alpaca_V -germline_db_D ",
        "database/alpaca_vhh_clean/alpaca_D -germline_db_J database/alpaca_vhh_clean/alpaca_J ",
        "-domain_system imgt -ig_seqtype Ig -organism alpaca -auxiliary_data optional_file/alpaca_gl.aux ",
        "-outfmt '7 std qseq sseq btop' ",
        paste("-query ", ncbi_directory,'/', f, "_atleast-1.fasta ", sep = ''),
        paste("-out ", final_dir, '/', f, ".fmt7", sep = ''), sep = '')
)
  
# Some elevator musi
print(f)

# So we can see what we're doing
print(
  paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 ",
        paste("MakeDb.py igblast -i ", sep = ''),
        paste("/data/Assembled", f, ".fmt7", sep = ''),
        paste(" -s ", sep = ''),
        paste("/data/Assembled", f, "_atleast-1.fasta", sep = ''), 
        paste(" -r /data/imgt-data/imgt_gapped_alpaca_IGHV.fasta ", 
              "/data/imgt-data/imgt_gapped_alpaca_IGHD.fasta ",
              "/data/imgt-data/imgt_gapped_alpaca_IGHJ.fasta ",
              " --outdir /data/Assembled/ --failed --format airr --extended ", sep = ''), sep = '')
  
)

# Tell the system
system(
  paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 ",
        paste("MakeDb.py igblast -i ", sep = ''),
        paste("data/Assembled/", f, ".fmt7", sep = ''),
        paste(" -s ", sep = ''),
        paste("data/Assembled/", f, "_atleast-1.fasta", sep = ''), 
        paste(" -r /data/imgt-data/imgt_gapped_alpaca_IGHV.fasta ", 
              "/data/imgt-data/imgt_gapped_alpaca_IGHD.fasta ",
              "/data/imgt-data/imgt_gapped_alpaca_IGHJ.fasta ",
              " --outdir /data/Assembled/ --failed --format airr --extended ", sep = ''), sep = '')
)
  
#}
  # Find prefixes of files again to ensure we don't run through the data twice
  completed_files <- gsub(".fmt7", "", list.files(final_dir))
  }
}
