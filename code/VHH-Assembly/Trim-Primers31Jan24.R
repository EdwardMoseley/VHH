# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files("~/VHH/data")[grepl("quality-pass", list.files("~/VHH/data"))]
## Pull their prefixes
assembled_files <- gsub("_quality-pass.fastq", "", files)
## Pull the ones we need
file_set <- assembled_files
#unassembled_files <- assembled_files[which(!(assembled_files %in% new_files))]
#file_set <- c("13_VecF1_R2_S13", "14_VecF2_R2_S14")
# 
#file_set <- setdiff(assembled_files, file_set)
# Undetermined_S0_quality-pass.fastq

## ###
## First is VPRIMER
## ###

## Iterate over each file
for (f in file_set){
  # Some elevator music
  print(f)

  # So we can see what we're doing
  print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
        paste("/data/", f, "_quality-pass.fastq", sep = '' ), 
        " -p /data/AlpVHH_R2_primers_05Feb24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
        paste(f, "_MPV", sep = ''), " --log ",
        paste("/data/", f, "_MPV_TRIM.log", sep = ''), sep = ''))
  
  # Tell the system
  system(
    print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
                paste("/data/", f, "_quality-pass.fastq", sep = '' ), 
                " -p /data/AlpVHH_R2_primers_05Feb24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
                paste(f, "_MPV", sep = ''), " --log ",
                paste("/data/", f, "_MPV_TRIM.log", sep = ''), sep = ''))
  )
}

## ###
## Second is CPRIMER
## ###


# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files("~/VHH/data")[grepl("_MPV_primers-pass.fastq", list.files("~/VHH/data"))]
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
        paste("/data/", f, "_MPV_primers-pass.fastq", sep = '' ), 
        " -p /data/AlpVHH_R1_primers_05Feb24.fasta --start 0 --mode cut --revpr --pf CPRIMER --outname ",
        paste(f, "_MPC", sep = ''), " --log ",
        paste("/data/", f, "_MPV_TRIM.log", sep = ''), sep = '')
        )
  
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
              paste("/data/", f, "_MPV_primers-pass.fastq", sep = '' ), 
              " -p /data/AlpVHH_R1_primers_05Feb24.fasta --start 0 --mode cut --revpr --pf CPRIMER --outname ",
              paste(f, "_MPC", sep = ''), " --log ",
              paste("/data/", f, "_MPC_TRIM.log", sep = ''), sep = '')
  )
}




