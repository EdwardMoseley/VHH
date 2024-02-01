# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files("~/VHH/data")[grepl("quality-pass", list.files("~/VHH/data"))]
## Pull their prefixes
assembled_files <- gsub("_quality-pass.fastq", "", files)
## Pull the ones we need
#unassembled_files <- assembled_files[which(!(assembled_files %in% new_files))]
file_set <- c("13_VecF1_R2_S13", "14_VecF2_R2_S14")
# 

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
        " -p /data/AlpVHH_R2_primers_07Jan24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
        paste(f, "_MPV", sep = ''), " --log ",
        paste("/data/", f, "_TRIM.log", sep = ''), sep = ''))
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
          paste("/data/", f, "_quality-pass.fastq", sep = '' ), 
          " -p /data/AlpVHH_R2_primers_07Jan24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
          paste(f, "_MPV", sep = ''), " --log ",
          paste("/data/", f, "_TRIM.log", sep = ''), sep = '')
  )
}

## ###
## Second is CPRIMER
## ###

#docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  \
#MaskPrimers.py score -s /data/13_VecF1_R2_S13-TRIM_primers-pass.fastq -p /data/AlpVHH_R1_primers_07Jan24.fasta \
#--start 0 --mode cut --revpr --pf CPRIMER --outname 13_VecF1_R2_S13-MPC-TRIM --log /data/MPC-13_VecF1_R2_S13-TRIM.log


# 13_VecF1_R2_S13_MPV_primers-pass.fastq

## Iterate over each file
for (f in file_set){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
              paste("/data/", f, "_MPV-primers-pass.fastq", sep = '' ), 
              " -p /data/AlpVHH_R1_primers_07Jan24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
              paste(f, "_MPC", sep = ''), " --log ",
              paste("/data/", f, "_MPC.log", sep = ''), sep = ''))
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0  MaskPrimers.py score -s ",
              paste("/data/", f, "_MPV-primers-pass.fastq", sep = '' ), 
              " -p /data/AlpVHH_R1_primers_07Jan24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
              paste(f, "_MPC", sep = ''), " --log ",
              paste("/data/", f, "_MPC.log", sep = ''), sep = '')
  )
}




