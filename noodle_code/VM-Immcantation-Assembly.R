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
for (f in new_files){
# Some elevator music
  print(f)
  print(paste(f, "_L001_R1_001.fastq", sep = ''))# %in% files)
  print(paste(f, "_L001_R2_001.fastq", sep = ''))
# So we can see what we're doing
  print(    paste("docker run -v /Users/emoseley/VHH/data/Assembly:/data immcantation/suite:4.4.0 AssemblePairs.py align -1 ", 
                  paste("/data/", f, R2, sep = ''),  
                  " -2 ", paste("/data/", f, R1, sep = ''),
                  " --coord illumina --rc tail --maxerror 0.2 --outname ", 
                  paste(f, sep = ''), sep = ''))

# Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 AssemblePairs.py align -1 ", 
          paste("/data/", f, R2, sep = ''),  
          " -2 ", paste("/data/", f, R1, sep = ''),
          " --coord illumina --rc tail --maxerror 0.2 --outname ", 
          paste(f, sep = ''), sep = '')
  )
}






!/usr/bin/python3 /root/.local/bin/AssemblePairs.py align -1 5_S3_L001_R1_001.fastq -2 5_S3_L001_R2_001.fastq \
--coord sra --rc tail --outname M1 --log AP.log