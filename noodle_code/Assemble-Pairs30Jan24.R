# List the files we need to assemble pairs for
files <- list.files("~/VHH/data")[grepl("001.fastq", list.files("~/VHH/data"))]
# Drop everything after the L001.*
new_files <- gsub("L001.*", "", files)
# We only need part of that
new_files <- unique(substr(new_files, 1, nchar(new_files)-1))

# Drop one we've already done
new_files <- new_files[new_files != "13_VecF1_R2_S13"]

## [22] "17_VecF3_R3_S17_L001_R2_001.fastq"    "18_VecF3_R4_S18_L001_R1_001.fastq"    "18_VecF3_R4_S18_L001_R2_001.fastq"   
## [25] "19_VecF1_R3_S19_L001_R1_001.fastq"    "19_VecF1_R3_S19_L001_R2_001.fastq"    "2_L1_lh2_S2_L001_R1_001.fastq"  

# Reads one and two
R1 <- "_L001_R1_001.fastq"
R2 <- "_L001_R2_001.fastq"

## Iterate over each file
for (f in new_files[1]){
  # Some elevator music
  print(f)
  print(paste(f, "_L001_R1_001.fastq", sep = ''))# %in% files)
  print(paste(f, "_L001_R2_001.fastq", sep = ''))
  # So we can see what we're doing
  print(    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 AssemblePairs.py align -1 ", 
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


docker run -v /Users/emoseley/VHH/data:/data  immcantation/suite:4.4.0  \
AssemblePairs.py align -1 /data/13_VecF1_R2_S13_L001_R2_001.fastq -2 /data/13_VecF1_R2_S13_L001_R1_001.fastq \
--coord illumina --rc tail --maxerror 0.2 \
--outname 13_VecF1


