# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files("~/VHH/data")[grepl("_MPC_primers-pass.fastq", list.files("~/VHH/data"))]
## Pull their prefixes
assembled_files <- gsub("_MPC_primers-pass.fastq", "", files)
## Pull the ones we need
file_set <- assembled_files

## Iterate over each file
for (f in file_set[1]){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 CollapseSeq.py -s ",
          paste("/data/", f, "_MPC_primers-pass.fastq", sep = '' ), 
          " -n 20 --inner --uf CPRIMER --cf VPRIMER --act set --outname ",
          paste(f, sep = ''), sep = '')
    
  )
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 CollapseSeq.py -s ",
          paste("/data/", f, "_MPC_primers-pass.fastq", sep = '' ), 
          " -n 20 --inner --uf CPRIMER --cf VPRIMER --act set --outname ",
          paste(f, sep = ''), sep = '')
  )
}

# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files("~/VHH/data")[grepl("_collapse-unique.fastq", list.files("~/VHH/data"))]
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
          paste("/data/", f, "_collapse-unique.fastq", sep = '' ), 
          " -f DUPCOUNT --num 1 --outname ",
          paste(f, sep = ''), sep = '')
    
  )
  
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 SplitSeq.py group -s ",
          paste("/data/", f, "_collapse-unique.fastq", sep = '' ), 
          " -f DUPCOUNT --num 1 --outname ",
          paste(f, sep = ''), sep = '')
  )
}