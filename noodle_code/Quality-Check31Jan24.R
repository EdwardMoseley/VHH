# Pull files we've already assembled
assembled_files <- list.files("~/VHH/data")[grepl("assemble-pass", list.files("~/VHH/data"))]
# Pull their prefixes
assembled_files <- gsub("_assemble-pass.fastq", "", assembled_files)
# Pull the ones we need
# unassembled_files <- assembled_files[which(!(assembled_files %in% new_files))]

suffix <- "_assemble-pass.fastq"

## Iterate over each file
for (f in assembled_files){
  # Some elevator music
  print(f)
  
  # So we can see what we're doing
  print(paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 FilterSeq.py quality -s ",
              paste("/data/", f, suffix, sep = '' ), 
              " -q 20 --outname ",
              f, sep = ''))
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 FilterSeq.py quality -s ",
          paste("/data/", f, suffix, sep = '' ), 
          " -q 20 --outname ",
          f, sep = '')
  )
}


