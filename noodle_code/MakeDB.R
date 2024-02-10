## ### 
## Make Database
## ###



# Copy back to data from ncbi
# cp *.fmt7 ~/VHH/data
# cp *.fasta ~/VHH/data

# List the files we need to assemble pairs for
## Pull files we've already assembled
files <- list.files("~/VHH/data")[grepl(".fmt7", list.files("~/VHH/data"))]
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
          paste("/data/", f, ".fmt7", sep = ''),
          paste(" -s ", sep = ''),
          paste("/data/", f, "_atleast-1.fasta", sep = ''), 
          paste(" -r /data/imgt-data/imgt_gapped_alpaca_IGHV.fasta ", 
                "/data/imgt-data/imgt_gapped_alpaca_IGHD.fasta ",
                "/data/imgt-data/imgt_gapped_alpaca_IGHJ.fasta ",
                " --outdir /data/ --failed --format airr --extended ", sep = ''), sep = '')
    
  )
  
  # Tell the system
  system(
    paste("docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 ",
          paste("MakeDb.py igblast -i ", sep = ''),
          paste("/data/", f, ".fmt7", sep = ''),
          paste(" -s ", sep = ''),
          paste("/data/", f, "_atleast-1.fasta", sep = ''), 
          paste(" -r /data/imgt-data/imgt_gapped_alpaca_IGHV.fasta ", 
                "/data/imgt-data/imgt_gapped_alpaca_IGHD.fasta ",
                "/data/imgt-data/imgt_gapped_alpaca_IGHJ.fasta ",
                " --outdir /data/ --failed --format airr --extended ", sep = ''), sep = '')
  )
}

