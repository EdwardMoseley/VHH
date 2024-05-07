## ###
## VHH Assembly Script

  #if(!require(ShortRead)){
  #  install.packages("ShortRead")
  #}
  
  ## Libraries
  #library(ShortRead)
  
  ## Assemble Pairs
  
  # List the files we need to assemble pairs for
  root_dir <- "~/VHH/data/Assembly"
  tmp_dir <- "~/VHH/data/tmp"
  final_dir <- "~/VHH/data/Assembled"
  files <- list.files(root_dir)
  
  
  #ncbi_directory <- "~/ncbi-igblast-1.17.1"
  ## Set the working directory
  #setwd(ncbi_directory)
  
  
  # files <- c("5_S3_L001_R1_001.fastq")
  
  # Unique because there's two for each
  files <- unique(gsub("*_L001_R[1|2]_001.fastq", "", files))
  # Drop extensions
  files <- files[!grepl("\\.", files)]
  
  
  
  # Find prefixes of files again to ensure we don't run through the data twice
  completed_files <- gsub(".fmt7", "", list.files(final_dir))
  
  # Reads one and two
  R1 <- "_L001_R1_001.fastq"
  R2 <- "_L001_R2_001.fastq"
  
  
  for (f in rev(files)){
    if (!(f %in% completed_files)){
      
      ## Iterate over each file
      #for (f in files){
      # Some elevator music
      print(f)
      print(paste(f, "_L001_R1_001.fastq", sep = ''))# %in% files)
      print(paste(f, "_L001_R2_001.fastq", sep = ''))
      # So we can see what we're doing
      print(    paste("python3.12 .local/bin/AssemblePairs.py align -1 ", 
                      paste("/data/Assembly/", f, R2, sep = ''),  
                      " -2 ", paste("/data/Assembly/", f, R1, sep = ''),
                      " --coord illumina --rc tail --maxerror 0.2 --outname ", 
                      paste(f, sep = ''), sep = ''))
      
      # Tell the system
      system(
        paste("python3.12 .local/bin/AssemblePairs.py align -1 ", 
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
      print(paste("python3.12 .local/bin/FilterSeq.py quality -s ",
                  paste("/data/Assembly/", f, suffix, sep = '' ), 
                  " -q 20 --outname ",
                  f, sep = ''))
      
      # Tell the system
      system(
        paste("python3.12 .local/bin/FilterSeq.py quality -s ",
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
      print(paste("python3.12 .local/bin/MaskPrimers.py score -s ",
                  paste("/data/Assembly/", f, suffix, sep = '' ), 
                  " -p /data/AlpVHH_R2_primers_05Feb24.fasta --start 0 --mode cut --pf VPRIMER --outname ",
                  paste(f, "_MPV", sep = ''), " --log ",
                  paste("/data/Assembly/", f, "_MPV_TRIM.log", sep = ''), sep = ''))
      
      # Tell the system
      system(
        print(paste("python3.12 .local/bin/MaskPrimers.py score -s ",
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
        paste("python3.12 .local/bin/MaskPrimers.py score -s ",
              paste("/data/Assembly/", f, suffix, sep = '' ), 
              " -p /data/AlpVHH_R1_primers_05Feb24.fasta --start 0 --mode cut --revpr --pf CPRIMER --outname ",
              paste(f, "_MPC", sep = ''), " --log ",
              paste("/data/Assembly/", f, "_MPV_TRIM.log", sep = ''), sep = '')
      )
      
      
      # Tell the system
      system(
        paste("python3.12 .local/bin/MaskPrimers.py score -s ",
              paste("/data/Assembly/", f, suffix, sep = '' ), 
              " -p /data/AlpVHH_R1_primers_05Feb24.fasta --start 0 --mode cut --revpr --pf CPRIMER --outname ",
              paste(f, "_MPC", sep = ''), " --log ",
              paste("/data/Assembly/", f, "_MPC_TRIM.log", sep = ''), sep = '')
      )
      #}
      
      
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
        paste("python3.12 .local/bin/CollapseSeq.py -s ",
              paste("/data/Assembly/", f, suffix, sep = '' ), 
              " -n 20 --inner --uf CPRIMER --cf VPRIMER --act set --outname ",
              paste(f, sep = ''), sep = '')
        
      )
      
      # Tell the system
      system(
        paste("python3.12 .local/bin/CollapseSeq.py -s ",
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
        paste("python3.12 .local/bin/SplitSeq.py group -s ",
              paste("/data/Assembly/", f, suffix, sep = '' ), 
              " -f DUPCOUNT --num 1 --outname ",
              paste(f, sep = ''), sep = '')
        
      )
      
      
      # Tell the system
      system(
        paste("python3.12 .local/bin/SplitSeq.py group -s ",
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
      
    }
  }
  
