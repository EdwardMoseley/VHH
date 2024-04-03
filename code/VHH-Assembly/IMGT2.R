
# Arguments
args_list <- "S1_atleast-2.fmt7"# list.files("~/VHH/data")[grepl(".fmt7", list.files("~/VHH/data"))]
# args = commandArgs(trailingOnly=TRUE)

# Directory
directory <- "~/VHH/data/"


for (i in 1:length(args_list)){
  file_i <- args_list[i]#paste(directory, args_list[i], sep = '')
  file_s <- gsub(".fmt7", ".fasta", file_i)
  results <- gsub(".fastq", '.fmt7', args_list)
  
  
  command <- paste('docker run -v /Users/emoseley/VHH/data:/data immcantation/suite:4.4.0 MakeDb.py igblast -i ', 
                   paste('/data/', file_i, sep = ''), ' -s ',
                   paste('/data/', file_s, sep = ''), ' -r /data/imgt-data/imgt_gapped_alpaca_IGHV.fasta /data/imgt-data/imgt_gapped_alpaca_IGHD.fasta  /data/imgt-data/imgt_gapped_alpaca_IGHJ.fasta --outdir /data/ --failed --format airr --extended', sep = '')
  
  print(command)
  #system(command)
}




#
setwd(directory)

files <- list.files(directory)[grepl("db-", list.files(directory))]

tmp_dat <- read.csv(paste(directory, files[4],sep = ''), sep = '\t', header = T, stringsAsFactors = F)

for (i in 2:length(files)){
  tmp_dat <- rbind(tmp_dat, 
                   read.csv(paste(directory, files[i], sep = ''), sep = '\t', header = T, stringsAsFactors = F))
}

length(unique(tmp_dat$cdr3))