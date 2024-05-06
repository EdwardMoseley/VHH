## ###
## QueryingtoSource.R
## ###

fuzzy_query_source <- function(query_seq, source_db, var_name){
  # Weights correspond to being 99x as expensive to insert/delete
  distances <- stringdist::stringdist(query_seq, source_db[[var_name]], method = "lv", 
                                   weight = c(0.99, 
                                              0.99,
                                              0.01,
                                              0.01))
  source_db$distances <- distances
  
  source_db <- source_db[order(source_db$distances, decreasing = F),]
  
  return(source_db)
}




#root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/"

root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/translated_dbs/culled_databases/"
list.files(root_dir)[grepl("Cleaned", list.files(root_dir))]


# Create a list to use as a dictionary
files <- list()

# Populate key/values with files and their keys
files[["S13"]] <- "13_VecF1_R2_S13-Cleaned19Apr24.csv"
files[["S14"]] <- "14_VecF1_R2_S14-Cleaned19Apr24.csv"

files[["S19"]] <- "19_VecF1_R3_S19-Cleaned19Apr24.csv"
files[["S20"]] <- "20_VecF2_R3_S20-Cleaned19Apr24.csv"
files[["S21"]] <- "21_VecF4_R3_S21-Cleaned19Apr24.csv" 

#max_clone_count <- 0
#max_id <- "S0"

for (nm in names(files)){
  files[[nm]] <- read.csv(paste(root_dir, files[[nm]], sep = ''), header = T, stringsAsFactors = F)
}




tmp <- files[["S13"]]
tmp <- fuzzy_query_source("AAAQTYRYSTRSADDFSD", tmp, "junction_aa")
table(tmp[tmp$junction_aa == "AAAQTYRYSTRSADDFSD",]$clone_id)
sum(tmp[tmp$distances == 0,]$duplicate_count)
head(tmp)


tmp <- files[["S14"]]
tmp <- fuzzy_query_source("AAAQTYRYSTRSADDFSD", tmp, "junction_aa")
table(tmp[tmp$junction_aa == "AAAQTYRYSTRSADDFSD",]$clone_id)
sum(tmp[tmp$distances == 0,]$duplicate_count)
head(tmp)

tmp <- files[["S13"]]
tmp <- fuzzy_query_source("GTVVSGASFDS", tmp)
table(tmp[tmp$junction_aa == "GTVVSGASFDS",]$clone_id)
head(tmp)
sum(tmp[tmp$distances == 0,]$duplicate_count)


table(tmp[tmp$clone_id == "7479",]$junction_aa)







tmp <- files[["S14"]]
tmp <- fuzzy_query_source("GTVVSGASFDS", tmp)
sum(tmp[tmp$distances == 0,]$duplicate_count)
table(tmp[tmp$junction_aa == "GTVVSGASFDS",]$clone_id)
head(tmp)







test <- clonal_statistics(files[["S13"]])

test[(test$median_CDR3 == "GTVVSGASFDS"),]

test_two <- clonal_statistics(files[["S14"]])























































root_dir <- "~/VHH/data/Cluster-DBs/phylogenetic_dbs/"
list.files(root_dir)[grepl("min_threshold", list.files(root_dir))]


# Create a list to use as a dictionary
files <- list()

# Populate key/values with files and their keys
files[["S13"]] <- "13_VecF1_R2_S13_min_threshold_clone-pass.csv"
files[["S14"]] <- "14_VecF1_R2_S14_min_threshold_clone-pass.csv"

files[["S19"]] <- "19_VecF1_R3_S19_min_threshold_clone-pass.csv"
files[["S20"]] <- "20_VecF2_R3_S20_min_threshold_clone-pass.csv"
files[["S21"]] <- "21_VecF4_R3_S21_min_threshold_clone-pass.csv" 

#max_clone_count <- 0
#max_id <- "S0"

for (nm in names(files)){
  files[[nm]] <- read.csv(paste(root_dir, files[[nm]], sep = ''), header = T, stringsAsFactors = F)
}

