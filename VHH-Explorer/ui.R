# Load the ggplot2 package which provides
# the 'mpg' dataset.
library(ggplot2)

tmp <- read.csv( "~/VHH/data/Cluster-DBs/cluster_database/9_L3_sh2_S9_min_threshold_clone-pass-consensuses15Mar24.csv", 
                 header = T, stringsAsFactors = F)

colnames(tmp) <- c("File", "clone_id", "CDR3-Consensus", "DNA-Seqs", "Unique-DNA-Seqs", "Unique-CDR3s", "CDR3-Entropy")

fluidPage(
  titlePanel("VHH-Explorer"),
  
  # Create a new Row in the UI for selectInputs
 fluidRow(
    column(4,
           selectInput("File",
                       "Illumina Run:",
                       c("All",
                         unique(as.character(tmp$File))))
           

    ),
    column(4),
    column(4, downloadButton("downloadData", "Download Data (Ignores Search/Order)"))

  ),
  # Create a new row for the table.
  DT::dataTableOutput("table")
)