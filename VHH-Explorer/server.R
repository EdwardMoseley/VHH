# Load the ggplot2 package which provides
# the 'mpg' dataset.
library(ggplot2)

#root_dir <- "~/VHH/data/Cluster-DBs/cluster_database/"


function(input, output) {
  tmp <- read.csv("~/VHH/data/Cluster-DBs/cluster_database/9_L3_sh2_S9_min_threshold_clone-pass-consensuses15Mar24.csv", 
                   header = T, stringsAsFactors = F)
  
  colnames(tmp) <- c("File", "clone_id", "CDR3-Consensus", "DNA-Seqs", "Unique-DNA-Seqs", "Unique-CDR3s", "CDR3-Entropy")
  
  # Filter data based on selections
  output$table <- DT::renderDataTable(DT::datatable({
    
    if (input$File != "All") {
      tmp <- tmp[tmp$File == input$File,]
    }
    tmp
  }))
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      "filtered_data.csv"  # Static filename
    },
    
    content = function(file) {
      if (input$File != "All") {
        write.csv(tmp[tmp$File == input$File,], file, row.names = F)
      } else {
        write.csv(tmp, file, row.names = FALSE)
      }
    }
  )
  
}