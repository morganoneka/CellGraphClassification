library("matlab")
library("stringr")

# iterate over the diagnoses
groups = c("Chronic Pancreatitis", "IPMN", "MCN", "PanIN", "PDAC")

for (group in groups){
  output_dir=paste("/Users/morganoneka/Box/My Stuff/GraphClassification/GiottoOutput_Binarized/", group, sep="")
  
  files <- list.files(path=output_dir, pattern="txt")
  
  # iterate over files in each subfolder
  for (file in files){
    # read in table
    data_in <- read.table(fullfile(output_dir,file), sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
    
    # identify significant depletion 
    data_lower <- data_in[which(data_in$p_lower_orig <= 0.05),]
    data_lower$InxType <- "Depletion"
    
    # identify significant enrichment
    data_higher <- data_in[which(data_in$p_higher_orig <= 0.05),]
    data_higher$InxType <- "Enrichment"
    
    # combine depletion and enrichment info
    data_combo <- rbind(data_lower, data_higher)
    
    # subset
    data_combo <- data_combo[,c("cell_1", "cell_2", "enrichm", "InxType")]
    
    colnames(data_combo) <- c("source", "target", "weight", "type")
    
    
    write.table(data_combo, fullfile(output_dir, paste(paste(str_split_fixed(file, "_", 3)[,1:2], collapse="_"), ".csv", sep="")),
                row.names=FALSE, quote=FALSE, sep=",")
    
  }
}





