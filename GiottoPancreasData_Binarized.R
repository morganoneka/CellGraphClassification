library("matlab")
library("optparse")
library("Giotto")
library("reshape2")
library("ggtern")
library("akima")
library("ggalluvial")
library("igraph")
library("stringr")

groups = c("Chronic Pancreatitis", "IPMN", "MCN", "PanIN", "PDAC")
rel_col = c("Path","Sample.Name","Phenotype","Cell.ID","Total.Cells","Cell.X.Position","Cell.Y.Position","CD8","FoxP3","PD.L1","CD4","Treg","APC","Epithelial","Tcell")



  # for each group
  for (group in groups){
    print(group)
    file_dir=paste("/Users/morganoneka/Box/Phenotype Spreadsheets for Each Dx/", group, sep="")
    output_dir=paste("/Users/morganoneka/Box/My Stuff/GraphClassification/GiottoOutput_Binarized/", group, sep="")
    files <- list.files(path=file_dir, pattern="txt")
    
    # patients <- unique(str_split_fixed(files, "_", 2)[,1])
    
    # read in files from each patient
    for (file in files){
      
      # print(paste(group, patient))
      patient_and_sample = paste(str_split_fixed(file,"_",3)[,1:2], collapse="_")
      
      if(file.exists(fullfile(output_dir, paste(patient_and_sample, "interaction_enrichment.txt", sep="_")))){
        print("File Exists")
        next
      }
      
      # for (pfile in patient_files){
      data_in <- read.table(fullfile(file_dir,file), header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
      
      all_markers_present = TRUE
      for (ct in c("CD8", "FoxP3", "PD.L1", "CD4", "Treg","APC","Epithelial", "Tcell")){
        if (! ct %in% colnames(data_in)){
          all_markers_present = FALSE
        }
      }
      
      if (!all_markers_present){
        print("Missing a marker")
        next
      }
      
      data_in <- data_in[,rel_col]
      
      #TODO: get phenotype 
      data_in$CellType <- unlist(lapply(1:nrow(data_in), function(x){
        markers = c()
        for (celltype in c("CD8", "FoxP3", "PD.L1", "CD4", "Treg","APC","Epithelial", "Tcell")){
          if (data_in[x,celltype] == "pos"){
            markers = c(markers, paste(celltype,"+", sep=""))
          }
          
          
        }
        
        if (length(markers) == 0){
          return("Negative")
        } else{
          return(paste(markers, sep=" ", collapse=" "))
          
        }
      }))
      
      # ggplot(data_in, aes(x=Cell.X.Position, y=Cell.Y.Position)) + geom_point() 
      # ggsave(fullfile(output_dir, paste(patient, "cellplot.png", sep="")), width=7, height=7, device = "png")
      
      non_gene_columns <- c("cell_ID", "Cell.X.Position", "Cell.Y.Position", "CellType")
      x_col = "Cell.X.Position"
      y_col = "Cell.Y.Position"
      
      data_in[data_in=="neg"] = 0
      data_in[data_in=="pos"] = 1
      
      data_in$CD8 <- as.numeric(data_in$CD8)
      data_in$FoxP3 <- as.numeric(data_in$FoxP3)
      data_in$PD.L1 <- as.numeric(data_in$PD.L1)
      data_in$CD4 <- as.numeric(data_in$CD4)
      data_in$Treg <- as.numeric(data_in$Treg)
      data_in$APC <- as.numeric(data_in$APC)
      data_in$Tcell <- as.numeric(data_in$Tcell)
      
      
      
      print("Making Giottto object")
      # go <- createGiottoObject(transpose(data_in[,colnames(data_in)[!(colnames(data_in) %in% non_gene_columns)]]), data_in[,c(x_col,y_col)])
      go <- createGiottoObject(transpose(data_in[,c("CD8", "FoxP3", "PD.L1", "CD4", "Treg","APC","Epithelial", "Tcell")]), data_in[,c(x_col,y_col)])
      print("Giotto object made")
      
      with_annot <- addGeneMetadata(go, data.frame(gene_ID = c("CD8", "FoxP3", "PD.L1", "CD4", "Treg","APC","Epithelial", "Tcell")))
      with_network <- createSpatialNetwork(with_annot)
      with_spatialgrid <- createSpatialGrid(with_network, sdimx_stepsize = 50, sdimy_stepsize = 50)
      print("Spatial Grid made")
      
      
      with_phenotype <- with_spatialgrid
      with_phenotype@cell_metadata$phenotype <- data_in$CellType
      print("Reassigned phenotype")
      
      # tryCatch(, error = next)
      cell_prox_enrich <- cellProximityEnrichment(with_phenotype, cluster_column="phenotype")
      
      
      # for (i in 1:10){
      #   tryCatch(cell_prox_enrich <- cellProximityEnrichment(with_phenotype, cluster_column="phenotype"),
      #            error = function(e){
      #              
      #            }
      #            )
      # }
      # 
      print("Proximity enrichment")
      
      # plot cell enrichment
      #TODO: make sure plotting is working
      cell_prox_enrich$enrichm$cell_1 <- str_split_fixed(cell_prox_enrich$enrichm$unified_int,"--",2)[,1]
      cell_prox_enrich$enrichm$cell_2 <- str_split_fixed(cell_prox_enrich$enrichm$unified_int,"--",2)[,2]
      
      to_write_out <- cell_prox_enrich$enrichm[,c("enrichm", "p_higher_orig", "p_lower_orig", "PI_value", "cell_1", "cell_2")]
      
      write.table(to_write_out[order(-abs(to_write_out$enrichm)),],file=fullfile(output_dir, paste(patient_and_sample, "interaction_enrichment.txt", sep="_")), sep="\t", row.names=FALSE, quote=FALSE)
    }
  }















