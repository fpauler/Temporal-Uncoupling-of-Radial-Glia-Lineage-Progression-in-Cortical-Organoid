# ========================================================================================
# step 12
# script to prepare diBella et al data to be used as reference for MADM-CloneSeq analysis
# ========================================================================================

library (Seurat)
library (data.table)
library (ggplot2)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up


# metaData_scDevSC.txt downloaded from single cell portal
# at the time of download (03/2022):
# https://singlecell.broadinstitute.org/single_cell/study/SCP1290/molecular-logic-of-cellular-diversification-in-the-mammalian-cerebral-cortex
meta_data <- fread(paste(nobackup_base, "DiBella/metaData_scDevSC.txt", sep=""))
meta_data <- as.data.frame(meta_data)
meta_data <- meta_data[2:nrow(meta_data),]

# remove spaces from ttype IDs - makes problems with saving intermediate files
meta_data$New_cellType <- gsub( pattern = " ", replacement = "_", x = meta_data$New_cellType)

# files come from GSE153162_RAW.tar downloaded from GEO
# link al files to developmental age
file_list <- c("GSM4635072_E11_5_filtered_gene_bc_matrices_h5.h5", "GSM4635073_E12_5_filtered_gene_bc_matrices_h5.h5",
               "GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5", "GSM4635075_E14_5_filtered_gene_bc_matrices_h5.h5",
               "GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5", "GSM4635077_E16_filtered_gene_bc_matrices_h5.h5",
                "GSM5277844_E17_5_filtered_feature_bc_matrix.h5", "GSM4635078_E18_5_S1_filtered_gene_bc_matrices_h5.h5",
               "GSM4635079_E18_S3_filtered_gene_bc_matrices_h5.h5", "GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5",
               "GSM4635081_P1_S2_filtered_gene_bc_matrices_h5.h5", "GSM5277843_E10_v1_filtered_feature_bc_matrix.h5" )

age_list <- c("E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18_S1", "E18_S3", "P1_S1", "P1", "E10")
ct_prefix <- c("E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E18", "P1", "P1", "E10")

cell_types_to_include <- c("Apical_progenitors", "Intermediate_progenitors", "Immature_neurons", 
                           "Migrating_neurons", "SCPN", "CThPN", "DL_CPN", "UL_CPN", "Layer_4", "NP", "Layer_6b")

seurat.obj.lst <- list()

# use selected developmental ages here
for ( i in c(12, 1:7, 9, 10 )) {

  fileName <- file_list[i]
  
  sc.data <- Read10X_h5(filename = paste(nobackup_base, "DiBella/", fileName, sep=""))

  sub_meta_data <- meta_data[ which( meta_data$biosample_id ==  age_list[i] & meta_data$New_cellType %in% cell_types_to_include ), ]
	
  # remove cell types with less than 100 cells
  cell_types_to_remove <- c( names ( which( table(sub_meta_data$New_cellType) < 50 ) ) )

  cellIDs <- sub_meta_data$NAME
	
  cellIDs <- gsub( pattern = "-1", replacement = "", x = cellIDs )
	
  cellIDs <- sapply(cellIDs, function (x) {
    tmp <- strsplit(x = x, split = "_")[[1]]
    return ( tmp[length(tmp)] )
  })
  
  cellIDs <- paste( cellIDs, 1, sep="-" )
  
  rownames(sub_meta_data) <- cellIDs
	
  colIDX <- which( colnames(sc.data) %in% cellIDs)
  sc.data <- sc.data[,colIDX]
  sub_meta_data <- sub_meta_data[colnames(sc.data),]
  sub_meta_data$age <- ct_prefix[i]
    
  message ( paste( age_list[i], fileName, nrow( sub_meta_data ), ncol(sc.data), nrow(sc.data), ct_prefix[i] ) )
	
  seurat.obj.lst[[ age_list[i] ]] <-  CreateSeuratObject(counts = sc.data, project = age_list[i], min.cells = 3, min.features = 200)
  seurat.obj.lst[[ age_list[i] ]] <- AddMetaData( object = seurat.obj.lst[[ age_list[i] ]], metadata = sub_meta_data[, c("age", "New_cellType", "Phase")] )
	
  seurat.obj.lst[[ age_list[i] ]] <- NormalizeData( seurat.obj.lst[[ age_list[i] ]] )
  seurat.obj.lst[[ age_list[i] ]] <- FindVariableFeatures( seurat.obj.lst[[ age_list[i] ]], nfeatures = 3000 )
  
  seurat.obj.lst[[ age_list[i] ]] <- ScaleData( seurat.obj.lst[[ age_list[i] ]] )
  seurat.obj.lst[[ age_list[i] ]] <- RunPCA( seurat.obj.lst[[ age_list[i] ]], assay = "RNA", npcs = 30)
  
  elbow_plot <- ElbowPlot( object = seurat.obj.lst[[ age_list[i] ]], ndims = 30 )
  
  seurat.obj.lst[[ age_list[i] ]] <- FindNeighbors(seurat.obj.lst[[ age_list[i] ]], dims = 1:15 )
  seurat.obj.lst[[ age_list[i] ]] <- FindClusters(seurat.obj.lst[[ age_list[i] ]], resolution = 0.8)
  
  seurat.obj.lst[[ age_list[i] ]] <- RunUMAP(seurat.obj.lst[[ age_list[i] ]], reduction = "pca", dims = 1:15, 
                           n.neighbors = 30, metric = "euclidean", 
                           min.dist = 0.3, return.model = F, seed.use = 2401)
  
  umap_plot <- DimPlot( seurat.obj.lst[[ age_list[i] ]], label =T ) + NoLegend()
  
  # save for sanity check
  overview_plot <- elbow_plot + umap_plot
  ggsave( filename = paste(main_base, "plots/QC/diBella_", age_list[i], ".pdf",  sep=""), plot = overview_plot, width = 6, height = 4)
}

saveRDS( object = seurat.obj.lst, file = paste(base, "RDS_files/DiBella.seurat.obj.list.E10_to_P1.RDS", sep=""))


