# ======================================================================================
# step 05
# this script is for preparing QC figures after Gfp detection step
# main task is to do a QC for Gfp marking relevant cell types
# this script had to be run multiple times (once for each dev. age) due to memory limits
# ======================================================================================

library (Seurat)
library (dplyr)
library (ggplot2)
library (clusterProfiler)
library (org.Mm.eg.db)
library( openxlsx )
library (pheatmap)

set.seed( 2401 )

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

# define clusters to include
# comes from: EB analysis, Step03_GFP_extract.R
clusters2include <- list()
clusters2include[["D8"]] <- c("2", "3", "9") 
clusters2include[["D13"]] <- c("0", "1", "2", "3", "4", "5")
clusters2include[["D20"]] <- c("0", "1", "2", "3", "4", "6", "8", "9", "10", "11")
clusters2include[["D25"]] <- c("0", "1", "2", "3", "4", "6", "7", "8", "10", "11",
                               "12", "13", "15" )

# cell cycle gene information
# read cell cycle genes to regress out cell cycle effects
cell_cycle_phase_gene <- read.table(file = paste(base, "other_data/Mus_musculus.csv", sep=""), sep = ",", header = T)
gene_conversion <- bitr(geneID = cell_cycle_phase_gene$geneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
cell_cycle_phase_gene <- merge(cell_cycle_phase_gene, gene_conversion, by.x="geneID", by.y="ENSEMBL")

# lookup table for Linnarsson annotation
# from http://mousebrain.org/development/celltypes.html obtained 09/2025
# NOTE: I removed one instance of NblM7, which was duplicated
linnarsson_annotation <- read.xlsx( xlsxFile = paste(base, "other_data/Linnarsson_cellTypes.xlsx", sep=""))

# prepare a lookup table and add Linnarsson annotation to Seurat object
# this is a broader cell type annotation
linnarsson_annotation_lookup <- linnarsson_annotation[ which(!is.na(linnarsson_annotation$Symbol)), c(1,3,4,5) ]
rownames(linnarsson_annotation_lookup) <- linnarsson_annotation_lookup$Symbol

# link EB age to embryonic age
# necessary to use the data saved from embryo data for this analysis
org_age_list <- c("D8", "D13", "D20", "D25")
names( org_age_list ) <- c("e11.0","e13.5", "e16.5", "e18.0")

organoid_cellID <- c("lineA1", "lineA2", "lineB1", "lineB2")

# determine the next dev age to analyse
# only necessary if script runs multiple times
for (i in 1:length(names( org_age_list ))) {
  age <- names( org_age_list )[i]
  if ( !file.exists( paste( nobackup_base,  "/RDS_files/MetaData.GfpDetect.", org_age_list[[ age ]],".RDS", sep="" ) ) ) {
    message( "Index ", i, " age: ", age, " needs to be processed next")
    break
  }
}

# change i for specific re-plotting needs!
# i <- 3
# age <- names( org_age_list )[i]

if ( i < 5 ) {
  
  #################
  ##
  # organoid data
  ##
  #################
  
  # reads the EB cell info for GFP classification and joining
  if (!file.exists( paste( nobackup_base,  "/RDS_files/seurat.", org_age_list[[ age ]], ".joint.RDS", sep="" ) )) {
    
    seurat.obj <- readRDS( file = paste( nobackup_base,  "/RDS_files/seurat.", org_age_list[[ age ]], ".RDS", sep="" ) )
    
    # prepare the base UMAP/clustering that GFP fitering is based on
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:25)
    seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:25)
    DimPlot(seurat.obj, group.by = "seurat_clusters", label=T ) + NoLegend()
    
    # add information about GFP expression 
    seurat.obj@meta.data$GFP_cluster <- "NO"
    seurat.obj@meta.data$GFP_cluster[ which( as.character(seurat.obj@meta.data$seurat_clusters) %in% clusters2include[[ org_age_list[[ age ]] ]] ) ] <- "YES"
    
    seurat.obj <- JoinLayers( seurat.obj )
    
    seurat.obj <- CellCycleScoring (
      object = seurat.obj, 
      g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
      s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
    )
    
    seurat.obj$CC.Difference <- seurat.obj$S.Score - seurat.obj$G2M.Score
    
    # some data sets are too big for my memory to process in one go - save and restart
    saveRDS( object = seurat.obj, file = paste( nobackup_base,  "/RDS_files/seurat.", org_age_list[[ age ]], ".joint.RDS", sep="" ) )
    
    stop("joining done - restart")
    
  } else {
    seurat.obj <- readRDS( file = paste( nobackup_base,  "/RDS_files/seurat.", org_age_list[[ age ]], ".joint.RDS", sep="" ) )
  }
  
  # since this analysis takes very long - allow option to skip step and add saved meta data directly
  if (!file.exists(filename = paste( nobackup_base,  "/RDS_files/MetaData.GfpDetect.", org_age_list[[ age ]],".RDS", sep="" ))) {
   
    seurat.obj <- SCTransform(seurat.obj, vars.to.regress = c("CC.Difference"), verbose = FALSE)
    
    seurat.obj <- RunPCA( seurat.obj, assay = "SCT", npcs = 30)
    
    # this data set comes from embryo analysis
    LaManno_Seurat <- readRDS( file = paste( nobackup_base,  "/RDS_files/LaManno.", age,".RDS", sep="" ) )
    
    anchors <- FindTransferAnchors(reference = LaManno_Seurat, query = seurat.obj, 
                                   normalization.method = "SCT",
                                   reference.reduction = "pca", dims = 1:20)
    
    predictions_tissue <- TransferData(anchorset = anchors, refdata = LaManno_Seurat$Tissue, dims = 1:20)
    predictions_clusterID <- TransferData(anchorset = anchors, refdata = LaManno_Seurat$clusterID, dims = 1:20)
    predictions_subclass <- TransferData(anchorset = anchors, refdata = LaManno_Seurat$subclass, dims = 1:20)
    
    predictions <- cbind( predictions_tissue[,1,drop=F], predictions_clusterID[,1,drop=F], predictions_subclass[,1,drop=F] )
    colnames(predictions) <- c("Tissue", "clusterID", "subclass")
    
    # add the prediction to the Seurat object
    seurat.obj <- AddMetaData(seurat.obj, metadata = predictions)
    
    # add class information based on Linnarsson annotation
    seurat.obj@meta.data$Class <- linnarsson_annotation_lookup[ seurat.obj$clusterID, "Class" ]
    
    # for next steps change back to RNA assay
    DefaultAssay( seurat.obj ) <- "RNA"
    
    if (age == "e11.0") {
      # only relevant for e11/D8
      # detailed analysis indicates that cluster 9 should be broken up - do here by increasing resolution
      seurat.obj <- FindClusters(seurat.obj, resolution = 0.6)
      # change GFP tag of cluster 19
      # this is a cluster of immature neurons that was initially not separated
      # however the analysis here clearly indicates that this should be done as cluster 19 appears to be 
      # GABAergic and NOT Emx1 lineage 
      length( which( seurat.obj$seurat_clusters == "19" ) )
      # 333 cells
      seurat.obj@meta.data$GFP_cluster[ which( seurat.obj$seurat_clusters == "19" ) ] <- "NO"
    }
    
    frac_Subclass <- data.frame()
    for (cluster in unique(seurat.obj$seurat_clusters)) {
      tmp_df <- as.data.frame( table( seurat.obj@meta.data$subclass[ which(seurat.obj$seurat_clusters == as.character(cluster)) ] ) )
      tmp_df$frac <- tmp_df$Freq / sum(tmp_df$Freq)
      tmp_df$cluster <- as.character(cluster) 
      frac_Subclass <- rbind( frac_Subclass, tmp_df[which(tmp_df$frac == max(tmp_df$frac)),] )
    }
    
    # lookup table to assign Subclass labels
    subclass_lookup <- sapply( as.character( unique(frac_Subclass$cluster) ), function(cluster){
      paste( frac_Subclass[ which(frac_Subclass$cluster == cluster), "Var1"], collapse = "\n" )
    })
    
    seurat.obj@meta.data$SubclassID <- subclass_lookup[ as.character(seurat.obj$seurat_clusters) ]
    
  } else {
    meta_data <- readRDS( file = paste( nobackup_base,  "/RDS_files/MetaData.GfpDetect.", org_age_list[[ age ]],".RDS", sep="" ) ) 
    columsToAdd <- setdiff( colnames(meta_data), colnames(seurat.obj@meta.data))
    seurat.obj <- AddMetaData( object = seurat.obj, metadata = meta_data[,columsToAdd])
  }
  
  # create a combined annotation
  seurat.obj@meta.data$SubclassID_class <- paste( seurat.obj@meta.data$SubclassID, seurat.obj@meta.data$Class, sep="_" )
  
  # for the GFP fraction plot - only use annotations with >1% abundance in whole data set
  subclassID_list <- names( which( table(seurat.obj$SubclassID) / nrow(seurat.obj@meta.data) > 0.01 ) )
  Class_list <- names( which( table(seurat.obj$Class) / nrow(seurat.obj@meta.data) > 0.01 ) )
  SubclassID_class_list <- names( which( table(seurat.obj$SubclassID_class) / nrow(seurat.obj@meta.data) > 0.01 ) )
  
  # investigate  GABAergic neurons for markers of OBNB
  if (age == "e18.0") {
    
    Idents( seurat.obj) <- "SubclassID"
    tmp_inh_neurons <- subset( x = seurat.obj, idents = "Forebrain GABAergic" )
    Idents(tmp_inh_neurons) <- "GFP_cluster"
    this_obnb_markers <- FindMarkers(object = tmp_inh_neurons, ident.1 = "YES", ident.2 = "NO")
    #set of markers from here: https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.843794/full
    # and here: https://anatomypubs.onlinelibrary.wiley.com/doi/full/10.1002/ar.22733
    deg2plot <- this_obnb_markers[c("Lhx2", "Arx", "Pax6", "Tbr1", "Eomes", "Gsx2"),]
    deg2plot$gene <- rownames( deg2plot )
    deg2plot$gene <- factor( deg2plot$gene, levels = deg2plot$gene[ order(deg2plot$avg_log2FC, decreasing = T)] )
    deg_plot <- ggplot( deg2plot, aes(x=gene, y=avg_log2FC)) + geom_bar( stat = "identity" ) + theme_classic() + ggtitle("OBNB markers, DEG Forebrain GABAergic GFP vs non-GFP")
    ggsave( filename = paste(base, "/plots/Sup_Fig_2_3/organoid_", org_age_list[[ age ]], "_QC_DEG_OBNB.pdf", sep=""), plot = deg_plot )
     
  }
  
  # do pseudobulk expression analysis of key genes - not used for the paper
  emx1Expr <- AggregateExpression( seurat.obj, features = unique( c("GFP", "Emx1", "Top2a", VariableFeatures(seurat.obj) ) ),
                                   group.by = "SubclassID_class", return.seurat = T, normalization.method = "LogNormalize", scale.factor = 10000)
  
  pheatmap( emx1Expr@assays$RNA$scale.data[c("Emx1", "Top2a", "GFP"), gsub( "_",  "-", SubclassID_class_list)], scale = "none", 
            filename = paste(base, "/plots/Sup_Fig_2_3/organoid_", org_age_list[[ age ]], "_QC_ExpressionHeatmap.pdf", sep=""))
  
  # prepare some UMAPs and Feature plots - not all used for figures
  gfp_class_umap <- DimPlot( seurat.obj, reduction = "umap", group.by = "GFP_cluster", label = F ) + ggtitle ("Gfp") + 
    scale_color_manual( values = c( "YES" = "black", "NO" = "grey80"))
  ggsave ( filename = paste(base, "/plots/Sup_Fig_2_3/organoid_", org_age_list[[ age ]], "_QC_gfp.pdf", sep=""), plot = gfp_class_umap, width = 8, height = 8 )
  
  emx1_expression_umap <- FeaturePlot( seurat.obj, features = "Emx1", reduction = "umap", order = T, min.cutoff = "q15" )
  ggsave ( filename = paste(base, "/plots/Sup_Fig_2_3/organoid_", org_age_list[[ age ]], "_QC_emx1.pdf", sep=""), plot = emx1_expression_umap, width = 8, height = 8 )
  
  tissue_umap <- DimPlot( seurat.obj, reduction = "umap", group.by = "SubclassID", label = T ) + ggtitle ("Tissue origin")
  ggsave ( filename = paste(base, "/plots/Sup_Fig_2_3/organoid_", org_age_list[[ age ]], "_QC_SubclassID.pdf", sep=""), plot = tissue_umap, width = 10, height = 8 )
  
  cellType_umap <- DimPlot( seurat.obj, reduction = "umap", group.by = "Class", label = T ) + ggtitle ("Class")
  ggsave ( filename = paste(base, "/plots/Sup_Fig_2_3/organoid_", org_age_list[[ age ]], "_QC_Class.pdf", sep=""), plot = cellType_umap, width = 10, height = 8 )
  
  seurat_cluster_umap <- DimPlot( seurat.obj, reduction = "umap", group.by = "seurat_clusters", label = T ) + 
    ggtitle ("Seurat clusters") + NoLegend()
  ggsave ( filename = paste(base, "/plots/Sup_Fig_2_3/organoid_", org_age_list[[ age ]], "_QC_SeuratCLusters.pdf", sep=""), plot = seurat_cluster_umap, width = 8, height = 8 )
  
  # save meta data for later usage
  if (!file.exists( paste( nobackup_base,  "/RDS_files/MetaData.GfpDetect.", org_age_list[[ age ]],".RDS", sep="" ) )) {
    saveRDS( object = seurat.obj@meta.data, file = paste( nobackup_base,  "/RDS_files/MetaData.GfpDetect.", org_age_list[[ age ]],".RDS", sep="" ) )
  }
  
} else {
  message( "all done!")
}

