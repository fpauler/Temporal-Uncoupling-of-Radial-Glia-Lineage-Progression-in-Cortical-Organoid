# =================================================================
# step 04
# this script is for preparing QC figures after Gfp detection
# main task is to do a QC for Gfp marking relevant cell types
# =================================================================


library (Seurat)
library (dplyr)
library (ggplot2)
library (clusterProfiler)
library (org.Mm.eg.db)
library (openxlsx)
library (pheatmap)

set.seed( 2401 )

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

# Gfp cluster information, comes from: step03
clusters2remove <- list()
clusters2remove[["E10"]] <- setdiff( c(0:21), c("10", "1", "6", "16", "11", "7", "9") )
clusters2remove[["E13"]] <- c("13", "14", "15", "17", "19", "20", "21", "22", "23")
clusters2remove[["E16"]] <- c("10", "12", "13", "14", "16", "17", "18", "19", "20")
clusters2remove[["P0"]] <- c("12", "13", "14", "15", "16", "17", "19", "20", "21")

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

# these are the samples I use from the LaManno data
# focus only on relevant samples to reduce variability
# NOTE: for E11, I included more brain regions as dissection might be less precise at that age
sampleID_list <- list()
sampleID_list[[ "e11.0" ]] <- c( "10X40_3", "10X40_4", "10X40_5", "10X40_6" )
sampleID_list[[ "e13.5" ]] <- c( "10X14_4" )
sampleID_list[[ "e16.5" ]] <- c( "10X17_4" )
sampleID_list[[ "e18.0" ]] <- c( "10X32_2" )
   
age_list <- c("E10", "E13", "E16", "P0")
names( age_list ) <- c("e11.0","e13.5", "e16.5", "e18.0")

  # read embryo data
  embryo_obj <- readRDS( file =  paste( nobackup_base,  "/RDS_files/seurat.embryo.RDS", sep="" ) )
  
  # load loom libraries
  library (loomR)
  # process loom file
  lfile <- connect(filename = paste( nobackup_base, "/dev_all.loom", sep=""), mode = "r+", skip.validate = T)
  
  # loop through developmental ages
  for ( age in names( sampleID_list )[2:4] ) {

    # extract reference cells
    cellIDX <- which( lfile$col.attrs$Age[] == age )
    cellIDX <- cellIDX[which(!(duplicated(lfile$col.attrs$CellID[cellIDX]) | duplicated(lfile$col.attrs$CellID[cellIDX], fromLast = T)))]
    
    sampleIDX <- which( lfile$col.attrs$SampleID[] %in% sampleID_list[[ age ]] )
    
    cellIDX <- intersect( cellIDX, sampleIDX)
    
    data.subset <- lfile[["matrix"]][ cellIDX, ]
    
    #give cellID as row
    rownames(data.subset) <- lfile$col.attrs$CellID[cellIDX]
    #give Gene name as col - remove duplicates genes
    colnames(data.subset) <- lfile$row.attrs$Gene[]
    
    data.subset <- data.subset[ , which(!(duplicated(colnames(data.subset)) | duplicated(colnames(data.subset), fromLast = T))) ]
    
    meta_data <- list()
    meta_data[["cell_type"]] <- lfile$col.attrs$Class[cellIDX]
    names(meta_data[["cell_type"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["Age"]] <- lfile$col.attrs$Age[cellIDX]
    names(meta_data[["Age"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["Tissue"]] <- lfile$col.attrs$Tissue[cellIDX]
    names(meta_data[["Tissue"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["Location"]] <- lfile$col.attrs$Location_E9_E11[cellIDX]
    names(meta_data[["Location"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["clusterID"]] <- lfile$col.attrs$ClusterName[cellIDX]
    names(meta_data[["clusterID"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["subclass"]] <- lfile$col.attrs$Subclass[ cellIDX ]
    names(meta_data[["subclass"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    # read embryo data from this study
    
    # focus on a single dev. time point
    # this junk of code reproduces the UMAP used to identify GFP positive clusters
    tmp_embryo_obj <- subset( embryo_obj, idents = age_list[[ age ]] )
    
    tmp_embryo_obj <- NormalizeData(tmp_embryo_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    tmp_embryo_obj <- FindVariableFeatures(tmp_embryo_obj, selection.method = "vst", nfeatures = 2000)
    
    tmp_embryo_obj <- ScaleData(tmp_embryo_obj)
    tmp_embryo_obj <- RunPCA(tmp_embryo_obj, features = VariableFeatures(object = tmp_embryo_obj))
    
    tmp_embryo_obj <- FindNeighbors(tmp_embryo_obj, dims = 1:25)
    tmp_embryo_obj <- FindClusters(tmp_embryo_obj, resolution = 1)
    
    small_clusters <- names ( which( table( tmp_embryo_obj@meta.data$seurat_clusters) < 100 ) )
    
    # add information about GFP expression for supplemental figure
    tmp_embryo_obj@meta.data$GFP_cluster <- "YES"
    tmp_embryo_obj@meta.data$GFP_cluster[ which( as.character(tmp_embryo_obj@meta.data$seurat_clusters) %in% c( clusters2remove[[ age_list[[ age ]] ]], small_clusters) ) ] <- "NO"
    
    # this is the base UMAP that was used to filter for GFP cells
    tmp_embryo_obj <- RunUMAP(tmp_embryo_obj, dims = 1:25)
    
    # get raw counts to prepare new objects with matching gene sets
    embryo_counts <- GetAssayData(object = tmp_embryo_obj, assay = "RNA", layer = "counts")
      
    # determin common genes - use of identical gene lists is necessary for later integration
    common_genes <- intersect(rownames(embryo_counts), colnames(data.subset))
      
    thisStudy_Seurat <- CreateSeuratObject(counts = embryo_counts[common_genes, ], project = paste(age, "_mTmG", sep=""), min.cells = 3, min.features = 200)
    LaManno_Seurat <- CreateSeuratObject(counts = as.matrix(t(data.subset[, common_genes ])), project = paste(age, "_LaManno", sep=""), min.cells = 3, min.features = 200)
      
    #add these labels to the Seurat object
    for (i in names(meta_data)){
      LaManno_Seurat <- AddMetaData(
        object = LaManno_Seurat,
        metadata = meta_data[[i]],
        col.name = i
      )  
    }
      
    # Process LaManno reference data
      
    # CellCycleScoring 
    LaManno_Seurat <- NormalizeData( LaManno_Seurat )
    LaManno_Seurat <- CellCycleScoring (
      object = LaManno_Seurat,
      g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
      s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
    )
      
    # regress out cell cycle effects
    LaManno_Seurat$CC.Difference <- LaManno_Seurat$S.Score - LaManno_Seurat$G2M.Score
    LaManno_Seurat <- SCTransform(LaManno_Seurat, vars.to.regress = c("CC.Difference"), verbose = FALSE)
    LaManno_Seurat <- RunPCA(LaManno_Seurat, assay = "SCT")
    
    LaManno_Seurat <- RunUMAP(LaManno_Seurat, reduction = "pca", dims = 1:15, reduction.name = "umap")
      
    # process data from this study for integration
    # CellCycleScoring 
    thisStudy_Seurat <- NormalizeData(thisStudy_Seurat)
    thisStudy_Seurat <- CellCycleScoring (
      object = thisStudy_Seurat, 
      g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
      s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
    )
      
    # regress out cell cycle effects
    thisStudy_Seurat$CC.Difference <- thisStudy_Seurat$S.Score - thisStudy_Seurat$G2M.Score
    thisStudy_Seurat <- SCTransform(thisStudy_Seurat, vars.to.regress = c("CC.Difference"), verbose = FALSE)
    thisStudy_Seurat <- RunPCA( thisStudy_Seurat, assay = "SCT", npcs = 30)
      
    # do the integration / label transfer
    anchors <- FindTransferAnchors(reference = LaManno_Seurat, query = thisStudy_Seurat, 
                                   normalization.method = "SCT",
                                   reference.reduction = "pca", dims = 1:15)
    
    # transfer different labels
    predictions_tissue <- TransferData(anchorset = anchors, refdata = LaManno_Seurat$Tissue, dims = 1:15)
    predictions_clusterID <- TransferData(anchorset = anchors, refdata = LaManno_Seurat$clusterID, dims = 1:15)
    predictions_subclass <- TransferData(anchorset = anchors, refdata = LaManno_Seurat$subclass, dims = 1:15)
      
    predictions <- cbind( predictions_tissue[,1,drop=F], predictions_clusterID[,1,drop=F], predictions_subclass[,1,drop=F] )
    colnames(predictions) <- c("Tissue", "clusterID", "subclass")
      
    # add the labels to the original Seurat object
    tmp_embryo_obj <- AddMetaData(tmp_embryo_obj, metadata = predictions)

    # assign a subclasslabel to each seurat cluster
    # this is done by the largest fraction of subclass label 
    frac_Subclass <- data.frame()
    for (cluster in unique(tmp_embryo_obj$seurat_clusters)) {
      tmp_df <- as.data.frame( table( tmp_embryo_obj@meta.data$subclass[ which(tmp_embryo_obj$seurat_clusters == as.character(cluster)) ] ) )
      tmp_df$frac <- tmp_df$Freq / sum(tmp_df$Freq)
      tmp_df$cluster <- as.character(cluster) 
      frac_Subclass <- rbind( frac_Subclass, tmp_df[which(tmp_df$frac == max(tmp_df$frac)),] )
    }
      
    # lookup table to assign Subclass labels
    subclass_lookup <- sapply( as.character( unique(frac_Subclass$cluster) ), function(cluster){
      paste( frac_Subclass[ which(frac_Subclass$cluster == cluster), "Var1"], collapse = "\n" )
    })
      
    tmp_embryo_obj@meta.data$SubclassID <- subclass_lookup[ as.character(tmp_embryo_obj$seurat_clusters) ]
    tmp_embryo_obj@meta.data$Class <- linnarsson_annotation_lookup[ tmp_embryo_obj$clusterID, "Class" ]
      
    # create a combined annotation
    tmp_embryo_obj@meta.data$SubclassID_class <- paste( tmp_embryo_obj@meta.data$SubclassID, tmp_embryo_obj@meta.data$Class, sep="_" )
      

    # investigate  GABAergic neurons for markers of OBNB
    if (age %in% c("e16.5", "e18.0")) {
      
      Idents( tmp_embryo_obj) <- "SubclassID"
      tmp_inh_neurons <- subset( x = tmp_embryo_obj, idents = "Forebrain GABAergic" )
      Idents(tmp_inh_neurons) <- "GFP_cluster"
      
      this_obnb_markers <- FindMarkers(object = tmp_inh_neurons, ident.1 = "YES", ident.2 = "NO")
      #set of markers from here: https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2022.843794/full
      # and here: https://anatomypubs.onlinelibrary.wiley.com/doi/full/10.1002/ar.22733
      deg2plot <- this_obnb_markers[c("Lhx2", "Arx", "Pax6", "Tbr1", "Eomes", "Gsx2"),]
      deg2plot$gene <- rownames( deg2plot )
      deg2plot$gene <- factor( deg2plot$gene, levels = deg2plot$gene[ order(deg2plot$avg_log2FC, decreasing = T)] )
      deg_plot <- ggplot( deg2plot, aes(x=gene, y=avg_log2FC)) + geom_bar( stat = "identity" ) + theme_classic() + 
      			  ggtitle("OBNB markers, DEG Forebrain GABAergic GFP vs non-GFP")
      ggsave( filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_QC_DEG_OBNB.pdf", sep=""), plot = deg_plot )
      
    }
    
    # for the GFP fraction plot - only use annotations with >1% abundance in whole data set
    subclassID_list <- names( which( table(tmp_embryo_obj$SubclassID) / nrow(tmp_embryo_obj@meta.data) > 0.01 ) )
    Class_list <- names( which( table(tmp_embryo_obj$Class) / nrow(tmp_embryo_obj@meta.data) > 0.01 ) )
    SubclassID_class_list <- names( which( table(tmp_embryo_obj$SubclassID_class) / nrow(tmp_embryo_obj@meta.data) > 0.01 ) )

    # do pseudobulk expression analysis of key genes - not used for the paper
    emx1Expr <- AggregateExpression( tmp_embryo_obj, features = unique( c("GFP", "Emx1", "Top2a", VariableFeatures(tmp_embryo_obj) ) ), 
                                     group.by = "SubclassID_class", return.seurat = T, normalization.method = "LogNormalize", scale.factor = 10000)
 
    pheatmap( emx1Expr@assays$RNA$scale.data[c("Emx1", "Top2a", "GFP"), gsub( "_",  "-", SubclassID_class_list)], scale = "none", 
              filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_QC_ExpressionHeatmap.pdf", sep=""))
    
    # prepare some UMAPs and Feature plots - not all used for figures
    gfp_class_umap <- DimPlot( tmp_embryo_obj, reduction = "umap", group.by = "GFP_cluster", label = F ) + ggtitle ("Gfp") + 
      scale_color_manual( values = c( "YES" = "black", "NO" = "grey80"))
    ggsave ( filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_QC_gfp.pdf", sep=""), plot = gfp_class_umap, width = 8, height = 8 )
    
    emx1_expression_umap <- FeaturePlot( tmp_embryo_obj, features = "Emx1", order = T )
    ggsave ( filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_QC_emx1.pdf", sep=""), plot = emx1_expression_umap, width = 8, height = 8 )
    
    tissue_umap <- DimPlot( tmp_embryo_obj, reduction = "umap", group.by = "SubclassID", label = T ) + ggtitle ("Tissue origin")
    ggsave ( filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_QC_SubclassID.pdf", sep=""), plot = tissue_umap, width = 10, height = 8 )
    
    cellType_umap <- DimPlot( tmp_embryo_obj, reduction = "umap", group.by = "Class", label = T ) + ggtitle ("Class")
    ggsave ( filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_QC_Class.pdf", sep=""), plot = cellType_umap, width = 10, height = 8 )
        
    seurat_cluster_umap <- DimPlot( tmp_embryo_obj, reduction = "umap", group.by = "seurat_clusters", label = T ) + 
      ggtitle ("Seurat clusters") + NoLegend()
    ggsave ( filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_QC_SeuratCLusters.pdf", sep=""), plot = seurat_cluster_umap, width = 8, height = 8 )
    
    # GFP cell frequency plot - SubclassID
    gfp_frac_df <- data.frame()
    for (subclass in subclassID_list ) {
      tmp_df <- as.data.frame( table( tmp_embryo_obj$GFP_cluster[ which (tmp_embryo_obj$SubclassID == subclass ) ] ) )
      tmp_df$subclass <- subclass
      gfp_frac_df <- rbind( gfp_frac_df, tmp_df)
    }
    
    # prepare an order based on cell type abundance
    sum_freq <- gfp_frac_df %>%
      group_by(subclass) %>%
      summarize(
        sum = sum(Freq)
      )
    
    # calculate the relative abundance
    total_vec <- sum_freq$sum
    names( total_vec ) <- sum_freq$subclass
    gfp_frac_df$total <- total_vec[ gfp_frac_df$subclass ]
    gfp_frac_df$rel <- gfp_frac_df$Freq / gfp_frac_df$total
    
    gfp_frac_df$subclass <- factor( gfp_frac_df$subclass, levels = sum_freq$subclass[ order(sum_freq$sum) ] )
    
    gfp_frac_plot <- ggplot( gfp_frac_df, aes(x=subclass, y=rel, fill=Var1, group = subclass)) + 
      geom_bar(stat="identity") +  scale_fill_manual( values = c( "YES" = "black", "NO" = "grey80")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) 
    ggsave ( plot = gfp_frac_plot, filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_GFPfrac.Subclass.pdf", sep=""), width = 8, height = 5 )
    
    # GFP cell frequency plot - Class
    gfp_frac_df <- data.frame()
    for (class in Class_list ) {
      tmp_df <- as.data.frame( table( tmp_embryo_obj$GFP_cluster[ which (tmp_embryo_obj$Class == class ) ] ) )
      tmp_df$class <- class
      gfp_frac_df <- rbind( gfp_frac_df, tmp_df)
    }
    
    # prepare an order based on cell type abundance
    sum_freq <- gfp_frac_df %>%
      group_by(class) %>%
      summarize(
        sum = sum(Freq)
      )
    
    # calculate the relative abundance
    total_vec <- sum_freq$sum
    names( total_vec ) <- sum_freq$class
    gfp_frac_df$total <- total_vec[ gfp_frac_df$class ]
    gfp_frac_df$rel <- gfp_frac_df$Freq / gfp_frac_df$total
    
    gfp_frac_df$class <- factor( gfp_frac_df$class, levels = sum_freq$class[ order(sum_freq$sum) ] )
    
    gfp_frac_plot <- ggplot( gfp_frac_df, aes(x=class, y=rel, fill=Var1, group = class)) + 
      geom_bar(stat="identity") +  scale_fill_manual( values = c( "YES" = "black", "NO" = "grey80")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) 
    ggsave ( plot = gfp_frac_plot, filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_GFPfrac.Class.pdf", sep=""), width = 8, height = 5 )
    
    # combined annotation
    tmp_counts <- table(f@meta.data$SubclassID_class)
    tmp_counts <- tmp_counts / sum(tmp_counts)
    group_vec <- names( which(tmp_counts > 0.03) )
    
    gfp_frac_df <- data.frame()
    for (group in group_vec ) {
      tmp_df <- as.data.frame( table( tmp_embryo_obj$GFP_cluster[ which (tmp_embryo_obj$SubclassID_class == group ) ] ) )
      tmp_df$group <- group
      gfp_frac_df <- rbind( gfp_frac_df, tmp_df)
    }
    
    # prepare an order based on cell type abundance
    sum_freq <- gfp_frac_df %>%
      group_by(group) %>%
      summarize(
        sum = sum(Freq)
      )
    
    # calculate the relative abundance
    total_vec <- sum_freq$sum
    names( total_vec ) <- sum_freq$group
    gfp_frac_df$total <- total_vec[ gfp_frac_df$group ]
    gfp_frac_df$rel <- gfp_frac_df$Freq / gfp_frac_df$total
    
    gfp_frac_df$group <- factor( gfp_frac_df$group, levels = sum_freq$group[ order(sum_freq$sum) ] )
    
    gfp_frac_plot <- ggplot( gfp_frac_df, aes(x=group, y=rel, fill=Var1, group = group)) + 
      geom_bar(stat="identity") +  scale_fill_manual( values = c( "YES" = "black", "NO" = "grey80")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) 
    ggsave ( plot = gfp_frac_plot, filename = paste(base, "/plots/Sup_Fig_5_6/embryo_", age, "_GFPfrac.Subclass_Class.pdf", sep=""), width = 8, height = 5 )
    
    # save for processing organoid data
    saveRDS( object = LaManno_Seurat, file = paste( nobackup_base,  "/RDS_files/LaManno.", age,".RDS", sep="" ) )
  }

