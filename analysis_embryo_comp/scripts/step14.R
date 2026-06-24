# ========================================================
# step 14
# integration and analysis of Reference, EB/organoid and PatchSeq
# ========================================================

library (Seurat)
library (ggplot2)
library (simspec)
library (RANN)
library (patchwork)
library (SeuratWrappers)
library (monocle3)
library (GenBinomApps)
library (dplyr)
library (ggbeeswarm)
library (openxlsx)
library (rgl)
library (pheatmap)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

set.seed(2401)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

if (!file.exists(paste(nobackup_base, "RDS_files/Seurat_Neurons.RDS", sep=""))) {
  
  Seurat_neuron <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.RDS", sep="") )
  
  DimPlot( Seurat_neuron, group.by = "cellType")
  
  Idents(Seurat_neuron) <- "cellType"
  Seurat_neuron <- subset(Seurat_neuron, idents = c("RGP", "IP", "iN"))
  Idents( Seurat_neuron ) <- "group"
  Seurat_neuron <- subset( Seurat_neuron, idents = "EB" )
  Idents( Seurat_neuron ) <- "age"
  Seurat_neuron <- subset( Seurat_neuron, idents = c("D13", "D20", "D25") )
  
  DimPlot( Seurat_neuron, group.by = "cellType")
  
  # An object of class Seurat 
  # 34420 features across 64647 samples within 2 assays 
  # Active assay: sketch (17210 features, 5000 variable features)
  # 2 layers present: counts, data
  # 1 other assay present: RNA
  # 3 dimensional reductions calculated: pca, css, umap
  
  Seurat_neuron@assays$RNA@layers$data <- NULL
  Seurat_neuron@assays$RNA@layers$scale.data <- NULL
  Seurat_neuron@reductions$pca <- NULL
  Seurat_neuron@reductions$css <- NULL
  Seurat_neuron@reductions$umap <- NULL
  
  DefaultAssay(Seurat_neuron) <- "RNA"
  
  Seurat_neuron@assays$sketch <- NULL
  
  # An object of class Seurat 
  # 17210 features across 64647 samples within 1 assay 
  # Active assay: RNA (17210 features, 5000 variable features)
  # 1 layer present: counts
  
  saveRDS( object = Seurat_neuron, file = paste(nobackup_base, "RDS_files/Seurat_Neurons.RDS", sep=""))
  
}

if (!file.exists(paste(nobackup_base, "RDS_files/Seurat_neuron_EB.RDS", sep="/"))) {
  
  Seurat_neuron <- readRDS( file = paste(nobackup_base, "RDS_files/Seurat_Neurons.RDS", sep=""))
  
  old_meta <- Seurat_neuron@meta.data
  
  Idents( Seurat_neuron ) <- "age"
  
  Seurat_neuron_EB <- list()
  for ( age in c("D13", "D20", "D25") ) {
    
    # start with a clean slate for this analysis
    Seurat_neuron_EB[[ age ]] <- CreateSeuratObject(counts = GetAssayData( object = subset( Seurat_neuron, idents = age ), layer = "counts"), 
                                                    project = "D13", min.cells = 3, min.features = 200)
    
    # add back relevant meta data
    Seurat_neuron_EB[[ age ]] <- AddMetaData(object = Seurat_neuron_EB[[ age ]], metadata = old_meta[,c("orig.ident", "cellType", "age", "group", "CC.Difference")])
    
    Seurat_neuron_EB[[ age ]] <- NormalizeData( Seurat_neuron_EB[[ age ]] )
    Seurat_neuron_EB[[ age ]] <- FindVariableFeatures( Seurat_neuron_EB[[ age ]], nfeatures = 5000 )
    
    Seurat_neuron_EB[[ age ]] <- ScaleData( Seurat_neuron_EB[[ age ]] )
    Seurat_neuron_EB[[ age ]] <- RunPCA( Seurat_neuron_EB[[ age ]], assay = "RNA", npcs = 30)
    
    elbow_plot <- ElbowPlot( Seurat_neuron_EB[[ age ]] )
    
    Seurat_neuron_EB[[ age ]] <- FindNeighbors(Seurat_neuron_EB[[ age ]], dims = 1:15, verbose = FALSE)
    Seurat_neuron_EB[[ age ]] <- FindClusters(Seurat_neuron_EB[[ age ]], resolution = 0.8, verbose = FALSE)
    
    Seurat_neuron_EB[[ age ]] <- RunUMAP(Seurat_neuron_EB[[ age ]], reduction = "pca", dims = 1:15, 
                                         n.neighbors = 30, metric = "cosine", 
                                         min.dist = 0.3, return.model = T, seed.use = 2401)
    
    umap_plot <- DimPlot( Seurat_neuron_EB[[ age ]], label =T ) + NoLegend()
    
    overview_plot <- elbow_plot + umap_plot
    ggsave( filename = paste(base, "plots/QC/Neurons_", age, ".pdf",  sep=""), plot = overview_plot, width = 6, height = 4)
    
  }
  
  # remove all unnecessary layers
  Seurat_neuron_EB_merged <- merge( x = Seurat_neuron_EB[[ 1 ]], y = Seurat_neuron_EB[ 2:3 ] )
  for (layer in Layers( Seurat_neuron_EB_merged )[ which( !grepl( pattern = "counts", x = Layers( Seurat_neuron_EB_merged ) ) )]) {
    Seurat_neuron_EB_merged@assays$RNA@layers[[ layer ]] <- NULL
  }
  
  saveRDS( object = Seurat_neuron_EB_merged, file = paste(nobackup_base, "RDS_files/Seurat_neuron_EB.RDS", sep="/") )
  
}


if (  !file.exists( paste(nobackup_base, "RDS_files/all_neurons_merged.RDS", sep="") ) ) {
  
  Seurat_neuron_EB_merged <- readRDS( file = paste(nobackup_base, "RDS_files/Seurat_neuron_EB.RDS", sep="/") )
  # An object of class Seurat 
  # 15651 features across 64647 samples within 1 assay 
  # Active assay: RNA (15651 features, 5000 variable features)
  # 3 layers present: counts.1, counts.2, counts.3
  
  colnames( Seurat_neuron_EB_merged@meta.data )
  # "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "cellType"        "age"             "group"          
  # "CC.Difference"   "RNA_snn_res.0.8" "seurat_clusters"
  
  diBella.obj.list <- readRDS( file = paste(nobackup_base, "RDS_files/DiBella.seurat.obj.list.E10_to_P1.RDS", sep=""))
  
  # names( diBella.obj.list )
  # "E10"    "E11"    "E12"    "E13"    "E14"    "E15"    "E16"    "E17"    "E18_S3" "P1_S1" 
  
  diBella.obj <- merge( x = diBella.obj.list[[3]], y = diBella.obj.list[ 4:length(diBella.obj.list)])
  # remove all unnecessary layers
  for (layer in Layers( diBella.obj )[ which( !grepl( pattern = "counts", x = Layers( diBella.obj ) ) )]) {
    diBella.obj@assays$RNA@layers[[ layer ]] <- NULL
  }
  
  PatchSeq_Seurat <- readRDS( file = paste(nobackup_base, "RDS_files/Seurat_PatchSeq.RDS", sep="/") )
  
  ######################################################
  # remove 2 cells with uncertain clone type assignment
  ######################################################
  
  Idents(PatchSeq_Seurat) <- "clone_type"
  PatchSeq_Seurat <- subset(PatchSeq_Seurat, idents = c("N", "SN"))
  
  DefaultAssay( PatchSeq_Seurat ) <- "RNA"
  PatchSeq_Seurat[[ "SCT" ]] <- NULL
  PatchSeq_Seurat@assays$RNA@layers$data <- NULL
  PatchSeq_Seurat@assays$RNA@layers$scale.data <- NULL
  
  PatchSeq_Seurat$age <- "PatchSeq"
  PatchSeq_Seurat@reductions$pca <- NULL
  PatchSeq_Seurat@reductions$umap <- NULL
  
  all_neurons_merged <- merge( x = Seurat_neuron_EB_merged, y = c(diBella.obj, PatchSeq_Seurat ) )
  
  all_neurons_merged <- JoinLayers( all_neurons_merged )
  saveRDS( object = all_neurons_merged, file = paste(nobackup_base, "RDS_files/all_neurons_merged.RDS", sep=""))
  
}

all_neurons_merged <- readRDS( file = paste(nobackup_base, "RDS_files/all_neurons_merged.RDS", sep=""))

# prepare a grouping variable that separates the PatchSeq plates - to account for batch effects
new_group <- sapply( 1:nrow(all_neurons_merged@meta.data), function (row) {
  if (all_neurons_merged@meta.data$age[row] == "PatchSeq") {
    return(paste("PatchSeq", all_neurons_merged@meta.data$plate[row], sep="."))
  } else {
    return (all_neurons_merged@meta.data$age[row])
  }
})

all_neurons_merged <- NormalizeData( all_neurons_merged )

all_neurons_merged <- FindVariableFeatures( all_neurons_merged, nfeatures = 3000 )

all_neurons_merged <- ScaleData( all_neurons_merged )
all_neurons_merged <- RunPCA( all_neurons_merged, assay = "RNA", npcs = 30)

all_neurons_merged@meta.data$mergeGroup <- new_group

# integrate the data
all_neurons_merged <- cluster_sim_spectrum( object = all_neurons_merged, 
                                            spectrum_type = "corr_ztransform", 
                                            label_tag = "mergeGroup", 
                                            cluster_col = "seurat_clusters" )

all_neurons_merged <- RunUMAP(all_neurons_merged, reduction = "css", dims = 1:ncol(Embeddings(all_neurons_merged, "css")), 
                         n.neighbors = 30, metric = "euclidean", 
                         min.dist = 0.3, return.model = F)

# UMAP separates cells mainly based on origin - diBella, EB, PatchSeq
DimPlot( all_neurons_merged, group.by="age")
ggsave( filename = paste(base, "/plots/QC/Initial_integrated_umap.png", sep="") )

# separate the different data sets into sub sets for further analysis
Idents ( all_neurons_merged ) <- "age"

# removing non-reference cells is easier than extracting diBella cells
reference_neurons <- subset( all_neurons_merged, idents = c("D13", "D20", "D25", "PatchSeq"), invert=T )
# extract EB cells
EB_neurons <- subset( all_neurons_merged, idents = c("D13", "D20", "D25") )
# extract Patch-Seq cells
PatchSeq_neurons <- subset( all_neurons_merged, idents = c("PatchSeq") )

# create a reference UMAP - use all css dimensions
reference_neurons <- RunUMAP(reference_neurons, reduction = "css", dims = 1:ncol(Embeddings(EB_neurons, "css")), 
                              n.neighbors = 30, metric = "euclidean", 
                              min.dist = 0.3, return.model = T, seed.use = 2401)

reference_neurons <- RunUMAP(reference_neurons, reduction = "css", dims = 1:ncol(Embeddings(EB_neurons, "css")), 
                                n.neighbors = 30, metric = "euclidean", n.components = 3,
                                min.dist = 0.3, return.model = T, seed.use = 2401, , reduction.name = "umap_3d")

# project the EB neurons on that reference - this will be ref.umap
EB_neurons <- ProjectUMAP( query = EB_neurons, query.reduction = "css", query.dims = 1:ncol(Embeddings(reference_neurons, "css")), 
                  reference = reference_neurons, reference.reduction = "css", reference.dims = 1:ncol(Embeddings(reference_neurons, "css")), 
                  reduction.model = "umap")

EB_neurons <- ProjectUMAP( query = EB_neurons, query.reduction = "css", query.dims = 1:ncol(Embeddings(reference_neurons, "css")), 
                           reference = reference_neurons, reference.reduction = "css", reference.dims = 1:ncol(Embeddings(reference_neurons, "css")), 
                           reduction.name = "umap_3d", reduction.model = "umap_3d")

# project the PatchSeq neurons on that reference
PatchSeq_neurons <- ProjectUMAP( query = PatchSeq_neurons, query.reduction = "css", query.dims = 1:ncol(Embeddings(reference_neurons, "css")), 
                           reference = reference_neurons, reference.reduction = "css", reference.dims = 1:ncol(Embeddings(reference_neurons, "css")), 
                           reduction.model = "umap")

PatchSeq_neurons <- ProjectUMAP( query = PatchSeq_neurons, query.reduction = "css", query.dims = 1:ncol(Embeddings(reference_neurons, "css")), 
                                 reference = reference_neurons, reference.reduction = "css", reference.dims = 1:ncol(Embeddings(reference_neurons, "css")), 
                                 reduction.name = "umap_3d", reduction.model = "umap_3d")

# prepare DiBella reference for display 
# define shorter cell type descriptions for display
short_cellType_label <- c( "AP","CThPN","DL_CPN","ImNeurons","IP","Layer_4","Layer_6b","MigrNeurons","NP","SCPN","UL_CPN") 
names( short_cellType_label ) <- sort(unique(reference_neurons$New_cellType))
reference_neurons@meta.data$short_cellType_label <- short_cellType_label[ reference_neurons@meta.data$New_cellType]

# here I use only 150 dimensions - worked best for me
reference_neurons <- FindNeighbors(reference_neurons, dims = 1:150, verbose = FALSE, reduction = "css")
reference_neurons <- FindClusters(reference_neurons, resolution = 1.6, verbose = FALSE)

reference_neurons@meta.data$broad_group <- NA
reference_neurons@meta.data$broad_group[ which(reference_neurons@meta.data$seurat_clusters %in% c("4","0", "11","8")) ] <- "UL"
reference_neurons@meta.data$broad_group[ which(reference_neurons@meta.data$seurat_clusters %in% c("21", "16", "14", "1", "22", "19")) ] <- "DL"
reference_neurons@meta.data$broad_group[ which(reference_neurons@meta.data$seurat_clusters %in% c("9","7","5","10","17", "15")) ] <- "migrating"
reference_neurons@meta.data$broad_group[ which(reference_neurons@meta.data$seurat_clusters %in% c("12","18", "13", "20")) ] <- "IP"
reference_neurons@meta.data$broad_group[ which(reference_neurons@meta.data$seurat_clusters %in% c("5","3","6", "2")) ] <- "AP"

# 3D plots for QC
# plot PatchSeq cells on the reference
umap_patch <- as.data.frame( Embeddings( PatchSeq_neurons, reduction = "umap_3d") )
umap_patch$group <- "patch"
umap_patch$cloneType <- PatchSeq_neurons$clone_type

umap_ref <- as.data.frame( Embeddings( reference_neurons, reduction = "umap_3d") )

umap_ref$broad_group <- reference_neurons$broad_group
umap_ref$short_cellType_label <- reference_neurons$short_cellType_label
umap_ref$seurat_clusters <- reference_neurons$seurat_clusters
umap_ref$short_cellType_label <- reference_neurons$short_cellType_label

# define a color palette for cell annotation
broad_group_col_list <- list()
broad_group_col_list["AP"] <- "grey80"
broad_group_col_list["IP"] <- "grey60"
broad_group_col_list["migrating"] <- "grey40"
broad_group_col_list["DL"] <- "goldenrod"
broad_group_col_list["UL"] <- "cornflowerblue"

# do the plot
plot3d( umap_ref, col = broad_group_col_list[ umap_ref$broad_group ], size = 4, alpha = 0.3, decorate = F)
legend3d("topright", legend = names(broad_group_col_list), col = as.character(broad_group_col_list), pch = 19)

# save the view point if necessary
if (!file.exists(file = paste(base, "/plots/Fig_4/view_matrix.RDS", sep=""))) {
  # PatchSeq cells need to be well visible later - set viewpoint with PatchSeq cells here
  points3d( umap_patch, col = "black", size = 7 )
  um <- par3d()$userMatrix
  saveRDS( object = um, file = paste(base, "/plots/Fig_4/view_matrix.RDS", sep=""))
} else {
  um <- readRDS( file = paste(base, "/plots/Fig_4/view_matrix.RDS", sep=""))
}

# set viewpoint
view3d(userMatrix = um)
# prepare UMAP with only reference cells
snapshot3d( paste(base, "/plots/QC/Fig_4_Reference_umap3d_broad_group.png", sep=""), fmt = "png", width = 4000, height = 4000)

# quick check that Reference annotation fits my broad categories
sapply( c("AP", "IP", "migrating", "DL", "UL"), function (broad_group){
  sapply( unique(umap_ref$short_cellType_label), function (ref_label){
    length(which(umap_ref$broad_group == broad_group & umap_ref$short_cellType_label == ref_label))
  })
}) -> broad_group_association_mat

pheatmap( apply(broad_group_association_mat, 1, function (x) x/sum(x) ), cluster_rows = F, filename = paste(base, "/plots/QC/Fig_4_ref_broad_group.pdf", sep="") )
dev.off()

# set factor levels for broad_group
reference_neurons$broad_group <- factor( reference_neurons$broad_group, levels = c("AP", "IP", "migrating", "DL", "UL"))

# define broad groups in the reference data - 2D for QC
ref_cellType_plot <- DimPlot( reference_neurons, group.by="short_cellType_label", reduction = "umap", label=F)
ref_age_plot <- DimPlot( reference_neurons, group.by="age", reduction = "umap", label=F)
ref_Satb2_plot <- FeaturePlot( reference_neurons, features = c("Satb2"), order=T, min.cutoff = "q10", reduction = "umap") + NoLegend()
ref_Ctip2_plot <- FeaturePlot( reference_neurons, features = c("Bcl11b"), order=T, min.cutoff = "q10", reduction = "umap") + NoLegend()
ref_cluster_plot <- DimPlot( reference_neurons, group.by="seurat_clusters", reduction = "umap", label=T) + NoLegend()
ref_broadGroup_plot <- DimPlot( reference_neurons, group.by="broad_group", reduction = "umap", label=T) + NoLegend()

ref_cellType_plot + ref_age_plot + ref_Satb2_plot + ref_Ctip2_plot + ref_cluster_plot + ref_broadGroup_plot +
  plot_layout(ncol = 2)
ggsave( filename = paste(base, "/plots/QC/Fig_4_Neuron_Reference_overview.pdf", sep=""), width=10, height=12)

#########################################
###
# transfer labels from reference to organoid
###
#########################################

# determine nearest neighbors for label transfer based on 2D UMAP embedding
nearest_neighbor <- nn2(Embeddings(reference_neurons, "umap_3d"), query = Embeddings(EB_neurons, "umap_3d") )

# identify the actual label
nn_cluster <- lapply( 1:nrow(nearest_neighbor$nn.idx), function (i) as.character( reference_neurons@meta.data$seurat_clusters[ nearest_neighbor$nn.idx[i,] ] ) )

# certainty of label is determined by the % of neighboring reference cells with the same label
max_certainty <- sapply( nn_cluster, function (x){
  tmp <- table(x)
  tmp <- tmp/sum(tmp)
  return ( unique( tmp[which(tmp==max(tmp))] ) )
})

# How many cells have a >80% of neighbors with the same label?
length(which(max_certainty >= 0.8)) / length(max_certainty)
# 0.865423

# cells with an 2 or more labels of equal abundance can't be labeled unambiguously
# determine % here
cluster_assign <- sapply( nn_cluster, function (x){
  tmp <- table(x)
  tmp <- tmp/sum(tmp)
  max_tmp <- tmp[which(tmp==max(tmp))]
  # more than one label is equally abundant in th 10 nearest neighbors
  if( length(max_tmp) > 1) {
    return(NA)
  } else {
    return(names(max_tmp))
  }
})

# unassigned cells
length(which(is.na(cluster_assign))) / length(cluster_assign)
# 0.01891812

# add assignment to meta data
EB_neurons@meta.data$nn_clusters_dB <- cluster_assign

# assign cell types as for the reference
EB_neurons@meta.data$broad_group <- NA
EB_neurons@meta.data$broad_group[ which(EB_neurons@meta.data$nn_clusters_dB %in% c("4","0", "11","8")) ] <- "UL"
EB_neurons@meta.data$broad_group[ which(EB_neurons@meta.data$nn_clusters_dB %in% c("21", "16", "14", "1", "22", "19")) ] <- "DL"
EB_neurons@meta.data$broad_group[ which(EB_neurons@meta.data$nn_clusters_dB %in% c("9","7","5","10","17", "15")) ] <- "migrating"
EB_neurons@meta.data$broad_group[ which(EB_neurons@meta.data$nn_clusters_dB %in% c("12","18", "13", "20")) ] <- "IP"
EB_neurons@meta.data$broad_group[ which(EB_neurons@meta.data$nn_clusters_dB %in% c("5","3","6", "2")) ] <- "AP"

# plot the 3D UMAP
umap_EB <- as.data.frame( EB_neurons@reductions$umap_3d@cell.embeddings )

umap_EB$broad_group <- EB_neurons$broad_group
umap_EB <- umap_EB[ which(!is.na(umap_EB$broad_group)),]

plot3d( umap_EB, col = broad_group_col_list[ umap_EB$broad_group ], size = 4, alpha = 0.3, decorate = F)
legend3d("topright", legend = names(broad_group_col_list), col = as.character(broad_group_col_list), pch = 19)
# set viewpoint
view3d(userMatrix = um)
snapshot3d( paste(base, "/plots/QC/Fig_4_EB_umap3d_broad_group.png", sep=""), fmt = "png", width = 4000, height = 4000)

# save EB neuron data for separate analysis
saveRDS( object = EB_neurons@meta.data, file = paste(nobackup_base, "/RDS_files/EB_neurons_meta.RDS", sep="") )

saveRDS( object = EB_neurons, file = paste(nobackup_base, "/RDS_files/EB_neurons.RDS", sep="") )

####################################################################
##
# QC of organoid neurons 
# determine the overlap of marker genes between reference and EBs
##
####################################################################

Idents(EB_neurons) <- "broad_group"
Idents(reference_neurons) <- "broad_group"

EB_layer_markers <- FindMarkers( object = EB_neurons, ident.1 = "UL", ident.2 = "DL")
EB_layer_markers$gene <- rownames(EB_layer_markers)
Ref_layer_markers <- FindMarkers( object = reference_neurons, ident.1 = "UL", ident.2 = "DL")
Ref_layer_markers$gene <- rownames(Ref_layer_markers)

lFC2plot <- EB_layer_markers[ which(EB_layer_markers$p_val_adj < 0.01), c("gene", "avg_log2FC", "p_val_adj")]
colnames(lFC2plot) <- c("gene", "EB_FC", "EB_padj")
length( which(EB_layer_markers$p_val_adj < 0.01) )
# [1] 4459

tmp <- Ref_layer_markers[ which(Ref_layer_markers$p_val_adj < 0.01), c("gene", "avg_log2FC", "p_val_adj")]
colnames(tmp) <- c("gene", "Ref_FC", "Ref_padj")
length( which(Ref_layer_markers$p_val_adj < 0.01) )
# [1] 4683

lFC2plot <- merge(lFC2plot, tmp, by="gene")
nrow(lFC2plot)
# 2106

# prepare simple overlap analysis
out_vec <- c()
for ( i in 1:nrow(lFC2plot)) {
  if ( lFC2plot$EB_FC[i] < 0 & lFC2plot$Ref_FC[i] < 0 ) {
    out_vec <- c(out_vec, "down_down")
  } else if ( lFC2plot$EB_FC[i] > 0 & lFC2plot$Ref_FC[i] > 0 ) {
    out_vec <- c(out_vec, "up_up")
  } else if ( lFC2plot$EB_FC[i] < 0 & lFC2plot$Ref_FC[i] > 0 ) {
    out_vec <- c(out_vec, "down_up")
  }else {
    out_vec <- c(out_vec, "up_down")
  }
}

lFC2plot$group <- out_vec

# prepare correlation plot
marker_dotPlot <- ggplot( lFC2plot, aes(x=EB_FC, y=Ref_FC, color=out_vec)) + geom_point() + theme_classic() + 
  ylab ("logFC UL/DL Ref") + xlab("logFC UL/DL EB") + scale_color_manual( values = c("up_down" = "grey60", "down_up" = "grey60", "up_up"="black", "down_down"="black"))
ggsave( filename = paste(base, "/plots/Sup_Fig_22/S22b_EB_ref_layer_marker_overlap.pdf", sep=""), plot = marker_dotPlot )

# save this list 
# write out DEG analysis
# this file contains the raw DEG output
DEG_out_list <- list()
DEG_out_list[["EB layer markers"]] <- EB_layer_markers
DEG_out_list[["Ref layer markers"]] <- Ref_layer_markers
DEG_out_list[["Shared layer markers"]] <- lFC2plot

wb <- createWorkbook()

for (name in names(DEG_out_list)) {
  sheetName <- name
  
  addWorksheet(wb, sheetName)
  
  writeData( wb, sheetName, DEG_out_list[[ name ]] )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol(DEG_out_list[[ name ]]) )
  setColWidths( wb, sheetName, cols = 1:ncol(DEG_out_list[[ name ]]), widths="auto" )
  
}

saveWorkbook(wb, file = paste(base, "/plots/Sup_Fig_22/EB_ref_layer_marker_overlap.xslx", sep=""), overwrite = T) 

# show number of genes in a table
out_table <- table(out_vec)
data <- matrix( out_table[ c("up_up", "down_up", "up_down", "down_down")], nrow=2 )

#      [,1] [,2]
# [1,]  638  178
# [2,]  212 1078

# significant association of marker overlap
chisq.test(data)

# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  data
# X-squared = 789.27, df = 1, p-value < 2.2e-16

group_counts <- table( lFC2plot$group ) / nrow(lFC2plot)
# down_down + up_up
group_counts[1] + group_counts[4]
# 0.8148148

df2plot <- as.data.frame( table(out_vec) / length(out_vec) )
#     out_vec       Freq
# 1 down_down 0.51187085
# 2   down_up 0.10066477
# 3   up_down 0.08452042
# 4     up_up 0.30294397

#########################################
###
# Pseudotime analysis
###
#########################################

# pseudotime analysis on reference to do QC later

cds <- as.cell_data_set( reference_neurons )
cds <- cluster_cells(cds)

cds <- learn_graph( cds, verbose = F, close_loop = F, use_partition = F )
cds <- order_cells(cds)

traj_plot_ct <- plot_cells(cds, label_groups_by_cluster = FALSE, label_cell_groups = FALSE,
                           label_leaves = T, label_branch_points = FALSE, label_roots = FALSE, 
                           color_cells_by = "pseudotime", group_label_size = 5, 
                           trajectory_graph_segment_size = 3, )

ggsave( plot = traj_plot_ct, filename = paste(base, "/plots/Sup_Fig_23/Neuron_reference_monocle3_trajectory.png", sep=""),
        width = 10, height = 10)

# plot pseudotime distribution of reference cell types
reference_neurons <- AddMetaData( object = reference_neurons, metadata = pseudotime( cds ), col.name = "pseudotime" )
pseudotime_ref_df <- reference_neurons@meta.data[, c("broad_group", "pseudotime")]
pseudotime_ref_df$broad_group <- factor( pseudotime_ref_df$broad_group, levels = c("AP", "IP", "migrating", "UL", "DL"))
ggplot( pseudotime_ref_df, aes(x= broad_group, y=pseudotime)) + geom_boxplot() + theme_classic()

# not used in the paper
ggsave( filename = paste(base, "/plots/Sup_Fig_23/Neuron_ref_pseudotime_broad_clusters.pdf", sep="") )

# pseudotime for EB cells based on nearest neighbor in the Reference data set
EB_neurons@meta.data$pseudotime <- sapply( 1:nrow(nearest_neighbor$nn.idx), function (i) mean( reference_neurons@meta.data$pseudotime[ nearest_neighbor$nn.idx[i,] ] ) )

pseudotime_EB_df <- EB_neurons@meta.data[, c("broad_group", "pseudotime")]
pseudotime_EB_df$broad_group <- factor( pseudotime_EB_df$broad_group, levels = c("AP", "IP", "migrating", "UL", "DL") )
pseudotime_EB_df <- pseudotime_EB_df[ which(!is.na(pseudotime_EB_df$broad_group)), ]
ggplot( pseudotime_EB_df, aes(x=broad_group, y=pseudotime)) + geom_boxplot() + theme_classic()

# not used in the paper
ggsave( filename = paste(base, "/plots/Sup_Fig_23/Neuron_pseudotime_broad_clusters.pdf", sep="") )

pseudotime_ref_df$origin <- "Ref"
pseudotime_EB_df$origin <- "EB"

comb_pt_df <- rbind(pseudotime_ref_df, pseudotime_EB_df)

# calculate distribution of pseudotime in the different cell types
comb_pt_df %>%
  group_by(origin, broad_group) %>%
  summarise(n=n(), mean = mean(pseudotime), sd = sd(pseudotime)) -> pt_summary

pt_summary <- as.data.frame( pt_summary )

# plot this distribution
ggplot( pt_summary, aes(x=origin, y=mean)) + geom_point() + 
  geom_errorbar( aes(ymin=mean-sd, ymax=mean+sd)) + 
  facet_grid(~broad_group) +
  theme_classic()

# used for reviewer's figure
ggsave( filename = paste(base, "/plots/Sup_Fig_23/Neuron_pseudotime_broad_clusters_EB_vs_Ref.pdf", sep="") )

#statistics of difference
pt_stats <- data.frame( EB_mean = pt_summary[which(pt_summary$origin == "EB"), "mean"], 
            Ref_mean = pt_summary[which(pt_summary$origin == "Ref"), "mean"], 
            Ref_sd = pt_summary[which(pt_summary$origin == "Ref"), "sd"] )
rownames(pt_stats) <- pt_summary$broad_group[1:5]

pt_stats$z_score <- (pt_stats$EB_mean - pt_stats$Ref_mean) / pt_stats$Ref_sd
pt_stats$p_value_from_z <- sapply(pt_stats$z_score, function (x) { 
  pvalue <- pnorm(q = x, lower.tail=FALSE) 
  pvalue <- formatC(pvalue, format = "e", digits = 2)
})

# write out table and stats
wb <- createWorkbook()

sheetName <- "pt stats"

addWorksheet(wb, sheetName)
writeData( wb, sheetName, pt_summary )
addFilter( wb, sheetName, row = 1, cols = 1:ncol(pt_summary) )
setColWidths( wb, sheetName, cols = 1:ncol(pt_summary), widths="auto" )

sheetName <- "pt diff stats"

addWorksheet(wb, sheetName)
writeData( wb, sheetName, pt_stats )
addFilter( wb, sheetName, row = 1, cols = 1:ncol(pt_stats) )
setColWidths( wb, sheetName, cols = 1:ncol(pt_stats), widths="auto" )

saveWorkbook(wb, file = paste(base, "/plots/Sup_Fig_23/Pseudotime_EB_vs_Ref.xlsx", sep=""), overwrite = T) 

####################################
# Emx1 expression during pseudotime - reviewer figure
####################################

emx1_pseudotime <- as.data.frame( GetAssayData( reference_neurons, layer = "data")["Emx1",] )
emx1_pseudotime$pseudotime <- cut( reference_neurons@meta.data$pseudotime, breaks = 10, labels = c(1:10))
colnames(emx1_pseudotime) <- c("Emx1", "ps_bin")
emx1_pseudotime$group <- "Reference"

tmp_emx1_pseudotime <- as.data.frame( GetAssayData( EB_neurons, layer = "data")["Emx1",] )
tmp_emx1_pseudotime$pseudotime <- cut( EB_neurons@meta.data$pseudotime, breaks = 10, labels = c(1:10))
colnames(tmp_emx1_pseudotime) <- c("Emx1", "ps_bin")
tmp_emx1_pseudotime$group <- "organoid"

emx1_pseudotime<- rbind( emx1_pseudotime, tmp_emx1_pseudotime)

emx1_pseudotime %>%
  group_by( ps_bin, group ) %>%
  summarize ( mean=mean(Emx1), sd=sd(Emx1)) -> Emx1_ps

ggplot(Emx1_ps, aes(x=ps_bin, y=mean, color=group, group=group)) + 
  geom_ribbon( aes(ymin=mean-sd, ymax=mean+sd, fill = group), alpha = 0.5 ) + 
  geom_line() + geom_point() + ylab("mean normalized Emx1 expression") + xlab("pseudotime bin") +
  theme_classic() + facet_wrap(~group)

ggsave( filename = paste(base, "/plots/reviewer_figure/Emx1_pseudotime.png", sep="") , width = 7)


#########################################
###
# PatchSeq analysis
###
#########################################

# assign PatchSeq cells to reference
# determine closest neighbors in UMAP space
nearest_neighbor <- nn2(Embeddings(reference_neurons, "umap_3d"), query = Embeddings(PatchSeq_neurons, "umap_3d") )

# cell type assignment based on nearest neighbors of the reference
nn_cluster <- lapply( 1:nrow(nearest_neighbor$nn.idx), function (i) as.character( reference_neurons@meta.data$broad_group[ nearest_neighbor$nn.idx[i,] ] ) )

max_certainty <- sapply( nn_cluster, function (x){
  tmp <- table(x)
  tmp <- tmp/sum(tmp)
  return ( unique( tmp[which(tmp==max(tmp))] ) )
})

length(which(max_certainty >= 0.8)) / length(max_certainty)
# 0.8638498

cluster_assign <- sapply( nn_cluster, function (x){
  tmp <- table(x)
  tmp <- tmp/sum(tmp)
  max_tmp <- tmp[which(tmp==max(tmp))]
  if( length(max_tmp) > 1) {
    return( NA )
  } else {
    return(names(max_tmp))
  }
})

# unassigned cells
length(which(is.na(cluster_assign))) / length(cluster_assign)
# 0.02816901
1-length(which(is.na(cluster_assign))) / length(cluster_assign)
# 0.971831

# add assignment to meta data
PatchSeq_neurons@meta.data$nn_clusters_dB <- cluster_assign
PatchSeq_neurons@meta.data$certainty <- max_certainty

# determine age of patched cells based on pseudotime
PatchSeq_neurons@meta.data$pseudotime <- sapply( 1:nrow(nearest_neighbor$nn.idx), function (i) mean( reference_neurons@meta.data$pseudotime[ nearest_neighbor$nn.idx[i,] ] ) )

# clone stats
PatchSeq_neurons@meta.data$nn_clusters_dB[ which( is.na( PatchSeq_neurons@meta.data$nn_clusters_dB ) ) ] <- "undefined"
meta_before_filtering <- PatchSeq_neurons@meta.data

# save this metadata
wb <- createWorkbook()
sheetName <- "MADM Clone-Seq"
addWorksheet(wb, sheetName)
writeData( wb, sheetName, meta_before_filtering )
addFilter( wb, sheetName, row = 1, cols = 1:ncol(meta_before_filtering) )
setColWidths( wb, sheetName, cols = 1:ncol(meta_before_filtering), widths="auto" )
saveWorkbook(wb, file = paste(base, "/plots/Sup_Fig_23/MADMCloneSeq_metaData.xlsx", sep=""), overwrite = T) 

table( meta_before_filtering$nn_clusters_dB )
#  DL        UL undefined 
# 128        79         6 

# remove PatchSeq cells with uncertain label
Idents(PatchSeq_neurons) <- "nn_clusters_dB"
PatchSeq_neurons <- subset( PatchSeq_neurons, idents="undefined", invert=T)

# An object of class Seurat 
# 37426 features across 207 samples within 1 assay 
# Active assay: RNA (37426 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 5 dimensional reductions calculated: pca, css, umap, ref.umap, umap_3d


####################################################
##
# calculate marker DEGs and overlap with reference
##
####################################################

Idents(PatchSeq_neurons) <- "nn_clusters_dB"
PatchSeq_layer_markers <- FindMarkers( object = PatchSeq_neurons, 
                                       ident.1 = "UL", ident.2 = "DL", test.use = "MAST", 
                                       features = rownames(Ref_layer_markers)[1:200], logfc.threshold = 0  )
PatchSeq_layer_markers$gene <- rownames( PatchSeq_layer_markers )

# perform a less conservative p-value adjustment
PatchSeq_layer_markers$p_val_adj <- p.adjust( PatchSeq_layer_markers$p_val )

# plot a few well known layer markers
lFC2plot <- PatchSeq_layer_markers[c("Bcl11b", "Foxp2", "Sox5", "Satb2", "Pou3f3", "Cux2"), c("gene", "avg_log2FC", "p_val_adj"), drop=F]
lFC2plot$layer <- c("DL", "DL", "DL", "UL", "UL", "UL")
lFC2plot$gene <- factor( lFC2plot$gene, levels = lFC2plot$gene[ order(lFC2plot$avg_log2FC)])
lFC2plot$sig <- sapply( lFC2plot$p_val_adj, function(p) ifelse(p<0.1,"*",""))

ggplot( lFC2plot, aes(x=gene, y=avg_log2FC, fill=layer)) + geom_bar(stat="identity") + 
  geom_text( aes(label = sig), size=10) + geom_text(aes(label=sig)) +
  scale_fill_manual(values=c("DL" = "cornflowerblue", "UL" = "goldenrod")) + theme_classic()
ggsave( filename = paste(base, "/plots/Sup_Fig_23/S23c_Neuron_PatchSeq_layerMarkers.pdf", sep="") )

# save this list 
# write out DEG analysis
# this file contains the raw DEG output
wb <- createWorkbook()
sheetName <- "layer markers"
addWorksheet(wb, sheetName)

writeData( wb, sheetName, PatchSeq_layer_markers )
addFilter( wb, sheetName, row = 1, cols = 1:ncol(PatchSeq_layer_markers) )
setColWidths( wb, sheetName, cols = 1:ncol(PatchSeq_layer_markers), widths="auto" )

saveWorkbook(wb, file = paste(base, "/plots/Sup_Fig_23/PatchSeq_layer_markers.xlsx", sep=""), overwrite = T) 

### do some analysis of clone composition
clone_composition <- PatchSeq_neurons@meta.data %>%
  group_by(clone, clone_type, clone_size_patched) %>%
  summarise(types = paste( sort(unique(nn_clusters_dB)), collapse=","), n=n() )

# the goal of the analysis is to identify lineage restriction
# clones with a single informative cell are not informative for this analysis
# remove here
clone_composition <- clone_composition[ which( clone_composition$n > 1), ]

# remove PatchSeq clones with only one cell
Idents(PatchSeq_neurons) <- "clone"
PatchSeq_neurons <- subset( PatchSeq_neurons, idents=clone_composition$clone, invert=F)

# An object of class Seurat 
# 37426 features across 195 samples within 1 assay 
# Active assay: RNA (37426 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 5 dimensional reductions calculated: pca, css, umap, ref.umap, umap_3d

# define colors for clone types
# clone type on reference
col2type <- list()
col2type[[ "N" ]] <- "blue"
col2type[[ "SN" ]] <- "turquoise"

# basic QC showing similar feature distribution for N and SN
VlnPlot( PatchSeq_neurons, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "clone_type") &
  scale_fill_manual( values = col2type)
ggsave( paste(base, "/plots/Sup_Fig_22/S22d_Patch_QC.pdf", sep="") )

# plot PatchSeq cells on the reference - separately for N and SN
# get embeddings for filtered cells
umap_patch <- as.data.frame( Embeddings( PatchSeq_neurons, reduction = "umap_3d") )
umap_patch$group <- "patch"
umap_patch$cloneType <- PatchSeq_neurons$clone_type

# first reference umap with axes information
# used for figure making
plot3d( umap_ref, col = broad_group_col_list[ umap_ref$broad_group ], size = 4, alpha = 0.3, decorate = T)
legend3d("topright", legend = names(broad_group_col_list), col = as.character(broad_group_col_list), pch = 19)
# set viewpoint
view3d(userMatrix = um)
snapshot3d( paste(base, "/plots/Fig_4/Ref_broad_group_axes.png", sep=""), fmt = "png", width = 4000, height = 4000)

rgl.close()

# umap without/with cell type information
plot3d( umap_ref, col = broad_group_col_list[ umap_ref$broad_group ], size = 4, alpha = 0.3, decorate = F)
points3d( umap_patch, col = "black", size = 10 )
# set viewpoint
view3d(userMatrix = um)
snapshot3d( paste(base, "/plots/Fig_4/F4b_Patch_on_Ref_broad_group_all.png", sep=""), fmt = "png", width = 4000, height = 4000)

# not used in figure
plot3d( umap_ref, col = broad_group_col_list[ umap_ref$broad_group ], size = 4, alpha = 0.3, decorate = F)
points3d( umap_patch[which(umap_patch$cloneType == "N"),], col = "black", size = 10 )
# set viewpoint
view3d(userMatrix = um)
snapshot3d( paste(base, "/plots/Fig_4/Patch_on_Ref_broad_group_N.png", sep=""), fmt = "png", width = 4000, height = 4000)

# not used in figure
plot3d( umap_ref, col = broad_group_col_list[ umap_ref$broad_group ], size = 4, alpha = 0.3, decorate = F)
points3d( umap_patch[which(umap_patch$cloneType == "SN"),], col = "black", size = 10 )
# set viewpoint
view3d(userMatrix = um)
snapshot3d( paste(base, "/plots/Fig_4/Patch_on_Ref_broad_group_SN.png", sep=""), fmt = "png", width = 4000, height = 4000)

rgl.close()

# then with pseudotime information
# color palette
colfunc <- colorRampPalette(c("blue", "orange"))

# draw out the legend for figure making
legend_image <- as.raster(matrix(colfunc(100)[1:100], ncol=1))
pdf( file = paste(base, "/plots/Fig_4/Pseudotime_color_legend.pdf", sep=""))
plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '') #Generates a blank plot
text(x=0.4, y = c(0.2,1), labels = c(round(max(reference_neurons@meta.data$pseudotime)), 0), cex = 1.5) #Creates the numeric labels on the scale
rasterImage(legend_image, 0.5, 0.2, 0.7,1) #Values can be modified here to alter where and how wide/tall the gradient is drawn in the plotting area
dev.off()

umap_ref$pseudotime <- cut( reference_neurons@meta.data$pseudotime, breaks = 100, labels = c(1:100))
plot3d( umap_ref, col = colfunc(100)[umap_ref$pseudotime], size = 4, alpha = 0.3, decorate = F)
points3d( umap_patch[which(umap_patch$cloneType == "N"),], col = "black", size = 10 )
# set viewpoint
view3d(userMatrix = um)
snapshot3d( paste(base, "/plots/Fig_4/F4e_Patch_on_Ref_pseudotime_N.png", sep=""), fmt = "png", width = 4000, height = 4000)

plot3d( umap_ref, col = colfunc(100)[umap_ref$pseudotime], size = 4, alpha = 0.3, decorate = F)
points3d( umap_patch[which(umap_patch$cloneType == "SN"),], col = "black", size = 10 )
# set viewpoint
view3d(userMatrix = um)
snapshot3d( paste(base, "/plots/Fig_4/F4f_Patch_on_Ref_pseudotime_SN.png", sep=""), fmt = "png", width = 4000, height = 4000)

# plot the clone types
plot3d( umap_ref, col = "grey80", size = 4, alpha = 0.3, decorate = F)
points3d( umap_patch[which(umap_patch$cloneType == "N"),], col = col2type[ umap_patch[which(umap_patch$cloneType == "N"), "cloneType" ]],  size = 10 )
# set viewpoint
view3d(userMatrix = um)
legend3d("topright", legend = names(col2type), col = as.character(col2type), pch = 19)
snapshot3d( paste(base, "/plots/Fig_4/F4c_Patch_on_Ref_cloneType_N.png", sep=""), fmt = "png", width = 4000, height = 4000)

plot3d( umap_ref, col = "grey80", size = 4, alpha = 0.3, decorate = F)
points3d( umap_patch[which(umap_patch$cloneType == "SN"),], col = "black", size = 12, alpha = 0.8, )
points3d( umap_patch[which(umap_patch$cloneType == "SN"),], col = col2type[ umap_patch[which(umap_patch$cloneType == "SN"), "cloneType" ]], size = 10)
# set viewpoint
view3d(userMatrix = um)
snapshot3d( paste(base, "/plots/Fig_4/F4d_Patch_on_Ref_cloneType_SN.png", sep=""), fmt = "png", width = 4000, height = 4000)

############################################
###
##
# distance within and between clone types
##
###
############################################

nn_intra <- nn2(umap_patch[which(umap_patch$cloneType == "SN"), c(1:3)], query = umap_patch[which(umap_patch$cloneType == "SN"), c(1:3)] )
nn_inter <- nn2(umap_patch[which(umap_patch$cloneType == "N"), c(1:3)], query = umap_patch[which(umap_patch$cloneType == "SN"), c(1:3)] )
dist2plot <- data.frame( dist = nn_intra$nn.dists[,2], cloneType = "SN", group = "intra")
dist2plot <- rbind(dist2plot, data.frame( dist = nn_inter$nn.dists[,1], cloneType = "SN", group = "inter") )

nn_intra <- nn2(umap_patch[which(umap_patch$cloneType == "N"), c(1:3)], query = umap_patch[which(umap_patch$cloneType == "N"), c(1:3)] )
nn_inter <- nn2(umap_patch[which(umap_patch$cloneType == "SN"), c(1:3)], query = umap_patch[which(umap_patch$cloneType == "N"), c(1:3)] )
dist2plot <- rbind(dist2plot, data.frame( dist = nn_intra$nn.dists[,2], cloneType = "N", group = "intra") )
dist2plot <- rbind(dist2plot, data.frame( dist = nn_inter$nn.dists[,1], cloneType = "N", group = "inter") )

ggplot( dist2plot, aes(x=group, y=dist)) + geom_boxplot() + theme_classic() + facet_wrap(~cloneType) + 
  geom_hline( yintercept = c(0.243, 0.199) )
ggsave( filename = paste(base, "/plots/Sup_Fig_23/S23b_Neuron_clone_umap_dist.pdf", sep="") )

# write out stats for the paper

dist2plot %>%
  group_by( clone_type, group) %>%
  summarize( mean = mean(dist),
             median = median(dist),
             sd = sd(dist),
             n=n()) -> Panel_b_stats.csv
write.csv (x = as.data.frame(Panel_b_stats.csv), file = paste(base, "plots/Sup_Fig_23/Panel_b_stats.csv", sep=""), row.names = T)

test_aov <- aov(dist2plot$dist ~ factor(dist2plot$cloneType) * factor(dist2plot$group))
summary(test_aov)

#                                                       Df Sum Sq Mean Sq F value Pr(>F)  
# factor(dist2plot$cloneType)                           1   0.50  0.5004   4.719 0.0304 *
# factor(dist2plot$group)                               1   0.11  0.1102   1.039 0.3087  
# factor(dist2plot$cloneType):factor(dist2plot$group)   1   0.43  0.4339   4.092 0.0438 *
# Residuals                                           386  40.93  0.1060                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

out <- TukeyHSD ( test_aov )

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = dist2plot$dist ~ factor(dist2plot$cloneType) * factor(dist2plot$group))
# 
# $`factor(dist2plot$cloneType)`
#             diff        lwr          upr     p adj
# SN-N -0.07214162 -0.1374373 -0.006845978 0.0304434
# 
# $`factor(dist2plot$group)`
#        diff         lwr       upr     p adj
# intra-interf -0.03361425 -0.09845411 0.0312256 0.3087086
# 
# $`factor(dist2plot$cloneType):factor(dist2plot$group)`
#                           diff         lwr         upr     p adj
# SN:inter-N:inter  -0.139318158 -0.26050353 -0.01813279 0.0168202
# N:intra-N:inter   -0.092867408 -0.20668178  0.02094696 0.1531849
# SN:intra-N:inter  -0.097832485 -0.21901785  0.02335288 0.1604673
# N:intra-SN:inter   0.046450750 -0.07473462  0.16763612 0.7557948
# SN:intra-SN:inter  0.041485672 -0.08664737  0.16961871 0.8375826
# SN:intra-N:intra  -0.004965077 -0.12615045  0.11622029 0.9995766

write.csv (x = as.data.frame(test_aov), file = paste(base, "plots/Sup_Fig_23/Panel_b_ANOVA_stats.csv", sep=""), row.names = T)
write.csv (x = as.data.frame(out[[3]]), file = paste(base, "plots/Sup_Fig_23/Panel_b_ANOVA_pw_stats.csv", sep=""), row.names = T)


# summary stats of clone data
PatchSeq_neurons@meta.data %>%
  group_by(clone, clone_type, clone_size_complete) %>%
  summarise(n=n()) -> clone_summary

clone_summary$clone_size_complete <- as.numeric(clone_summary$clone_size_complete)

clone_summary$coverage <- clone_summary$n / clone_summary$clone_size_complete

clone_summary %>%
  group_by(clone_type) %>%
  summarize (mean = mean(coverage)) -> clone_summary_df

# A tibble: 2 × 2
#   clone_type  mean
#   <chr>      <dbl>
#   1 N          0.672
#   2 SN         0.734

ggplot( ) + 
  ylim(0,1.5) + 
  geom_bar( data = clone_summary_df, aes(x=clone_type, y=mean), stat="identity", position="dodge") + 
  geom_beeswarm( data = clone_summary, aes(x=clone_type, y=coverage)) + 
  theme_classic() + ggtitle("clone coverage") + ylim(0,1)
ggsave( filename = paste(base, "/plots/Sup_Fig_23/S23a_Clone_coverage.pdf", sep=""), width = 2 )

# cell type in each clone type
clone_layer_df <- as.data.frame( table( PatchSeq_neurons@meta.data$nn_clusters_dB[ which( PatchSeq_neurons@meta.data$clone_type == "SN" ) ] ) )
clone_layer_df$frac <- clone_layer_df$Freq / sum( clone_layer_df$Freq )
clone_layer_df$group <- "SN"

tmp_clone_layer_df <- as.data.frame( table( PatchSeq_neurons@meta.data$nn_clusters_dB[ which( PatchSeq_neurons@meta.data$clone_type == "N" ) ] ) )
tmp_clone_layer_df$frac <- tmp_clone_layer_df$Freq / sum( tmp_clone_layer_df$Freq )
tmp_clone_layer_df$group <- "N"

# not used in figure
clone_layer_df <- rbind( clone_layer_df, tmp_clone_layer_df )
ggplot( clone_layer_df, aes(x=Var1, y=frac)) + 
  geom_bar( stat = "identity", position = "dodge", color="black") + 
  # scale_fill_manual( values = col2type ) +
  facet_wrap(~group) + theme_classic()
ggsave( filename = paste(base, "/plots/Sup_Fig_22/Clone_layer_association.pdf", sep=""), width = 2 )

# count clone number
df2plot <- as.data.frame( table( clone_summary$clone_type ) )
ggplot( df2plot, aes(x=Var1, y=Freq)) + geom_bar(stat="identity", position="dodge") + theme_classic()
ggsave( filename = paste(base, "/plots/Sup_Fig_22/Clone_number plot.pdf", sep=""), width = 2 )
#   Var1 Freq
# 1    N   21
# 2   SN   34

# prepare clone type groupings
clone_composition$simple_type <- "ML"
clone_composition$simple_type[ which(clone_composition$types == "DL")] <- "DL"
clone_composition$simple_type[ which(clone_composition$types == "UL")] <- "UL"
clone_composition$simple_type_restriction <- paste(clone_composition$clone_type, clone_composition$simple_type, sep="_")

# add this grouping to Seurat object - note that all single cell clones do not get a label
PatchSeq_neurons@meta.data$restriction_type <- sapply( PatchSeq_neurons$clone, function (clone) {
  if ( length( which( clone_composition$clone == clone ) ) > 0 ) {
    return ( clone_composition$simple_type[ which(clone_composition$clone == clone) ] )
  } else {
    return ( NA )
  }
} )

PatchSeq_neurons@meta.data$clone_restriction_type <- paste(PatchSeq_neurons@meta.data$clone_type, PatchSeq_neurons@meta.data$restriction_type, sep="_")

clone_PS_plot1 <- ggplot( PatchSeq_neurons@meta.data, aes(x=clone_type, y=pseudotime)) + 
  geom_boxplot() + theme_classic() + ylim(0, max(PatchSeq_neurons@meta.data$pseudotime))

clone_PS_plot2 <- ggplot( PatchSeq_neurons@meta.data, aes( x = restriction_type, y = pseudotime )) + 
  geom_boxplot( outlier.size = 0 ) + geom_beeswarm( size = 0.5 ) + facet_wrap( ~ clone_type) + theme_classic() + 
  ylim(0, max(PatchSeq_neurons@meta.data$pseudotime))

clone_PS_plot1 + clone_PS_plot2
ggsave( filename = paste(base, "/plots/Fig_4/F4gh_Pseudotime_lineageRestr_clones.pdf", sep="") )

###################################################
##
# save the Patch-Seq meta data used for plotting
##
###################################################

Patch_metadata_2_save <- PatchSeq_neurons@meta.data
Patch_metadata_2_save <- cbind( Patch_metadata_2_save, umap_patch[ rownames(Patch_metadata_2_save), 1:3] )
Ref_metadata_2_save <- reference_neurons@meta.data
Ref_metadata_2_save <- cbind( Ref_metadata_2_save, umap_ref[ rownames(Ref_metadata_2_save), 1:3] )

wb <- createWorkbook()

sheetName <- "PatchSeq metadata"
addWorksheet(wb, sheetName)
writeData( wb, sheetName, Patch_metadata_2_save )
addFilter( wb, sheetName, row = 1, cols = 1:ncol(Patch_metadata_2_save) )
setColWidths( wb, sheetName, cols = 1:ncol(Patch_metadata_2_save), widths="auto" )

sheetName <- "Reference metadata"
addWorksheet(wb, sheetName)
writeData( wb, sheetName, Ref_metadata_2_save )
addFilter( wb, sheetName, row = 1, cols = 1:ncol(Ref_metadata_2_save) )
setColWidths( wb, sheetName, cols = 1:ncol(Ref_metadata_2_save), widths="auto" )

saveWorkbook(wb, file = paste(base, "/plots/Fig_4/PatchSeq_metadata.xlsx", sep=""), overwrite = T) 

#########################
###
# Pseudotime
###
#########################

PatchSeq_neurons@meta.data %>%
  group_by( clone_type ) %>%
  summarize( mean = mean(pseudotime),
             median = median(pseudotime),
             sd = sd(pseudotime))
# # A tibble: 2 × 4
#   clone_type  mean median    sd
#   <chr>      <dbl>  <dbl> <dbl>
# 1 N           56.2   58.8  5.78
# 2 SN          58.2   58.9  3.91

PatchSeq_neurons@meta.data %>%
  group_by( clone_type, restriction_type) %>%
  summarize( mean = mean(pseudotime),
             median = median(pseudotime),
             sd = sd(pseudotime))

# `summarise()` has grouped output by 'clone_type'. You can override using the `.groups` argument.
# # A tibble: 5 × 5
# # Groups:   clone_type [2]
#   clone_type restriction_type  mean median    sd
#   <chr>      <chr>            <dbl>  <dbl> <dbl>
# 1 N          DL                55.0   59.4  7.98
# 2 N          ML                56.5   58.8  5.17
# 3 SN         DL                58.4   59.3  5.36
# 4 SN         ML                58.2   58.9  3.48
# 5 SN         UL                57.7   58.6  1.59


#========================================
# tests / stats reported in the paper:
#========================================

# only compare N and SN
t.test(PatchSeq_neurons@meta.data$pseudotime[which(PatchSeq_neurons@meta.data$clone_type == "N")], 
       PatchSeq_neurons@meta.data$pseudotime[which(PatchSeq_neurons@meta.data$clone_type == "SN")])

# Welch Two Sample t-test
# 
# data:  clone_data$pseudotime[which(clone_data$clone_type == "N")] and clone_data$pseudotime[which(clone_data$clone_type == "SN")]
# t = -2.8704, df = 188.93, p-value = 0.004567
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -3.3688052 -0.6245023
# sample estimates:
#   mean of x mean of y 
# 56.18996  58.18661 

###########################
##
# restriction_type stats
##
###########################

PatchSeq_neurons@meta.data %>%
  group_by( clone_type, restriction_type ) %>%
  summarize( mean = mean(pseudotime),
             median = median(pseudotime),
             sd = sd(pseudotime)) -> Panel_h_stats
write.csv (x = as.data.frame(Panel_h_stats), file = paste(base, "plots/Fig_4/Panel_h_stats.csv", sep=""), row.names = T)


SN_data <- PatchSeq_neurons@meta.data[which(PatchSeq_neurons@meta.data$clone_type == "SN"),]
test_aov <- aov(SN_data$pseudotime ~ factor(SN_data$restriction_type))
out <- summary(test_aov)
out <- as.data.frame( out[[1]] )
write.csv (x = out, file = paste(base, "plots/Fig_4/Panel_h_ANOVA_stats.csv", sep=""), row.names = T)

# plot clone size and coverage for each restriction type
# SN clones are smaller by definition

# not used for figure
ggplot( clone_composition, aes(x=simple_type_restriction, y=n)) + 
  geom_boxplot() + geom_beeswarm() + 
  theme_classic()+ ylim(0,max(clone_composition$n)+1)
ggsave( filename = paste(base, "/plots/Sup_Fig_23/Clone_size_lineageRestr.pdf", sep="") )

# restriction type histogram
clone_summary <- clone_composition %>%
  group_by(clone_type, simple_type) %>%
  summarise( n=n() )

tot_cloneType <- table( clone_composition$clone_type )
clone_summary$total <- tot_cloneType[ as.character(clone_summary$clone_type) ]
clone_summary$rel <- clone_summary$n / clone_summary$total

# confidence intervals are calculated here
conf_int <- t(sapply (1:nrow(clone_summary), function (x) {
  tmp <- clopper.pearson.ci(clone_summary$n[x], clone_summary$total[x], alpha = 0.05, CI = "two.sided")
}))

conf_int <- as.data.frame( conf_int )

clone_summary <- cbind( clone_summary, conf_int )
clone_summary$Lower.limit <- as.numeric( clone_summary$Lower.limit )
clone_summary$Upper.limit <- as.numeric( clone_summary$Upper.limit )

clone_summary <- as.data.frame(clone_summary)

ggplot( clone_summary, aes(x=simple_type, y=rel)) + 
  geom_bar(stat="identity", position="dodge") + 
  geom_errorbar( aes(ymin = Lower.limit, ymax = Upper.limit)) +
  facet_grid(~clone_type) + theme_classic()
ggsave( filename = paste(base, "/plots/Sup_Fig_23/S23d_Group_count_lineageRestr.pdf", sep="") )

wb <- createWorkbook()

sheetName <- "lineage restricted clones"
addWorksheet(wb, sheetName)
writeData( wb, sheetName, clone_summary )
addFilter( wb, sheetName, row = 1, cols = 1:ncol(clone_summary) )
setColWidths( wb, sheetName, cols = 1:ncol(clone_summary), widths="auto" )

saveWorkbook(wb, file = paste(base, "/plots/Sup_Fig_23/PatchSeq_lineage_restr.xlsx", sep=""), overwrite = T) 

sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.5 LTS
# 
# Matrix products: default
# BLAS:   /opt/R/4.3.2/lib/R/lib/libRblas.so 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=de_AT.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] shiny_1.8.0                 pheatmap_1.0.12             rgl_1.3.1                   openxlsx_4.2.5.2           
# [5] ggbeeswarm_0.7.2            dplyr_1.1.4                 GenBinomApps_1.2.1          monocle3_1.3.7             
# [9] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0 GenomicRanges_1.54.1        GenomeInfoDb_1.38.8        
# [13] IRanges_2.36.0              S4Vectors_0.40.2            MatrixGenerics_1.14.0       matrixStats_1.2.0          
# [17] Biobase_2.62.0              BiocGenerics_0.48.1         SeuratWrappers_0.3.5        patchwork_1.2.0            
# [21] RANN_2.6.1                  simspec_0.0.0.9000          ggplot2_3.5.0               Seurat_5.0.1               
# [25] SeuratObject_5.0.1          sp_2.1-3                   
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.3.2           later_1.3.2             bitops_1.0-7            tibble_3.2.1           
# [6] R.oo_1.26.0             polyclip_1.10-6         fastDummies_1.7.3       lifecycle_1.0.4         globals_0.16.2         
# [11] lattice_0.21-9          MASS_7.3-60             MAST_1.28.0             magrittr_2.0.3          sass_0.4.8             
# [16] limma_3.58.1            plotly_4.10.4           jquerylib_0.1.4         remotes_2.5.0           httpuv_1.6.14          
# [21] sctransform_0.4.1       spam_2.10-0             zip_2.3.1               spatstat.sparse_3.0-3   reticulate_1.44.1      
# [26] cowplot_1.1.3           pbapply_1.7-2           minqa_1.2.6             RColorBrewer_1.1-3      abind_1.4-5            
# [31] zlibbioc_1.48.2         Rtsne_0.17              presto_1.0.0            purrr_1.0.2             R.utils_2.12.3         
# [36] RCurl_1.98-1.14         GenomeInfoDbData_1.2.11 ggrepel_0.9.5           irlba_2.3.5.1           listenv_0.9.1          
# [41] spatstat.utils_3.1-0    goftest_1.2-3           RSpectra_0.16-1         spatstat.random_3.2-2   fitdistrplus_1.1-11    
# [46] parallelly_1.37.0       leiden_0.4.3.1          codetools_0.2-19        DelayedArray_0.28.0     tidyselect_1.2.1       
# [51] farver_2.1.1            viridis_0.6.5           lme4_1.1-35.3           base64enc_0.1-3         spatstat.explore_3.2-6 
# [56] jsonlite_1.8.8          ellipsis_0.3.2          progressr_0.14.0        ggridges_0.5.6          survival_3.5-7         
# [61] systemfonts_1.0.6       progress_1.2.3          tools_4.3.2             ragg_1.3.0              ica_1.0-3              
# [66] Rcpp_1.0.12             glue_1.7.0              gridExtra_2.3           SparseArray_1.2.4       xfun_0.42              
# [71] withr_3.0.0             BiocManager_1.30.22     fastmap_1.1.1           boot_1.3-28.1           fansi_1.0.6            
# [76] digest_0.6.34           rsvd_1.0.5              R6_2.5.1                mime_0.12               textshaping_0.3.7      
# [81] colorspace_2.1-0        scattermore_1.2         tensor_1.5              spatstat.data_3.0-4     R.methodsS3_1.8.2      
# [86] utf8_1.2.4              tidyr_1.3.1             generics_0.1.3          data.table_1.15.0       prettyunits_1.2.0      
# [91] httr_1.4.7              htmlwidgets_1.6.4       S4Arrays_1.2.1          uwot_0.1.16             pkgconfig_2.0.3        
# [96] gtable_0.3.4            lmtest_0.9-40           XVector_0.42.0          htmltools_0.5.7         dotCall64_1.1-1        
# [101] scales_1.3.0            png_0.1-8               knitr_1.45              rstudioapi_0.16.0       reshape2_1.4.4         
# [106] nlme_3.1-163            nloptr_2.0.3            cachem_1.0.8            proxy_0.4-27            zoo_1.8-12             
# [111] stringr_1.5.1           KernSmooth_2.23-22      parallel_4.3.2          miniUI_0.1.1.1          vipor_0.4.7            
# [116] ggrastr_1.0.2           pillar_1.9.0            grid_4.3.2              vctrs_0.6.5             promises_1.2.1         
# [121] xtable_1.8-4            cluster_2.1.4           beeswarm_0.4.0          packrat_0.9.2           cli_3.6.2              
# [126] compiler_4.3.2          rlang_1.1.3             crayon_1.5.2            leidenbase_0.1.27       future.apply_1.11.1    
# [131] labeling_0.4.3          plyr_1.8.9              stringi_1.8.3           viridisLite_0.4.2       deldir_2.0-2           
# [136] assertthat_0.2.1        munsell_0.5.0           lazyeval_0.2.2          spatstat.geom_3.2-8     Matrix_1.6-5           
# [141] RcppHNSW_0.6.0          hms_1.1.3               future_1.33.1           statmod_1.5.0           ROCR_1.0-11            
# [146] memoise_2.0.1           igraph_2.1.4            bslib_0.6.1
