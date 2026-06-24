# ========================================================
# step 07
# script to combine and compare embryo and organoid data
#=========================================================

library (Seurat)
library (uwot)
library (ggplot2)
library (simspec)
library (GeneOverlap)
library(openxlsx)

set.seed(2401)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

# prepare a merged Seurat object to integrate 
if (!file.exists( paste(nobackup_base, "RDS_files/all_seurat_merged.RDS", sep="") )) {
  
  # read embryo data, add an age tag
  seurat_obj_list <- readRDS( file = paste( nobackup_base, "RDS_files/embryo_GFP_perAge.RDS", sep = "") )
  for (i in 1:length(seurat_obj_list)) {
    seurat_obj_list[[ i ]]@meta.data$age <- names( seurat_obj_list )[ i ]
  }
  
  seurat_merged <- merge( x = seurat_obj_list[[1]], y = seurat_obj_list[2:4]  )

  seurat_merged <- subset(seurat_merged, subset = nCount_RNA < 20000 )
  
  # remove unnecessary layers to make object smaller
  seurat_merged@assays$RNA@layers$data.E10.1 <- NULL
  seurat_merged@assays$RNA@layers$scale.data.1 <- NULL
  seurat_merged@assays$RNA@layers$data.1 <- NULL
  seurat_merged@assays$RNA@layers$data.E13.2 <- NULL
  seurat_merged@assays$RNA@layers$scale.data.2 <- NULL
  seurat_merged@assays$RNA@layers$data.2 <- NULL
  seurat_merged@assays$RNA@layers$data.E16.3 <- NULL
  seurat_merged@assays$RNA@layers$scale.data.3 <- NULL
  seurat_merged@assays$RNA@layers$data.3 <- NULL
  seurat_merged@assays$RNA@layers$data.P0.4 <- NULL
  seurat_merged@assays$RNA@layers$scale.data.4 <- NULL
  seurat_merged@assays$RNA@layers$data.4 <- NULL
  
  
  seurat_merged <- NormalizeData( seurat_merged )
  seurat_merged <- FindVariableFeatures( seurat_merged, nfeatures = 5000 )
  
  seurat_merged <- ScaleData( seurat_merged )
  seurat_merged <- RunPCA( seurat_merged, assay = "RNA", npcs = 30)
  
  seurat_merged <- FindNeighbors(seurat_merged, reduction = "pca", dims = 1:15 )
  seurat_merged <- FindClusters(seurat_merged, resolution = 0.25)
  
  seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:15, 
                           n.neighbors = 30, metric = "cosine", 
                           min.dist = 0.3, return.model = T, seed.use = 2401)
  
  # sanity check - not saved
  DimPlot( seurat_merged, group.by = "age" )
  DimPlot( seurat_merged, group.by = "seurat_clusters", label = T )
  
  FeaturePlot( object = seurat_merged, features = c("Top2a", "Eomes", "Bcl11b", 
                                                    "Satb2", "Nhlh2", "Dlx1", 
                                                    "Aldh1l1", "Olig2", "Gfap"), 
               order = T, min.cutoff = "q15", reduction = "umap" ) & NoLegend()
  
  # remove unnecessary data to reduce object size
  
  seurat_merged@assays$RNA@layers$data.E10.1 <- NULL
  seurat_merged@assays$RNA@layers$data.E13.2 <- NULL
  seurat_merged@assays$RNA@layers$data.E16.3 <- NULL
  seurat_merged@assays$RNA@layers$data.P0.4 <- NULL
  seurat_merged@assays$RNA@layers$scale.data <- NULL
  
  seurat_merged@reductions$pca <- NULL
  seurat_merged@reductions$umap <- NULL

  # read the organoid data
  EB.obj <- readRDS( file = paste( nobackup_base, "RDS_files/seurat_obj.stress.CC_correct.RDS", sep = "") )
  
  # remove unnecessary data to reduce object size
  
  EB.obj@assays$RNA@layers$data.D8 <- NULL
  EB.obj@assays$RNA@layers$data.D13 <- NULL
  EB.obj@assays$RNA@layers$data.D20 <- NULL
  EB.obj@assays$RNA@layers$data.D25 <- NULL
  EB.obj@reductions$pca <- NULL
  EB.obj@reductions$umap <- NULL
  
  EB.obj@meta.data$group <- "EB"
  seurat_merged@meta.data$group <- "embryo"
  
  seurat_merged <- merge( seurat_merged, EB.obj )
  
  # save the object on disk
  saveRDS ( object = seurat_merged, file = paste(nobackup_base, "RDS_files/all_seurat_merged.RDS", sep=""))
  
}

# restart to free memory

if (!file.exists(paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.downsample.RDS", sep=""))) {
  
  # read the merged seurat object
  seurat_merged <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.RDS", sep=""))
  
  seurat_merged <- JoinLayers( seurat_merged )
  seurat_merged <- NormalizeData( seurat_merged )
  seurat_merged <- FindVariableFeatures( seurat_merged, nfeatures = 5000 )
  seurat_merged <- ScaleData( seurat_merged )
  seurat_merged <- RunPCA( seurat_merged, assay = "RNA", npcs = 30)
  
  ElbowPlot( seurat_merged, ndims = 30)
  
  # perform integration and UMAP based on optimized parameters
  dims <- 25
  n_neighbors <- 10
  min_dist <- 0.2
  metric <- "cosine"
  
  seurat_merged <- cluster_sim_spectrum( object = seurat_merged, label_tag = "group", 
                                         dims_use = 1:dims,
                                         cluster_col = "seurat_clusters", 
                                         use_scale = F, var_genes = VariableFeatures( seurat_merged ) )
  
  seurat_merged <- RunUMAP(seurat_merged, reduction = "css", dims = 1:ncol(Embeddings(seurat_merged, "css")), 
                           n.neighbors = n_neighbors, metric = metric, n.components = 2,
                           min.dist = min_dist, return.model = F, seed.use = 2401)
  
  seurat_merged <- FindNeighbors(seurat_merged, reduction = "css", dims = 1:ncol(Embeddings(seurat_merged, "css") ) )
  seurat_merged <- FindClusters(seurat_merged, resolution = 1)
  
  DimPlot( object = seurat_merged, group.by = "seurat_clusters", label = T, label.size = 6) + NoLegend()
  ggsave( paste(base, "plots/Sup_Fig_7/UMAP_overview.png", sep=""), width = 10, height = 8)
  
  table( seurat_merged$group )
  #     EB embryo 
  # 103252  24942 
  
  # downsample to get aprox. same number of EB and embryo cells 
  # for comparable visualisation and DEG analysis
  seurat_merged <- SketchData(
    object = seurat_merged,
    ncells = 35000,
    method = "LeverageScore",
    sketched.assay = "sketch",
    seed = 2401
  )
  
  # prepare a downsampled data set with comparable numbers of EB and embryo cells
  # get cell IDs from sketched EB cells
  sketch_cells <- colnames( seurat_merged@assays$sketch )
  EB_cells <- rownames( seurat_merged@meta.data[ which( seurat_merged@meta.data$group == "EB"), ] )
  EB_cells <- intersect(sketch_cells, EB_cells)
  length( EB_cells)
  # 28568
  
  # extract sketched EB cells
  EB_sketch <- subset( seurat_merged, cells = EB_cells)
  # NOTE: later subsetting will remove the UMAP embeddings - save here and add later
  umap_coord <- Embeddings( seurat_merged, reduction = "umap")
  # subset embryo cells and add sketched EB cells
  Idents( seurat_merged ) <- "group"
  embryo_sketch <- subset( seurat_merged, idents = "embryo")
  
  # prepare a merged dataset with comparable numbers of EB and embryo cells
  seurat_sketch <- merge( embryo_sketch, EB_sketch)
  # add back UMAP embedding
  seurat_sketch[["umap"]] <- CreateDimReducObject(embeddings = umap_coord[ colnames(seurat_sketch), ], key = "umap_", assay = "RNA")
  
  seurat_sketch@meta.data$age <- factor( seurat_sketch@meta.data$age, levels = c("D8", "D13", "D20", "D25", "E10", "E13", "E16", "P0"))
  
  # define broad cell types based on marker expression
  # sanity check for cell type assignment - not saved
  DimPlot( object = seurat_sketch, group.by = "seurat_clusters", label = T, label.size = 6) + NoLegend()
  
  cluster2cellType <- c( rep("iN", 9), rep("IP", 6), rep("RGP", 5), rep("aNSC", 7), 
                         "OBNB", "oligo", "oligo", "astro", "astro","CR" )
  names( cluster2cellType ) <- as.character( c(17,10,1,0,7,20,27,15,11,
                                               14,19,18,23,29,8,
                                               21,9,13,2,3,
                                               12,30,22,24,28,25,16,
                                               4, 5, 31, 6, 32, 26) )
  
  
  seurat_sketch@meta.data$cellType <- cluster2cellType[ as.character(seurat_sketch@meta.data$seurat_clusters) ]
  
  # investigate cells in D8 sample that slipped filtering
  D8_meta_data <- readRDS( file = paste( nobackup_base,  "/RDS_files/MetaData.GfpDetect.D8.RDS", sep="" ) )
  cells2remove <- rownames(D8_meta_data) [ which( D8_meta_data$GFP_cluster == "NO" ) ]
  
  seurat_sketch@meta.data$toremove <- "No"
  seurat_sketch@meta.data$toremove[ which( rownames(seurat_sketch@meta.data) %in% paste(cells2remove, "_1", sep="" ) ) ] <- "Yes"
  
  table( seurat_sketch@meta.data$toremove )
  
  #    No   Yes 
  # 53444    65
  
  DimPlot( seurat_sketch, group.by = "toremove" )
  ggsave( paste(base, "plots/QC/D8_65cells_removed.png", sep=""), width = 4, height = 3)
  

  # save the objects
  seurat_merged@meta.data$cellType <- cluster2cellType[ as.character(seurat_merged@meta.data$seurat_clusters) ]
  saveRDS ( object = seurat_sketch@meta.data, file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.meta_data.RDS", sep=""))
  saveRDS ( object = seurat_merged, file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.RDS", sep=""))
  saveRDS ( object = seurat_sketch, file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.downsample.RDS", sep=""))
}


# read the down sampled data set
seurat_sketch <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.downsample.RDS", sep=""))

# remove GFP negative cells missed in first round of removal
Idents( seurat_sketch ) <- "toremove"
seurat_sketch <- subset( seurat_sketch, idents="No" )

# plot some UMAPS for figures
map_plot <- DimPlot( object = seurat_sketch, group.by = "cellType", label = T) + NoLegend()
ggsave( paste(base, "plots/Fig_1/F1g_sketch_UMAP_cellTypes.png", sep=""), width = 4, height = 3, plot = map_plot)

umap_plot <- DimPlot( object = seurat_sketch, group.by = "cellType", split.by="age", 
                      label = T, ncol = 4) + NoLegend()
ggsave( paste(base, "plots/Fig_1/F1j_sketch_UMAP_cellTypes_perAge.png", sep=""), width = 8, height = 6, plot = umap_plot)

# define broad cell types based on marker expression
DimPlot( object = seurat_sketch, group.by = "seurat_clusters", label = T, label.size = 6) + NoLegend()
ggsave( paste(base, "plots/Sup_Fig_7/S7a_sketch_UMAP_overview.png", sep=""), width = 8, height = 6)

# further analysis of the cell types between EB and embryo

# find similarity of cluster markers
# test how similar the molecular features of the cell types in EB and embryo are

EB_cells <- rownames( seurat_sketch@meta.data [ seurat_sketch@meta.data$group == "EB", ] )
Embryo_cells <- rownames( seurat_sketch@meta.data [ seurat_sketch@meta.data$group == "embryo", ] )

Idents(seurat_sketch) <- "cellType"
seurat_sketch <- JoinLayers( seurat_sketch )

seurat_sketch_EB <- subset( seurat_sketch, cells = EB_cells )
seurat_sketch_embryo <- subset( seurat_sketch, cells = Embryo_cells )

cellType_markers_EB <- FindAllMarkers( seurat_sketch_EB, only.pos = T )
cellType_markers_Embryo <- FindAllMarkers( seurat_sketch_embryo, only.pos = T )

sapply ( unique(cellType_markers_EB$cluster), function (EB_cluster) {
  sapply ( unique(cellType_markers_Embryo$cluster), function (Embryo_cluster) {
    
    EB_markers <- cellType_markers_EB$gene[ which( cellType_markers_EB$cluster == EB_cluster )][1:100]
    Embryo_markers <- cellType_markers_Embryo$gene[ which( cellType_markers_Embryo$cluster == Embryo_cluster)][1:100]
    
    return (length(intersect(EB_markers, Embryo_markers)) / length (unique(EB_markers, Embryo_markers)))

  })
}) -> marker_overlap_mat

colnames(marker_overlap_mat) <- unique(cellType_markers_EB$cluster)
rownames(marker_overlap_mat) <- unique(cellType_markers_Embryo$cluster)

marker_overlap_df <- reshape2::melt( marker_overlap_mat )

sample_order <- c("RGP", "IP", "aNSC", "astro", "OBNB", "oligo", "CR", "iN")

marker_overlap_df$Var1 <- factor( marker_overlap_df$Var1, levels =sample_order )
marker_overlap_df$Var2 <- factor( marker_overlap_df$Var2, levels =sample_order )

colnames( marker_overlap_df ) <- c("Embryo", "EB", "frac" )

ggplot( marker_overlap_df, aes( x=Embryo, y=EB, fill=frac)) + 
  geom_tile() + scale_fill_gradient( low = "white", high = "black") + theme_classic()
ggsave( filename = paste(base, "plots/Fig_1/F1h_marker_overlap.pdf", sep=""))

# heatmap with overlap genes
sapply( sample_order, function(cellType){
  intersect( cellType_markers_EB$gene[ which( cellType_markers_EB$cluster == cellType ) ][1:20],
             cellType_markers_Embryo$gene[ which( cellType_markers_Embryo$cluster == cellType)][1:20] )
}) -> common_markers

common_markers <- unlist(common_markers) [ which( unlist(common_markers) %in% names(which(table(unlist(common_markers)) == 1)) ) ]
seurat_sketch_EB$cellType <- factor( seurat_sketch_EB$cellType, levels = sample_order )
seurat_sketch_EB <- ScaleData( seurat_sketch_EB, features = common_markers )
seurat_sketch_embryo <- ScaleData( seurat_sketch_embryo, features = common_markers )
seurat_sketch_embryo$cellType <- factor( seurat_sketch_embryo$cellType, levels = sample_order )

organoid_merker_heamap <- DoHeatmap( object = seurat_sketch_EB, features = common_markers, group.by = "cellType") + ggtitle("organoid")
embryo_merker_heamap <- DoHeatmap( object = seurat_sketch_embryo, features = common_markers, group.by = "cellType") + ggtitle("embryo")
comb_heatmap <- organoid_merker_heamap + embryo_merker_heamap
ggsave(  filename = paste(base, "plots/Sup_Fig_7/S7b_marker_overlap_heatmap.pdf", sep=""), width = 14, height = 8 )

# cell type marker genes analysis
cluster2cellType <- as.character( c(17,10,1,0,7,20,27,15,11,
                                    14,19,18,23,29,8,
                                    21,9,13,2,3,
                                    12,30,22,24,28,25,16,
                                    4, 5, 31, 6, 32, 26) )

names( cluster2cellType ) <- c( rep("iN", 9), rep("IP", 6), rep("RGP", 5), rep("aNSC", 7), 
                                "OBNB", "oligo", "oligo", "astro", "astro","CR" )

# collect overlapping marker genes for reporting
all_markers <- data.frame()
for ( cell_type in unique( cellType_markers_EB$cluster ) ) {
  EB_markers <- cellType_markers_EB$gene[ which( cellType_markers_EB$cluster == cell_type )][1:100]
  Embryo_markers <- cellType_markers_Embryo$gene[ which( cellType_markers_Embryo$cluster == cell_type)][1:100]
  
  gene_list <- intersect( EB_markers, Embryo_markers)
  
  all_markers <- rbind(all_markers, data.frame( cell_type = cell_type, genes = gene_list))
}

# prepare a data frame of seurat cluster to cell type association for reporting
cluster_to_cellType <- as.data.frame( cluster2cellType )
cluster_to_cellType$cellType <- names( cluster2cellType )
colnames( cluster_to_cellType ) <- c("seurat_cluster", "cell type")

wb <- createWorkbook()

  sheetName <- "cluster to cellType"
  
  addWorksheet(wb, sheetName)
  
  writeData( wb, sheetName, cluster_to_cellType )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol(cluster_to_cellType) )
  setColWidths( wb, sheetName, cols = 1:ncol(cluster_to_cellType), widths="auto" )

  sheetName <- "marker genes"
  
  addWorksheet(wb, sheetName)
  
  writeData( wb, sheetName, all_markers )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol(all_markers) )
  setColWidths( wb, sheetName, cols = 1:ncol(all_markers), widths="auto" )

saveWorkbook(wb, file = paste(base, "/plots/Sup_Fig_7/cellType_markers.xlsx", sep=""), overwrite = T) 

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] openxlsx_4.2.5.2   GeneOverlap_1.38.0 simspec_0.0.0.9000 ggplot2_3.5.0      uwot_0.1.16        Matrix_1.6-5       Seurat_5.0.1      
# [8] SeuratObject_5.0.1 sp_2.1-3          
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     rstudioapi_0.16.0      jsonlite_1.8.8         magrittr_2.0.3         spatstat.utils_3.1-0  
# [6] farver_2.1.1           ragg_1.3.0             vctrs_0.6.5            ROCR_1.0-11            spatstat.explore_3.2-6
# [11] htmltools_0.5.7        sctransform_0.4.1      parallelly_1.37.0      KernSmooth_2.23-22     htmlwidgets_1.6.4     
# [16] ica_1.0-3              plyr_1.8.9             plotly_4.10.4          zoo_1.8-12             igraph_2.1.4          
# [21] mime_0.12              lifecycle_1.0.4        pkgconfig_2.0.3        R6_2.5.1               fastmap_1.1.1         
# [26] fitdistrplus_1.1-11    future_1.33.1          shiny_1.8.0            digest_0.6.34          colorspace_2.1-0      
# [31] patchwork_1.2.0        tensor_1.5             RSpectra_0.16-1        irlba_2.3.5.1          textshaping_0.3.7     
# [36] labeling_0.4.3         progressr_0.14.0       fansi_1.0.6            spatstat.sparse_3.0-3  httr_1.4.7            
# [41] polyclip_1.10-6        abind_1.4-5            compiler_4.3.2         withr_3.0.0            fastDummies_1.7.3     
# [46] gplots_3.1.3.1         MASS_7.3-60            gtools_3.9.5           caTools_1.18.2         tools_4.3.2           
# [51] lmtest_0.9-40          zip_2.3.1              httpuv_1.6.14          future.apply_1.11.1    goftest_1.2-3         
# [56] glue_1.7.0             nlme_3.1-163           promises_1.2.1         grid_4.3.2             Rtsne_0.17            
# [61] cluster_2.1.4          reshape2_1.4.4         generics_0.1.3         gtable_0.3.4           spatstat.data_3.0-4   
# [66] tidyr_1.3.1            data.table_1.15.0      utf8_1.2.4             spatstat.geom_3.2-8    RcppAnnoy_0.0.22      
# [71] ggrepel_0.9.5          RANN_2.6.1             pillar_1.9.0           stringr_1.5.1          limma_3.58.1          
# [76] spam_2.10-0            RcppHNSW_0.6.0         later_1.3.2            splines_4.3.2          dplyr_1.1.4           
# [81] lattice_0.21-9         survival_3.5-7         deldir_2.0-2           tidyselect_1.2.1       miniUI_0.1.1.1        
# [86] pbapply_1.7-2          gridExtra_2.3          scattermore_1.2        statmod_1.5.0          matrixStats_1.2.0     
# [91] stringi_1.8.3          lazyeval_0.2.2         codetools_0.2-19       tibble_3.2.1           cli_3.6.2             
# [96] xtable_1.8-4           reticulate_1.44.1      systemfonts_1.0.6      munsell_0.5.0          Rcpp_1.0.12           
# [101] globals_0.16.2         spatstat.random_3.2-2  png_0.1-8              parallel_4.3.2         ellipsis_0.3.2        
# [106] presto_1.0.0           dotCall64_1.1-1        bitops_1.0-7           listenv_0.9.1          viridisLite_0.4.2     
# [111] scales_1.3.0           ggridges_0.5.6         leiden_0.4.3.1         packrat_0.9.2          purrr_1.0.2           
# [116] rlang_1.1.3            cowplot_1.1.3  
