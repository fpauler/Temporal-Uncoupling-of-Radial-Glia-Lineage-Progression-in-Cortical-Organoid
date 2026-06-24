# ===================================
# step 03
# assign GFP tag for later filtering
# ===================================

library (Seurat)
library (ggplot2)
library (SeuratObject)
library (org.Mm.eg.db)
library (clusterProfiler)
library (simspec)
library (patchwork)

set.seed(2401)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

clusters2remove <- list()

# at E10 the GFP positive cells are the minority
# p-value heatmap reveals GFP positive clusters
clusters2remove[["E10"]] <- setdiff( c(0:21), c("10", "1", "6", "16", "11", "7", "9") )

# at E13, E16, P0 GFP positive cells are the majority
# p-value heatmap reveals GFP negative cells - also remove small cell populations
clusters2remove[["E13"]] <- c("13", "14", "15", "17", "19", "20", "21", "22", "23")
clusters2remove[["E16"]] <- c("10", "12", "13", "14", "16", "17", "18", "19", "20")
clusters2remove[["P0"]] <- c("12", "13", "14", "15", "16", "17", "19", "20", "21")
  
seurat.obj.all <- readRDS( file =  paste( nobackup_base,  "/RDS_files/seurat.embryo.RDS", sep="" ) )
Idents(seurat.obj.all) <- "orig.ident"
  
seurat_obj_list <- list()
  
for ( ID in c("E10", "E13", "E16", "P0") ) {
    
    message ( ID )
    
    seurat.obj <- subset( seurat.obj.all, idents = ID )
    
    seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
    
    seurat.obj <- ScaleData(seurat.obj)
    seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
    elbowPlot <- ElbowPlot( object = seurat.obj, ndims = 30)
    
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:25)
    seurat.obj <- FindClusters(seurat.obj, resolution = 1)
    
    small_clusters <- names ( which( table( seurat.obj@meta.data$seurat_clusters) < 100 ) )
    
    # add information about GFP expression for supplemental figure
    seurat.obj@meta.data$GFP_cluster <- "YES"
    seurat.obj@meta.data$GFP_cluster[ which( as.character(seurat.obj@meta.data$seurat_clusters) %in% c( clusters2remove[[ ID ]], small_clusters) ) ] <- "NO"
    
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:25)
    
    plot1 <- DimPlot( seurat.obj, reduction = "umap", group.by = "GFP_cluster", label = F)
    plot2 <- FeaturePlot( seurat.obj, features = c("GFP"), order = T, min.cutoff = "q25") + NoLegend()
    comb_plot <- plot1+plot2
    ggsave( plot = comb_plot, file = paste ( base, "plots/QC/", ID, ".GFP.annotation.pdf", sep="" ), width = 20, height = 10 )
    
    Idents(seurat.obj) <- "GFP_cluster"
    #extract GFP positive cells
    seurat_obj_list[[ID]] <- subset( seurat.obj, idents = "YES" )
  
    umap_plot <- DimPlot( seurat_obj_list[[ID]], reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend()
    ggsave( plot = umap_plot, file = paste ( base, "plots/QC/", ID, ".GFP.extracted.pdf", sep="" ), width = 10, height = 10 )
  
  }
  
saveRDS( object = seurat_obj_list, file = paste( nobackup_base, "RDS_files/embryo_GFP_perAge.RDS", sep = "") )

# $E10
# An object of class Seurat 
# 15629 features across 3658 samples within 1 assay 
# Active assay: RNA (15629 features, 2000 variable features)
# 4 layers present: counts.E10, data.E10, scale.data, data
# 2 dimensional reductions calculated: pca, umap
# 
# $E13
# An object of class Seurat 
# 15151 features across 9456 samples within 1 assay 
# Active assay: RNA (15151 features, 1916 variable features)
# 4 layers present: counts.E13, data.E13, scale.data, data
# 2 dimensional reductions calculated: pca, umap
# 
# $E16
# An object of class Seurat 
# 14700 features across 3637 samples within 1 assay 
# Active assay: RNA (14700 features, 1866 variable features)
# 4 layers present: counts.E16, data.E16, scale.data, data
# 2 dimensional reductions calculated: pca, umap
# 
# $P0
# An object of class Seurat 
# 15623 features across 9009 samples within 1 assay 
# Active assay: RNA (15623 features, 1933 variable features)
# 4 layers present: counts.P0, data.P0, scale.data, data
# 2 dimensional reductions calculated: pca, umap

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
#   [1] patchwork_1.2.0        simspec_0.0.0.9000     clusterProfiler_4.10.1 org.Mm.eg.db_3.18.0    AnnotationDbi_1.64.1  
# [6] IRanges_2.36.0         S4Vectors_0.40.2       Biobase_2.62.0         BiocGenerics_0.48.1    ggplot2_3.5.0         
# [11] Seurat_5.0.1           SeuratObject_5.0.1     sp_2.1-3              
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.3.2           later_1.3.2             ggplotify_0.1.2         bitops_1.0-7           
# [6] tibble_3.2.1            polyclip_1.10-6         fastDummies_1.7.3       lifecycle_1.0.4         globals_0.16.2         
# [11] lattice_0.21-9          MASS_7.3-60             magrittr_2.0.3          plotly_4.10.4           httpuv_1.6.14          
# [16] sctransform_0.4.1       spam_2.10-0             spatstat.sparse_3.0-3   reticulate_1.35.0       cowplot_1.1.3          
# [21] pbapply_1.7-2           DBI_1.2.2               RColorBrewer_1.1-3      abind_1.4-5             zlibbioc_1.48.2        
# [26] Rtsne_0.17              purrr_1.0.2             ggraph_2.2.1            RCurl_1.98-1.14         yulab.utils_0.1.4      
# [31] tweenr_2.0.3            GenomeInfoDbData_1.2.11 enrichplot_1.22.0       ggrepel_0.9.5           irlba_2.3.5.1          
# [36] listenv_0.9.1           spatstat.utils_3.1-0    tidytree_0.4.6          goftest_1.2-3           RSpectra_0.16-1        
# [41] spatstat.random_3.2-2   fitdistrplus_1.1-11     parallelly_1.37.0       leiden_0.4.3.1          codetools_0.2-19       
# [46] DOSE_3.28.2             ggforce_0.4.2           tidyselect_1.2.1        aplot_0.2.2             farver_2.1.1           
# [51] viridis_0.6.5           matrixStats_1.2.0       spatstat.explore_3.2-6  jsonlite_1.8.8          ellipsis_0.3.2         
# [56] tidygraph_1.3.1         progressr_0.14.0        ggridges_0.5.6          survival_3.5-7          systemfonts_1.0.6      
# [61] tools_4.3.2             ragg_1.3.0              treeio_1.26.0           ica_1.0-3               Rcpp_1.0.12            
# [66] glue_1.7.0              gridExtra_2.3           qvalue_2.34.0           GenomeInfoDb_1.38.8     dplyr_1.1.4            
# [71] withr_3.0.0             fastmap_1.1.1           fansi_1.0.6             digest_0.6.34           gridGraphics_0.5-1     
# [76] R6_2.5.1                mime_0.12               textshaping_0.3.7       colorspace_2.1-0        scattermore_1.2        
# [81] GO.db_3.18.0            tensor_1.5              spatstat.data_3.0-4     RSQLite_2.3.6           utf8_1.2.4             
# [86] tidyr_1.3.1             generics_0.1.3          data.table_1.15.0       graphlayouts_1.1.1      httr_1.4.7             
# [91] htmlwidgets_1.6.4       scatterpie_0.2.2        uwot_0.1.16             pkgconfig_2.0.3         gtable_0.3.4           
# [96] blob_1.2.4              lmtest_0.9-40           XVector_0.42.0          shadowtext_0.1.3        htmltools_0.5.7        
# [101] dotCall64_1.1-1         fgsea_1.28.0            scales_1.3.0            png_0.1-8               ggfun_0.1.4            
# [106] rstudioapi_0.16.0       reshape2_1.4.4          nlme_3.1-163            cachem_1.0.8            zoo_1.8-12             
# [111] stringr_1.5.1           KernSmooth_2.23-22      parallel_4.3.2          miniUI_0.1.1.1          HDO.db_0.99.1          
# [116] pillar_1.9.0            grid_4.3.2              vctrs_0.6.5             RANN_2.6.1              promises_1.2.1         
# [121] xtable_1.8-4            cluster_2.1.4           packrat_0.9.2           cli_3.6.2               compiler_4.3.2         
# [126] rlang_1.1.3             crayon_1.5.2            future.apply_1.11.1     labeling_0.4.3          plyr_1.8.9             
# [131] fs_1.6.3                stringi_1.8.3           viridisLite_0.4.2       deldir_2.0-2            BiocParallel_1.36.0    
# [136] munsell_0.5.0           Biostrings_2.70.3       lazyeval_0.2.2          spatstat.geom_3.2-8     GOSemSim_2.28.1        
# [141] Matrix_1.6-5            RcppHNSW_0.6.0          bit64_4.0.5             future_1.33.1           KEGGREST_1.42.0        
# [146] shiny_1.8.0             ROCR_1.0-11             igraph_2.1.4            memoise_2.0.1           ggtree_3.10.1          
# [151] fastmatch_1.1-4         bit_4.0.5               gson_0.1.0              ape_5.8   
