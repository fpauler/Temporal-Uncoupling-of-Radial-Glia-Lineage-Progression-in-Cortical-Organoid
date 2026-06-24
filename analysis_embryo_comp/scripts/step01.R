# ===================================================
# step 01
# initial QC and filtering - embryo
# ===================================================

library (Seurat)
library (ggplot2)
library (openxlsx)

set.seed(2401)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

#read the sample list
sample_list <- read.xlsx( paste(base, "other_data/sample_list.xlsx", sep=""), sheet = 2)
sample_list <- sample_list[ 17:20 , ]

for ( ID in unique(sample_list$ID) ) {
  
  tmp <- sample_list[ which(sample_list$ID == ID), ]
  poolID <- unique( tmp$poolID )
  
  tmp_df <- read.csv( file = paste(base, "other_data/4plex_", poolID, "_Multiplex_config.csv", sep="") )
  
  conv_vec <- tmp_df[ c(11:14), 3]
  names(conv_vec) <- tmp_df[ c(11:14), 1] 
  
  seurat_obj_list <- list()
  
  path <- paste( nobackup_base, "h5_files/", sep="" )
  
  # expression data for each set of replicates is read into a list for later merging
  # that way I create Seurat objects with layers for integration
  for ( sampleID in as.character(tmp$sampleID) ) {
    
    message( ID, " ", sampleID)
    
    h5_file_path <- paste( nobackup_base, "h5_files/", sampleID, "_sample_filtered_feature_bc_matrix.h5", sep="" )
    message( h5_file_path )  
    
    counts <- Read10X_h5( filename = h5_file_path )
    seurat_obj_list[[sampleID]] <- CreateSeuratObject(counts = counts, project = conv_vec[ sampleID ], min.cells = 3, min.features = 200)
    
  }
  
  message("merging ages")
  
  seurat.obj <- merge( x = seurat_obj_list[[ 1 ]], y = c(seurat_obj_list[[ 2 ]], seurat_obj_list[[ 3 ]], seurat_obj_list[[ 4 ]] ) )
  
  message("preparing QC plots")
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
  
  # Visualize QC metrics as a violin plot
  vln_plot <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave( plot=vln_plot, file=paste (base, "plots/QC/", ID, ".Vln.raw.pdf", sep=""), width = 10, height = 7)
  
  message("filtering")
  seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
  vln_plot <- VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  ggsave( plot=vln_plot, file=paste (base, "plots/QC/", ID, ".Vln.filtered.pdf", sep=""), width = 10, height = 7)
  
  plot1 <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  QC_plot <- plot1 + plot2
  
  ggsave( plot=QC_plot, file=paste (base, "plots/QC/", ID, ".QC.pdf", sep=""), width = 10, height = 7)
  
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  
  all.genes <- rownames(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj, features = all.genes)
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  elbowPlot <- ElbowPlot( object = seurat.obj, ndims = 30)
  
  ggsave( plot = elbowPlot, file = paste (base, "plots/QC/", ID, ".elbow.pdf", sep=""), width = 7, height = 7)
  
  saveRDS( seurat.obj, paste( nobackup_base,  "/RDS_files/seurat.", ID, ".RDS", sep="") )
  
 }

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
#   [1] openxlsx_4.2.5.2   ggplot2_3.5.0      Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-3          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-2           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.3            magrittr_2.0.3        
# [6] RcppAnnoy_0.0.22       spatstat.geom_3.2-8    matrixStats_1.2.0      ggridges_0.5.6         compiler_4.3.2        
# [11] png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         hdf5r_1.3.10           stringr_1.5.1         
# [16] pkgconfig_2.0.3        fastmap_1.1.1          ellipsis_0.3.2         labeling_0.4.3         utf8_1.2.4            
# [21] promises_1.2.1         ggbeeswarm_0.7.2       bit_4.0.5              purrr_1.0.2            jsonlite_1.8.8        
# [26] goftest_1.2-3          later_1.3.2            spatstat.utils_3.1-0   irlba_2.3.5.1          parallel_4.3.2        
# [31] cluster_2.1.4          R6_2.5.1               ica_1.0-3              stringi_1.8.3          RColorBrewer_1.1-3    
# [36] spatstat.data_3.0-4    reticulate_1.35.0      parallelly_1.37.0      lmtest_0.9-40          scattermore_1.2       
# [41] Rcpp_1.0.12            tensor_1.5             future.apply_1.11.1    zoo_1.8-12             sctransform_0.4.1     
# [46] httpuv_1.6.14          Matrix_1.6-5           splines_4.3.2          igraph_2.1.4           tidyselect_1.2.1      
# [51] rstudioapi_0.16.0      abind_1.4-5            spatstat.random_3.2-2  codetools_0.2-19       miniUI_0.1.1.1        
# [56] spatstat.explore_3.2-6 listenv_0.9.1          lattice_0.21-9         tibble_3.2.1           plyr_1.8.9            
# [61] withr_3.0.0            shiny_1.8.0            ROCR_1.0-11            ggrastr_1.0.2          Rtsne_0.17            
# [66] future_1.33.1          fastDummies_1.7.3      survival_3.5-7         polyclip_1.10-6        zip_2.3.1             
# [71] fitdistrplus_1.1-11    pillar_1.9.0           packrat_0.9.2          KernSmooth_2.23-22     plotly_4.10.4         
# [76] generics_0.1.3         RcppHNSW_0.6.0         munsell_0.5.0          scales_1.3.0           globals_0.16.2        
# [81] xtable_1.8-4           glue_1.7.0             lazyeval_0.2.2         tools_4.3.2            data.table_1.15.0     
# [86] RSpectra_0.16-1        RANN_2.6.1             leiden_0.4.3.1         dotCall64_1.1-1        cowplot_1.1.3         
# [91] grid_4.3.2             tidyr_1.3.1            colorspace_2.1-0       nlme_3.1-163           patchwork_1.2.0       
# [96] beeswarm_0.4.0         vipor_0.4.7            cli_3.6.2              spatstat.sparse_3.0-3  spam_2.10-0           
# [101] fansi_1.0.6            viridisLite_0.4.2      dplyr_1.1.4            uwot_0.1.16            gtable_0.3.4          
# [106] digest_0.6.34          progressr_0.14.0       ggrepel_0.9.5          farver_2.1.1           htmlwidgets_1.6.4     
# [111] htmltools_0.5.7        lifecycle_1.0.4        httr_1.4.7             mime_0.12              bit64_4.0.5           
# [116] MASS_7.3-60   
