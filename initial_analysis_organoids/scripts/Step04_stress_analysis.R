# ========================================
# Step 4 
# prepare data for further analysis
# ========================================

library (Seurat)
library (ggplot2)
library (clusterProfiler)
library (org.Mm.eg.db)
library (SeuratWrappers)
library (monocle3)
library (future)
library (dplyr)
library (viridis)

set.seed(2401)

plan("multicore", workers = 6)
options(future.globals.maxSize= 4000000000)

base <- "~/initial_analysis_organoids/"
nobackup_base <- "~/noSave/"

# ggplot color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

if ( !file.exists(paste( nobackup_base, "RDS_files/seurat_obj.stress.CC_correct.RDS", sep = "")) ) {
  
  # read the meta data with the cell cycle scores and add to the object
  seurat.obj <- readRDS( paste( nobackup_base, "RDS_files/seurat_obj.stress.removed.RDS", sep = "") )
  CC_meta_data <- readRDS ( file = paste( nobackup_base, "RDS_files/seurat_obj.stress.CC_added.metaData.RDS", sep = "") )
  CC_meta_data$CC.Difference <- CC_meta_data$S.Score - CC_meta_data$G2M.Score
  
  seurat.obj <- AddMetaData(object = seurat.obj, metadata = CC_meta_data[,c("S.Score", "G2M.Score", "CC.Difference")])
  
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 4000)
  
  seurat.obj <- ScaleData( seurat.obj, vars.to.regress = "CC.Difference" )
  
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  
  seurat.obj <- FindNeighbors(seurat.obj, reduction = "pca", dims = 1:20)
  
  # Note: this clustering was used for similarity spectrum integration
  seurat.obj <- FindClusters(seurat.obj, resolution = 0.5, cluster.name = "clusters_unintegrated")
  
  seurat.obj <- RunUMAP(seurat.obj, reduction = "pca", dims = 1:20, reduction.name = "umap")
  
  # removing scaled data reduces object size
  seurat.obj@assays$RNA$scale.data <- NULL
  
  saveRDS( object = seurat.obj, file = paste( nobackup_base, "RDS_files/seurat_obj.stress.CC_correct.RDS", sep = "") )
  
}

# analysis continues in analysis_embryo_comp folder

sessionInfo()

# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /opt/R/4.3.2/lib/R/lib/libRblas.so 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] dplyr_1.1.4                 future_1.33.1               monocle3_1.3.7              SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0
# [6] GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         MatrixGenerics_1.14.0       matrixStats_1.2.0           SeuratWrappers_0.3.5       
# [11] org.Mm.eg.db_3.18.0         AnnotationDbi_1.64.1        IRanges_2.36.0              S4Vectors_0.40.2            Biobase_2.62.0             
# [16] BiocGenerics_0.48.1         clusterProfiler_4.10.1      ggplot2_3.5.0               Seurat_5.0.1                SeuratObject_5.0.1         
# [21] sp_2.1-3                   
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.3                spatstat.sparse_3.0-3   bitops_1.0-7            enrichplot_1.22.0       HDO.db_0.99.1           httr_1.4.7             
# [7] RColorBrewer_1.1-3      tools_4.3.2             sctransform_0.4.1       utf8_1.2.4              R6_2.5.1                lazyeval_0.2.2         
# [13] uwot_0.1.16             withr_3.0.0             gridExtra_2.3           progressr_0.14.0        cli_3.6.2               textshaping_0.3.7      
# [19] spatstat.explore_3.2-6  fastDummies_1.7.3       scatterpie_0.2.2        labeling_0.4.3          spatstat.data_3.0-4     ggridges_0.5.6         
# [25] pbapply_1.7-2           systemfonts_1.0.6       yulab.utils_0.1.4       gson_0.1.0              DOSE_3.28.2             R.utils_2.12.3         
# [31] harmony_1.2.0           parallelly_1.37.0       limma_3.58.1            rstudioapi_0.16.0       RSQLite_2.3.6           generics_0.1.3         
# [37] gridGraphics_0.5-1      ica_1.0-3               spatstat.random_3.2-2   GO.db_3.18.0            Matrix_1.6-5            ggbeeswarm_0.7.2       
# [43] fansi_1.0.6             abind_1.4-5             R.methodsS3_1.8.2       lifecycle_1.0.4         qvalue_2.34.0           SparseArray_1.2.4      
# [49] Rtsne_0.17              grid_4.3.2              blob_1.2.4              promises_1.2.1          crayon_1.5.2            miniUI_0.1.1.1         
# [55] lattice_0.21-9          cowplot_1.1.3           KEGGREST_1.42.0         pillar_1.9.0            fgsea_1.28.0            boot_1.3-28.1          
# [61] future.apply_1.11.1     codetools_0.2-19        fastmatch_1.1-4         leiden_0.4.3.1          glue_1.7.0              packrat_0.9.2          
# [67] ggfun_0.1.4             data.table_1.15.0       remotes_2.5.0           vctrs_0.6.5             png_0.1-8               treeio_1.26.0          
# [73] spam_2.10-0             gtable_0.3.4            cachem_1.0.8            S4Arrays_1.2.1          mime_0.12               tidygraph_1.3.1        
# [79] survival_3.5-7          statmod_1.5.0           ellipsis_0.3.2          fitdistrplus_1.1-11     ROCR_1.0-11             nlme_3.1-163           
# [85] ggtree_3.10.1           bit64_4.0.5             RcppAnnoy_0.0.22        irlba_2.3.5.1           vipor_0.4.7             KernSmooth_2.23-22     
# [91] colorspace_2.1-0        DBI_1.2.2               ggrastr_1.0.2           tidyselect_1.2.1        bit_4.0.5               compiler_4.3.2         
# [97] DelayedArray_0.28.0     plotly_4.10.4           shadowtext_0.1.3        scales_1.3.0            lmtest_0.9-40           stringr_1.5.1          
# [103] digest_0.6.34           goftest_1.2-3           presto_1.0.0            spatstat.utils_3.0-4    minqa_1.2.6             XVector_0.42.0         
# [109] RhpcBLASctl_0.23-42     htmltools_0.5.7         pkgconfig_2.0.3         lme4_1.1-35.3           fastmap_1.1.1           rlang_1.1.3            
# [115] htmlwidgets_1.6.4       shiny_1.8.0             farver_2.1.1            zoo_1.8-12              jsonlite_1.8.8          BiocParallel_1.36.0    
# [121] GOSemSim_2.28.1         R.oo_1.26.0             RCurl_1.98-1.14         magrittr_2.0.3          GenomeInfoDbData_1.2.11 ggplotify_0.1.2        
# [127] dotCall64_1.1-1         patchwork_1.2.0         munsell_0.5.0           Rcpp_1.0.12             ape_5.8                 viridis_0.6.5          
# [133] reticulate_1.35.0       stringi_1.8.3           ggraph_2.2.1            zlibbioc_1.48.2         MASS_7.3-60             plyr_1.8.9             
# [139] parallel_4.3.2          listenv_0.9.1           ggrepel_0.9.5           deldir_2.0-2            Biostrings_2.70.3       graphlayouts_1.1.1     
# [145] splines_4.3.2           tensor_1.5              igraph_2.0.2            spatstat.geom_3.2-8     RcppHNSW_0.6.0          reshape2_1.4.4         
# [151] BiocManager_1.30.22     nloptr_2.0.3            tweenr_2.0.3            httpuv_1.6.14           RANN_2.6.1              tidyr_1.3.1            
# [157] purrr_1.0.2             polyclip_1.10-6         scattermore_1.2         ggforce_0.4.2           rsvd_1.0.5              xtable_1.8-4           
# [163] RSpectra_0.16-1         tidytree_0.4.6          later_1.3.2             viridisLite_0.4.2       ragg_1.3.0              tibble_3.2.1           
# [169] aplot_0.2.2             memoise_2.0.1           beeswarm_0.4.0          cluster_2.1.4           globals_0.16.2 
