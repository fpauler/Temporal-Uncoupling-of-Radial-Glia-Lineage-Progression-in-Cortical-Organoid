# =====================================================
# step 02
# calculate significance of GFP expression differences
# =====================================================

library (Seurat)
library (dplyr)
library (ggplot2)
library (pheatmap)

set.seed(2401)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

seurat.obj.all <- readRDS( file =  paste( nobackup_base,  "/RDS_files/seurat.embryo.RDS", sep="" ) )
Idents(seurat.obj.all) <- "orig.ident"

for ( ID in unique(seurat.obj.all@meta.data$orig.ident) ) {
  
  message ( ID )
  
  seurat.obj <- subset( seurat.obj.all, idents = ID )
  
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)

  seurat.obj <- ScaleData(seurat.obj)
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  elbowPlot <- ElbowPlot( object = seurat.obj, ndims = 30)
  
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:25)
  seurat.obj <- FindClusters(seurat.obj, resolution = 1)
  
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:25)
  
  umap_plot <- DimPlot( seurat.obj, reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend()
  feature_plot <- FeaturePlot( seurat.obj, features = c("GFP", "Gad1", "Rbfox1", "Aldh1l1", "Hes5", "Mki67"), order = T, min.cutoff = "q25")
  
  ggsave( plot = umap_plot, file = paste ( base, "plots/QC/", ID, ".GFP_map.umap.pdf", sep="" ) )
  ggsave( plot = feature_plot, file = paste ( base, "plots/QC/", ID, ".GFP_map.feature_plot.pdf", sep="" ), width = 10, height = 10 )
  
  cluster_vec <- sort( unique(seurat.obj@meta.data$seurat_clusters) ) 
  
  all_comb <- combn( as.character(cluster_vec), 2)
  
  all_p_df <- data.frame()
  for (idx in 1:ncol(all_comb)) {
    tmp <- FindMarkers( object = seurat.obj, ident.1 = all_comb[1,idx], ident.2 = all_comb[2,idx], features = "GFP", test.use = "bimod" )
    if (nrow(tmp) > 0) {
      if ( abs(tmp$avg_log2FC) > 1) {
        pval  <- tmp$p_val_adj
      } else {
        pval  <- 1
      }
    } else {
      pval  <- 1
    }
    all_p_df <- rbind( all_p_df,data.frame( x=all_comb[1,idx], y=all_comb[2,idx], pval = pval ) )
  }
  
  sig_mat <- sapply( cluster_vec, function (ident1){
     sapply( cluster_vec, function (ident2){
       tmp <- all_p_df[ which(all_p_df$x == ident1 & all_p_df$y == ident2), ]
       if ( nrow(tmp) == 0 ) {
         tmp <- all_p_df[ which(all_p_df$x == ident2 & all_p_df$y == ident1), ]  
       }
       
       if ( nrow(tmp) == 0 ) {
         return(1)  
       } else {
         return(tmp$pval)
       }
       
     })
   })
   
  # check reproducibility
  # sig_mat_old <- readRDS ( file = paste( nobackup_base, "RDS_files/", ID, ".GFP.pvalue.mar.RDS", sep = "") )
  # identical(sig_mat, sig_mat_old)
  
   saveRDS (object = sig_mat, file = paste( nobackup_base, "RDS_files/", ID, ".GFP.pvalue.mar.RDS", sep = "") )

   sig_mat <- log10(sig_mat) *-1
   rownames(sig_mat) <- as.character( 0 : ( ncol(sig_mat) - 1 ) )
   colnames(sig_mat) <- as.character( 0 : ( ncol(sig_mat) - 1 ) )
   
   # exchange infinite value with some reasonable value to plot on the heatmap
   sig_mat[ which(is.infinite(sig_mat)) ] <- max (sig_mat[ which(!is.infinite(sig_mat)) ])
   pheatmap( sig_mat, scale="none", cluster_rows = T, cluster_cols = T, 
             filename = paste( base, "plots/QC/", ID,".eGFP.pvalue.pdf", sep="" ), main = paste( ID, "p-value plot" ))
   
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
#   [1] pheatmap_1.0.12    ggplot2_3.5.0      dplyr_1.1.4        Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-3          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-2           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.3            magrittr_2.0.3        
# [6] RcppAnnoy_0.0.22       spatstat.geom_3.2-8    matrixStats_1.2.0      ggridges_0.5.6         compiler_4.3.2        
# [11] png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3       
# [16] fastmap_1.1.1          ellipsis_0.3.2         labeling_0.4.3         utf8_1.2.4             promises_1.2.1        
# [21] purrr_1.0.2            jsonlite_1.8.8         goftest_1.2-3          later_1.3.2            spatstat.utils_3.1-0  
# [26] irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.4          R6_2.5.1               ica_1.0-3             
# [31] stringi_1.8.3          RColorBrewer_1.1-3     spatstat.data_3.0-4    reticulate_1.35.0      parallelly_1.37.0     
# [36] lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.12            tensor_1.5             future.apply_1.11.1   
# [41] zoo_1.8-12             sctransform_0.4.1      httpuv_1.6.14          Matrix_1.6-5           splines_4.3.2         
# [46] igraph_2.1.4           tidyselect_1.2.1       rstudioapi_0.16.0      abind_1.4-5            spatstat.random_3.2-2 
# [51] codetools_0.2-19       miniUI_0.1.1.1         spatstat.explore_3.2-6 listenv_0.9.1          lattice_0.21-9        
# [56] tibble_3.2.1           plyr_1.8.9             withr_3.0.0            shiny_1.8.0            ROCR_1.0-11           
# [61] Rtsne_0.17             future_1.33.1          fastDummies_1.7.3      survival_3.5-7         polyclip_1.10-6       
# [66] fitdistrplus_1.1-11    pillar_1.9.0           packrat_0.9.2          KernSmooth_2.23-22     plotly_4.10.4         
# [71] generics_0.1.3         RcppHNSW_0.6.0         munsell_0.5.0          scales_1.3.0           globals_0.16.2        
# [76] xtable_1.8-4           glue_1.7.0             lazyeval_0.2.2         tools_4.3.2            data.table_1.15.0     
# [81] RSpectra_0.16-1        RANN_2.6.1             leiden_0.4.3.1         dotCall64_1.1-1        cowplot_1.1.3         
# [86] grid_4.3.2             tidyr_1.3.1            colorspace_2.1-0       nlme_3.1-163           patchwork_1.2.0       
# [91] cli_3.6.2              spatstat.sparse_3.0-3  spam_2.10-0            fansi_1.0.6            viridisLite_0.4.2     
# [96] uwot_0.1.16            gtable_0.3.4           digest_0.6.34          progressr_0.14.0       ggrepel_0.9.5         
# [101] farver_2.1.1           htmlwidgets_1.6.4      htmltools_0.5.7        lifecycle_1.0.4        httr_1.4.7            
# [106] mime_0.12              MASS_7.3-60
