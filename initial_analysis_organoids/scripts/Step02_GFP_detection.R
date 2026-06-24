# ===============================================================
# Step 2
# detection of GFP containing clusters via pairwise comparisons
# ===============================================================

library (Seurat)
library (dplyr)
library (ggplot2)
library (pheatmap)

set.seed(2401)

# define working folder
base <- "~/initial_analysis_organoids/"
nobackup_base <- "~/noSave/"

for ( ID in c("D8", "D13", "D20", "D25" ) ) {
  
  message ( ID )
  
  seurat.obj <- readRDS( file =  paste( nobackup_base,  "/RDS_files/seurat.", ID, ".RDS", sep="" ) )
  
  seurat.obj <- JoinLayers( seurat.obj )
  
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:25)
  seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)
  
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:25)
  
  umap_plot <- DimPlot( seurat.obj, reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend()
  feature_plot <- FeaturePlot( seurat.obj, features = c("GFP", "Gad1", "Rbfox1", "Aldh1l1", "Hes5", "Mki67"), order = T, min.cutoff = "q25")
  
  ggsave( plot = umap_plot, file = paste ( base, "plots/QC/", ID, ".GFP_map.umap.pdf", sep="" ) )
  ggsave( plot = feature_plot, file = paste ( base, "plots/QC/", ID, ".GFP_map.feature_plot.pdf", sep="" ), width = 10, height = 10 )
  
  cluster_vec <- sort( unique(seurat.obj@meta.data$seurat_clusters) ) 
  
  # determine all pairwise comparisons to do
  all_comb <- combn( as.character(cluster_vec), 2)
  
  # prepare pairwise GFP expression comparison
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
  
  # transform the dataframe from above into a matrix
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
  
   saveRDS (object = sig_mat, file = paste( nobackup_base, "RDS_files/", ID, ".GFP.pvalue.mar.RDS", sep = "") )

}

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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
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
#   [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8          magrittr_2.0.3          spatstat.utils_3.0-4    farver_2.1.1            ragg_1.3.0             
# [8] zlibbioc_1.48.2         vctrs_0.6.5             ROCR_1.0-11             memoise_2.0.1           spatstat.explore_3.2-6  RCurl_1.98-1.14         htmltools_0.5.7        
# [15] sctransform_0.4.1       parallelly_1.37.0       KernSmooth_2.23-22      htmlwidgets_1.6.4       ica_1.0-3               plyr_1.8.9              plotly_4.10.4          
# [22] zoo_1.8-12              cachem_1.0.8            igraph_2.0.2            mime_0.12               lifecycle_1.0.4         pkgconfig_2.0.3         Matrix_1.6-5           
# [29] R6_2.5.1                fastmap_1.1.1           GenomeInfoDbData_1.2.11 fitdistrplus_1.1-11     future_1.33.1           shiny_1.8.0             digest_0.6.34          
# [36] colorspace_2.1-0        patchwork_1.2.0         AnnotationDbi_1.64.1    S4Vectors_0.40.2        tensor_1.5              RSpectra_0.16-1         irlba_2.3.5.1          
# [43] textshaping_0.3.7       RSQLite_2.3.6           labeling_0.4.3          progressr_0.14.0        fansi_1.0.6             spatstat.sparse_3.0-3   httr_1.4.7             
# [50] polyclip_1.10-6         abind_1.4-5             compiler_4.3.2          bit64_4.0.5             withr_3.0.0             DBI_1.2.2               fastDummies_1.7.3      
# [57] MASS_7.3-60             tools_4.3.2             lmtest_0.9-40           httpuv_1.6.14           future.apply_1.11.1     goftest_1.2-3           glue_1.7.0             
# [64] nlme_3.1-163            promises_1.2.1          grid_4.3.2              Rtsne_0.17              cluster_2.1.4           reshape2_1.4.4          generics_0.1.3         
# [71] gtable_0.3.4            spatstat.data_3.0-4     tidyr_1.3.1             data.table_1.15.0       XVector_0.42.0          utf8_1.2.4              BiocGenerics_0.48.1    
# [78] spatstat.geom_3.2-8     RcppAnnoy_0.0.22        ggrepel_0.9.5           RANN_2.6.1              pillar_1.9.0            stringr_1.5.1           spam_2.10-0            
# [85] RcppHNSW_0.6.0          later_1.3.2             splines_4.3.2           lattice_0.21-9          survival_3.5-7          bit_4.0.5               deldir_2.0-2           
# [92] tidyselect_1.2.1        Biostrings_2.70.3       miniUI_0.1.1.1          pbapply_1.7-2           gridExtra_2.3           IRanges_2.36.0          scattermore_1.2        
# [99] stats4_4.3.2            Biobase_2.62.0          matrixStats_1.2.0       stringi_1.8.3           lazyeval_0.2.2          codetools_0.2-19        tibble_3.2.1           
# [106] cli_3.6.2               uwot_0.1.16             systemfonts_1.0.6       xtable_1.8-4            reticulate_1.35.0       munsell_0.5.0           GenomeInfoDb_1.38.8    
# [113] Rcpp_1.0.12             globals_0.16.2          spatstat.random_3.2-2   png_0.1-8               parallel_4.3.2          ellipsis_0.3.2          blob_1.2.4             
# [120] dotCall64_1.1-1         bitops_1.0-7            listenv_0.9.1           viridisLite_0.4.2       scales_1.3.0            ggridges_0.5.6          crayon_1.5.2           
# [127] leiden_0.4.3.1          packrat_0.9.2           purrr_1.0.2             rlang_1.1.3             KEGGREST_1.42.0         cowplot_1.1.3          
