# =====================================================
# step 17
# determine developmental trajectories using slingshot
# =====================================================

library (Seurat)
library (slingshot)
library (Matrix)
library (ggplot2)
library (patchwork)
library (dplyr)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

EB_neurons <- readRDS( file = paste(nobackup_base, "/RDS_files/EB_neurons.RDS", sep=""))

Idents( EB_neurons ) <- "broad_group"
organoid_neurons <- subset( EB_neurons, cells = rownames(EB_neurons@meta.data[which(!is.na(EB_neurons$broad_group)),]))
organoid_neurons@meta.data$broad_group[ which(organoid_neurons@meta.data$broad_group %in% c("UL", "DL"))] <- "Neurons"
organoid_neurons@meta.data$broad_group <- factor( organoid_neurons@meta.data$broad_group, levels = c("AP", "IP", "migrating", "Neurons"))

# prepare the basic clustering of the data
# slingshot is cluster based and I combine all neuronal clusters to 1
# therefore I need to automatically annotate clusters when testing different parameters
# here I define a base clustering, that is used as a reference
EB_neurons <- FindNeighbors(EB_neurons, dims = 1:ncol(Embeddings(EB_neurons, "css")), verbose = FALSE, reduction = "css")
EB_neurons <- FindClusters(EB_neurons, resolution = 1, verbose = FALSE)

EB_neurons@meta.data$base_clustering <- EB_neurons@meta.data$seurat_clusters 

# note that this clustering comes from the initial integrated data, this is where the seemingly misplaced cells come from
DimPlot( EB_neurons, group.by = "seurat_clusters", label=T)
ggsave( filename = paste(base, "/plots/Fig_4/Neuron_traj/Dev_Traj_initial_UMAP_clusters.png", sep="") )

i <- 1
for (n_neighbors in seq (10, 30, 10)) {
  for (min_dist in seq(0.1, 0.3, 0.1)) {
    
    EB_neurons <- RunUMAP(EB_neurons, reduction = "css", dims = 1:ncol(Embeddings(EB_neurons, "css")), 
                          n.neighbors = n_neighbors, metric = "euclidean", 
                          min.dist = min_dist, return.model = T, seed.use = 2401)
    
    for (res in c(1,2,4)) {
      EB_neurons <- FindClusters(EB_neurons, resolution = res)
      
      # combine all neurons into one cluster
      # do this via the initial clustering, where I identified the clusters associated with neurons
      age_assoc <- sapply( unique( as.character(EB_neurons@meta.data$seurat_clusters) ), function (cluster){
        freq <- table( EB_neurons@meta.data$base_clustering[ which(EB_neurons@meta.data$seurat_clusters == cluster)] )
        rel_freq <- freq / sum(freq)
        if (any(!is.na(rel_freq))) {
          tmp_names <- names( which( rel_freq == max( rel_freq, na.rm = T ) ) )
          # catch exception where 2 clusters are equally abundant
          if (length(tmp_names) == 1) {
            return(tmp_names)
          } else {
            return(NA)
          }
          
        } else {
          return (NA)
        }
      })
      
      age_assoc <- age_assoc[ which (!is.na(age_assoc)) ]
      neuron_clusters <- names(age_assoc)[ which( age_assoc %in% as.character(c(8,5,15,13,12,6,1,2))) ] 
      
      EB_neurons@meta.data$seurat_clusters <- as.character ( EB_neurons@meta.data$seurat_clusters )
      EB_neurons@meta.data$seurat_clusters[ which( as.character(EB_neurons@meta.data$seurat_clusters) %in% neuron_clusters ) ] <- "99"
        
      umap_plot <- DimPlot( EB_neurons, group.by = "seurat_clusters", label = T) + NoLegend()
      ggsave( filename = paste(base, "/plots/Fig_4/Neuron_traj/EB_neuron_umap_", i, ".", res, ".", n_neighbors, ".", min_dist, ".pdf", sep=""), plot = umap_plot)
      
      meta_filter <- EB_neurons@meta.data
      meta_filter <- cbind(meta_filter, Embeddings(EB_neurons, "umap"))
  
      mm <- sparse.model.matrix(~ 0 + factor(meta_filter$cellType))
      colnames(mm) <- levels(factor(meta_filter$cellType))
      centroids2d <- as.matrix(t(t(EB_neurons@reductions$umap@cell.embeddings[ rownames(meta_filter), ]) %*% mm) / Matrix::colSums(mm))
      
      lineages <- getLineages(
        data           = EB_neurons@reductions$umap@cell.embeddings[ rownames(meta_filter), ],
        clusterLabels  = meta_filter$seurat_clusters,
        dist.method    = "slingshot", # It can be: "simple", "scaled.full", "scaled.diag", "slingshot" or "mnn"
        omega = T,
        omega_scale = 2,
      ) 
      
      col_list <- gg_color_hue ( length(unique(EB_neurons@meta.data$cellType)) )
      names(col_list) <- as.character( unique(EB_neurons@meta.data$cellType) )
      
      pdf( paste(base, "/plots/Fig_4/Neuron_traj/EB_neuron_trajectory_", i, ".", res, ".", n_neighbors, ".", min_dist, ".pdf", sep="") )
      plot(EB_neurons@reductions$umap@cell.embeddings[rownames(meta_filter),], 
           col = col_list[ as.character(EB_neurons@meta.data$cellType)], cex = .5, pch = 16)
      lines( as.SlingshotDataSet(lineages), lwd = 1, col = "black", cex = 2)
      dev.off()
      i<- i+1
      
    }
  }
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.2.0             ggplot2_3.5.0               Matrix_1.6-5                slingshot_2.10.0           
# [5] TrajectoryUtils_1.10.1      SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0 Biobase_2.62.0             
# [9] GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0              S4Vectors_0.40.2           
# [13] BiocGenerics_0.48.1         MatrixGenerics_1.14.0       matrixStats_1.2.0           princurve_2.1.6            
# [17] Seurat_5.0.1                SeuratObject_5.0.1          sp_2.1-3                   
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8          magrittr_2.0.3         
# [5] spatstat.utils_3.1-0    farver_2.1.1            ragg_1.3.0              zlibbioc_1.48.2        
# [9] vctrs_0.6.5             ROCR_1.0-11             spatstat.explore_3.2-6  RCurl_1.98-1.14        
# [13] S4Arrays_1.2.1          htmltools_0.5.7         SparseArray_1.2.4       sctransform_0.4.1      
# [17] parallelly_1.37.0       KernSmooth_2.23-22      htmlwidgets_1.6.4       ica_1.0-3              
# [21] plyr_1.8.9              plotly_4.10.4           zoo_1.8-12              igraph_2.1.4           
# [25] mime_0.12               lifecycle_1.0.4         pkgconfig_2.0.3         R6_2.5.1               
# [29] fastmap_1.1.1           GenomeInfoDbData_1.2.11 fitdistrplus_1.1-11     future_1.33.1          
# [33] shiny_1.8.0             digest_0.6.34           colorspace_2.1-0        tensor_1.5             
# [37] RSpectra_0.16-1         irlba_2.3.5.1           textshaping_0.3.7       labeling_0.4.3         
# [41] progressr_0.14.0        fansi_1.0.6             spatstat.sparse_3.0-3   httr_1.4.7             
# [45] polyclip_1.10-6         abind_1.4-5             compiler_4.3.2          withr_3.0.0            
# [49] fastDummies_1.7.3       MASS_7.3-60             DelayedArray_0.28.0     tools_4.3.2            
# [53] lmtest_0.9-40           httpuv_1.6.14           future.apply_1.11.1     goftest_1.2-3          
# [57] glue_1.7.0              nlme_3.1-163            promises_1.2.1          grid_4.3.2             
# [61] Rtsne_0.17              cluster_2.1.4           reshape2_1.4.4          generics_0.1.3         
# [65] gtable_0.3.4            spatstat.data_3.0-4     tidyr_1.3.1             data.table_1.15.0      
# [69] utf8_1.2.4              XVector_0.42.0          spatstat.geom_3.2-8     RcppAnnoy_0.0.22       
# [73] ggrepel_0.9.5           RANN_2.6.1              pillar_1.9.0            stringr_1.5.1          
# [77] spam_2.10-0             RcppHNSW_0.6.0          later_1.3.2             splines_4.3.2          
# [81] dplyr_1.1.4             lattice_0.21-9          survival_3.5-7          deldir_2.0-2           
# [85] tidyselect_1.2.1        miniUI_0.1.1.1          pbapply_1.7-2           gridExtra_2.3          
# [89] scattermore_1.2         stringi_1.8.3           lazyeval_0.2.2          codetools_0.2-19       
# [93] tibble_3.2.1            cli_3.6.2               uwot_0.1.16             systemfonts_1.0.6      
# [97] xtable_1.8-4            reticulate_1.35.0       munsell_0.5.0           Rcpp_1.0.12            
# [101] globals_0.16.2          spatstat.random_3.2-2   png_0.1-8               parallel_4.3.2         
# [105] ellipsis_0.3.2          dotCall64_1.1-1         bitops_1.0-7            listenv_0.9.1          
# [109] viridisLite_0.4.2       scales_1.3.0            ggridges_0.5.6          crayon_1.5.2           
# [113] leiden_0.4.3.1          packrat_0.9.2           purrr_1.0.2             rlang_1.1.3            
# [117] cowplot_1.1.3        
