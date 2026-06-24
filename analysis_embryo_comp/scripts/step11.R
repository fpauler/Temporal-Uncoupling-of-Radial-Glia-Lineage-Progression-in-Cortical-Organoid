# ====================================================================================================
# step 11
# script to analyse age matching of EBs based on reference data
# done by comparing transcriptomes between key cell types of reference and our data
# ====================================================================================================

library (ggplot2)
library (Seurat)
library (pheatmap)

set.seed( 2401 )

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

# =======================
# define functions here
# =======================
get_pseudobulk <- function(sobj, by1 = "celltype", by2 = "stage"){
  expr <- GetAssayData(sobj, layer = "data")   # normalized data

  meta <- sobj@meta.data[, c(by1, by2)]
  meta[,1] <- as.character( meta[,1] )
  meta[,2] <- as.character( meta[,2] )
  
  combos <- unique(meta)
  
  pb <- lapply(1:nrow(combos), function(i){
    
    cells <- rownames(meta[which(meta[ by1 ] == combos[i,by1] & meta[ by2 ] == combos[i,by2]),])
    
    if(length(cells) > 0){
      rowMeans(expr[, cells, drop = FALSE])
    } else {
      NULL
    }
  })
  pb <- do.call(cbind, pb)
  colnames(pb) <- paste(combos[[by1]], combos[[by2]], sep = "_")
  return(pb)
}

# function based on pheatmap 
scale_row <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  return((x - m) / s)
}

# our data
seurat.org <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.downsample.RDS", sep=""))

# extract organoid data
Idents( seurat.org ) <- "group"
seurat.org <- subset( seurat.org, idents = "EB")

# focus on key cell types
Idents( seurat.org ) <- "cellType"
seurat.org <- subset( seurat.org, idents = c("RGP", "IP", "iN", "astro", "oligo"))

seurat.org$cellType[which(seurat.org$cellType == "RGP")] <- "Radial glia"
seurat.org$cellType[which(seurat.org$cellType == "IP")] <- "Intermediate Prog"
seurat.org$cellType[which(seurat.org$cellType == "iN")] <- "Neuron"
seurat.org$cellType[which(seurat.org$cellType %in% c("astro", "oligo"))] <- "Glia"

pb_org <- get_pseudobulk( sobj = seurat.org, by1 = "cellType", by2 = "age")

# reference data from LaManno was analysed on ISTA cluster
# all files and descriptions in other_data folder
pb_ref <- readRDS ( file = paste(nobackup_base, "RDS_files/LaManno_pseudobulk.rds", sep=""))

# match cell type labels between studies and figure
colnames( pb_ref ) <- gsub( replacement = "Intermediate Prog", pattern = "Neuroblast",  x = colnames( pb_ref ) )
common_genes <- intersect(rownames(pb_org), rownames(pb_ref))
pb_org <- pb_org[common_genes,]
pb_ref <- pb_ref[common_genes,]

age_ref_list <- unique( sapply( colnames(pb_ref), function (x) strsplit(x, "_")[[1]][2] ) )
age_ref_list <- gsub( pattern = "e9.0", replacement = "e09.0", x = age_ref_list)
age_ref_list <- sort( as.character(age_ref_list) )
age_ref_list <- gsub( pattern = "e09.0", replacement = "e9.0", x = age_ref_list)

comp_mat_list <- list()
for (cell_type in c("Radial glia", "Neuron")) {
  
  sapply( age_ref_list, function (age_ref) {
    sapply( c("D8", "D13", "D20", "D25"), function (age_org) {
      
      ref_sample_id <- paste(cell_type, age_ref, sep="_")
      org_sample_id <- paste(cell_type, age_org, sep="_")
      
      if ( ref_sample_id %in% colnames(pb_ref) & org_sample_id %in% colnames(pb_org) ){
        tmp_ref <- pb_ref[ , ref_sample_id]
        tmp_org <- pb_org[ , org_sample_id]
        message ( age_ref, " ", age_org, " ", identical(names(tmp_ref), names(tmp_org) ) )
        if(identical(names(tmp_ref), names(tmp_org) )) {
          cor( tmp_ref, tmp_org )
        } else {
          return(NA)
        }
        
      } else {
        return (NA)
      }
    })
  }) -> sim_mat
  
  # scaling is done as in pheatmap package
  if (cell_type == "Radial glia") {
    scale_mat <- t(scale_row(t(sim_mat[c("D8", "D13", "D20"),])))
  }
  if (cell_type == "Neuron") {
    scale_mat <- t(scale_row(t(sim_mat[c("D13", "D20", "D25"),])))
  }
  if (cell_type == "Intermediate Prog") {
    scale_mat <- t(scale_row(t(sim_mat[c("D13", "D20", "D25"),])))
  }
  if (cell_type == "Glia") {
    scale_mat <- t(scale_row(t(sim_mat[c("D20", "D25"),])))
  }
  # plot differences - 0 means no difference
  # in other words: the higher the better, 0 is the max
  ct_out <- gsub( pattern = " ", replacement = "_", x = cell_type)
  comp_mat <- apply( scale_mat, 1, function (x) x - max(x, na.rm=T) )
  
  pheatmap(comp_mat, 
           scale="none", cluster_cols = F, cluster_rows = F,
           main = cell_type, filename = paste(base, "/plots/Sup_Fig_4/cellType_sim_", ct_out, ".pdf", sep=""))
  
  comp_mat_list[[ct_out]] <- comp_mat
}

# determine the best matching age
sapply(comp_mat_list, function (ct){
  apply(ct, 2, function(x) rownames(ct)[which(x == max(x, na.rm = T))])
})

#     Radial_glia Neuron
# D8  "e10"       "e10" 
# D13 "e12"       "e16" 
# D20 "e18"       "e18" 

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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12    Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-3           ggplot2_3.5.2     
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-2           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.3           
# [5] magrittr_2.0.3         RcppAnnoy_0.0.22       spatstat.geom_3.2-8    matrixStats_1.2.0     
# [9] ggridges_0.5.6         compiler_4.3.2         png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.1.1         
# [17] ellipsis_0.3.2         utf8_1.2.4             promises_1.2.1         purrr_1.0.2           
# [21] jsonlite_1.8.8         goftest_1.2-3          later_1.3.2            spatstat.utils_3.1-0  
# [25] irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.4          R6_2.5.1              
# [29] ica_1.0-3              stringi_1.8.3          RColorBrewer_1.1-3     spatstat.data_3.0-4   
# [33] reticulate_1.44.1      parallelly_1.37.0      lmtest_0.9-40          scattermore_1.2       
# [37] Rcpp_1.0.12            tensor_1.5             future.apply_1.11.1    zoo_1.8-12            
# [41] sctransform_0.4.1      httpuv_1.6.14          Matrix_1.6-5           splines_4.3.2         
# [45] igraph_2.1.4           tidyselect_1.2.1       rstudioapi_0.16.0      abind_1.4-5           
# [49] spatstat.random_3.2-2  codetools_0.2-19       miniUI_0.1.1.1         spatstat.explore_3.2-6
# [53] listenv_0.9.1          lattice_0.21-9         tibble_3.2.1           plyr_1.8.9            
# [57] shiny_1.8.0            withr_3.0.0            ROCR_1.0-11            Rtsne_0.17            
# [61] future_1.33.1          fastDummies_1.7.3      survival_3.5-7         polyclip_1.10-6       
# [65] fitdistrplus_1.1-11    pillar_1.9.0           packrat_0.9.2          KernSmooth_2.23-22    
# [69] plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0         scales_1.4.0          
# [73] globals_0.16.2         xtable_1.8-4           glue_1.7.0             lazyeval_0.2.2        
# [77] tools_4.3.2            data.table_1.15.0      RSpectra_0.16-1        RANN_2.6.1            
# [81] leiden_0.4.3.1         dotCall64_1.1-1        cowplot_1.1.3          grid_4.3.2            
# [85] tidyr_1.3.1            colorspace_2.1-0       nlme_3.1-163           patchwork_1.2.0       
# [89] cli_3.6.2              spatstat.sparse_3.0-3  spam_2.10-0            fansi_1.0.6           
# [93] viridisLite_0.4.2      dplyr_1.1.4            uwot_0.1.16            gtable_0.3.6          
# [97] digest_0.6.34          progressr_0.14.0       ggrepel_0.9.5          htmlwidgets_1.6.4     
# [101] farver_2.1.1           htmltools_0.5.7        lifecycle_1.0.4        httr_1.4.7            
# [105] mime_0.12              MASS_7.3-60
