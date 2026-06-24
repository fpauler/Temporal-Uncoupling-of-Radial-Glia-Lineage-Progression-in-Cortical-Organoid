# ==========================================
# step 19
# compare RGPs between in vitro and in vivo
# ==========================================

library (Seurat)
library (clusterProfiler)
library (org.Mm.eg.db)
library (GOSemSim)
library (simplifyEnrichment)
library (openxlsx)
library (ggplot2)
library (pheatmap)
library (ggVennDiagram)
library (cowplot)

set.seed(2401)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

seurat_sketch <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.downsample.RDS", sep=""))

Idents( seurat_sketch) <- "cellType"
DefaultAssay( seurat_sketch ) <- "RNA"

seurat_rgp <- subset( seurat_sketch, idents = "RGP")
seurat_rgp <- JoinLayers( seurat_rgp )

table(seurat_rgp$age)
#   D8  D13  D20  D25  E10  E13  E16   P0 
# 1797 3517   95    3 2473 3182   69    1 

DEG_stats <- list()
# we have 2 clearly distinct RGP groups: D8 and D23, separately 
Idents( seurat_rgp) <- "age"
DEG_stats[["D8/E10"]] <- FindMarkers( seurat_rgp, ident.1 = "D8", ident.2 = "E10")
DEG_stats[["D8/E10"]]$SYMBOL <- rownames( DEG_stats[["D8/E10"]] )
DEG_stats[["D13/E13"]] <- FindMarkers( seurat_rgp, ident.1 = "D13", ident.2 = "E13")
DEG_stats[["D13/E13"]]$SYMBOL <- rownames( DEG_stats[["D13/E13"]] )

DEG_list <- list()
DEG_list[["D08 org spec"]] <- DEG_stats[["D8/E10"]][ which(DEG_stats[["D8/E10"]]$p_val_adj < 0.01 & DEG_stats[["D8/E10"]]$avg_log2FC > 1), "SYMBOL"]
DEG_list[["D13 org spec"]] <- DEG_stats[["D13/E13"]][ which(DEG_stats[["D13/E13"]]$p_val_adj < 0.01 & DEG_stats[["D13/E13"]]$avg_log2FC > 1), "SYMBOL"]
DEG_list[["E10 embryo spec"]] <- DEG_stats[["D8/E10"]][ which(DEG_stats[["D8/E10"]]$p_val_adj < 0.01 & DEG_stats[["D8/E10"]]$avg_log2FC < -1), "SYMBOL"]
DEG_list[["E13 embryo spec"]] <- DEG_stats[["D13/E13"]][ which(DEG_stats[["D13/E13"]]$p_val_adj < 0.01 & DEG_stats[["D13/E13"]]$avg_log2FC < -1), "SYMBOL"]

ggVennDiagram( DEG_list ) + scale_fill_gradient(low="grey90",high = "red")
ggsave( filename = paste(base, "plots/Sup_Fig_24/S24c_DEG_overlap.pdf", sep=""))

# plot the distribution of fold-changes to illustrate the cutoffs used
hist_plot <- list()
for (age in c("D8/E10", "D13/E13")) {
  hist_df <- DEG_stats[[ age ]][ which( DEG_stats[[ age ]]$p_val_adj < 0.01),]
  hist_df$cat <- "ns"
  hist_df$cat[which(hist_df$avg_log2FC > 1)] <- "org_high"
  hist_df$cat[which(hist_df$avg_log2FC < -1)] <- "embryo_high"
  
  hist_plot[[ age ]] <- ggplot( hist_df, aes(x=avg_log2FC, fill=cat)) + 
    geom_histogram( breaks = seq(-15, 15, 1)) +
    theme_classic() + scale_fill_manual( values = c("org_high" = "goldenrod", "embryo_high" = "cornflowerblue", "ns" = "grey50")) +
    ggtitle(age, "log2 fold change of DEGs with adjusted p < 0.01") + ylim(0,2100) + NoLegend()
  
}

plot_grid( plotlist = hist_plot )
ggsave( filename = paste(base, "plots/Sup_Fig_24/DEG_lFC_hist.pdf", sep=""), width = 8)

# determine actual DEG numbers with different cut offs:
deg_n_df <- data.frame()
for (age in c("D8/E10", "D13/E13")) {
  for (cutoff in c(0,1)) {
    for (dir in c("org spec", "embryo spec")) {
      if (dir == "org spec") {
        n <- length( which( DEG_stats[[ age ]]$p_val_adj < 0.01 & DEG_stats[[ age ]]$avg_log2FC > cutoff ) )
      } else {
        n <- length( which( DEG_stats[[ age ]]$p_val_adj < 0.01 & DEG_stats[[ age ]]$avg_log2FC < cutoff*-1 ) )
      }
      deg_n_df <- rbind( deg_n_df, data.frame( age = age, dir = dir, cutoff = cutoff, n = n) )
    }
  }
}

deg_n_df$age <- factor( deg_n_df$age, levels = c("D8/E10", "D13/E13"))

ggplot( deg_n_df[which(deg_n_df$cutoff == 1),], aes(x=dir, y=n, fill=dir)) + 
  geom_bar(stat = "identity", position = "dodge", color="black") +
  theme_classic() + scale_fill_manual( values = c("org spec" = "yellow", "embryo spec" = "blue", "ns" = "grey50")) +
  facet_wrap(~age) +
  ggtitle("number of DEGs with adjusted p < 0.01") + ylim(0,2100) + NoLegend()
ggsave( filename = paste(base, "plots/Sup_Fig_24/S24b_DEG_number.pdf", sep=""), width = 4)

eGO <- list() 
for ( name in names(DEG_list)) {
  gene_conv <- bitr( geneID = DEG_list[[ name ]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) 
  for ( ont in c("BP")) {
    eGO[[ paste(name, ont, sep="_") ]] <- enrichGO( gene = gene_conv$ENTREZID, OrgDb = org.Mm.eg.db, ont = ont, readable = T)
  }
}

list_of_plots <- list()
for (name in names(eGO)) {
  name_parts <- strsplit( name, " ")[[1]]
  tmp_df <- eGO[[ name ]]@result[1:20, c("Description", "p.adjust")]
  tmp_df$age <- name_parts[1]
  tmp_df$type <- name_parts[2]
  tmp_df$score <- log10(tmp_df$p.adjust) * -1
  tmp_df$Description <- factor(tmp_df$Description, levels = rev(tmp_df$Description))
  list_of_plots[[ name_parts[1] ]] <- ggplot(tmp_df, aes(x=score, y=Description, fill = type)) + 
    geom_bar(stat="identity", position ="dodge", color="black") + 
    scale_fill_manual( values = c("org" = "yellow", "embryo" = "blue")) +
    ggtitle ( name_parts[1] ) + theme_classic()
  
}

list_of_plots[[ 1 ]] + list_of_plots[[ 2 ]] + list_of_plots[[ 3 ]] + list_of_plots[[ 4 ]]
ggsave( filename = paste(base, "plots/Sup_Fig_24/S24defg_TopGO_barplot.pdf", sep=""), width=15, height=8)

wb <- createWorkbook()

for ( name in names(DEG_stats)) {
  
  sheetName <- paste("DEG", gsub(pattern = "/", replacement = " vs ", x = name))
  tmp_go_df <- DEG_stats[[ name ]]
  
  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, tmp_go_df)
  addFilter(wb, sheetName, row = 1, cols = 1:ncol( tmp_go_df ))
  setColWidths(wb, sheetName, cols = 1:ncol( tmp_go_df ), widths=10)
}

for ( name in names(eGO)) {
  sheetName <- paste("GO", strsplit(name, " ")[[1]][1])
  tmp_go_df <- eGO[[ name ]]@result
  
  addWorksheet(wb, sheetName)
  writeData(wb, sheetName, tmp_go_df)
  addFilter(wb, sheetName, row = 1, cols = 1:ncol( tmp_go_df ))
  setColWidths(wb, sheetName, cols = 1:ncol( tmp_go_df ), widths=10)
}

saveWorkbook(wb, file = paste(base, "plots/Sup_Fig_24/DEG_GO.xlsx", sep=""), overwrite = T) 

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
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cowplot_1.1.3             ggVennDiagram_1.5.2       pheatmap_1.0.12          
# [4] ggplot2_3.5.0             openxlsx_4.2.5.2          simplifyEnrichment_1.12.0
# [7] GOSemSim_2.28.1           org.Mm.eg.db_3.18.0       AnnotationDbi_1.64.1     
# [10] IRanges_2.36.0            S4Vectors_0.40.2          Biobase_2.62.0           
# [13] BiocGenerics_0.48.1       clusterProfiler_4.10.1    Seurat_5.0.1             
# [16] SeuratObject_5.0.1        sp_2.1-3                 
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.3.2           later_1.3.2             bitops_1.0-7           
# [5] ggplotify_0.1.2         tibble_3.2.1            polyclip_1.10-6         fastDummies_1.7.3      
# [9] lifecycle_1.0.4         doParallel_1.0.17       NLP_0.3-0               globals_0.16.2         
# [13] lattice_0.21-9          MASS_7.3-60             magrittr_2.0.3          limma_3.58.1           
# [17] plotly_4.10.4           httpuv_1.6.14           sctransform_0.4.1       zip_2.3.1              
# [21] spam_2.10-0             spatstat.sparse_3.0-3   reticulate_1.44.1       pbapply_1.7-2          
# [25] DBI_1.2.2               RColorBrewer_1.1-3      abind_1.4-5             zlibbioc_1.48.2        
# [29] Rtsne_0.17              presto_1.0.0            purrr_1.0.2             ggraph_2.2.1           
# [33] RCurl_1.98-1.14         yulab.utils_0.1.4       tweenr_2.0.3            circlize_0.4.16        
# [37] GenomeInfoDbData_1.2.11 enrichplot_1.22.0       tm_0.7-14               ggrepel_0.9.5          
# [41] irlba_2.3.5.1           listenv_0.9.1           spatstat.utils_3.1-0    tidytree_0.4.6         
# [45] goftest_1.2-3           RSpectra_0.16-1         spatstat.random_3.2-2   fitdistrplus_1.1-11    
# [49] parallelly_1.37.0       leiden_0.4.3.1          codetools_0.2-19        xml2_1.3.6             
# [53] DOSE_3.28.2             ggforce_0.4.2           shape_1.4.6.1           tidyselect_1.2.1       
# [57] aplot_0.2.2             farver_2.1.1            viridis_0.6.5           matrixStats_1.2.0      
# [61] spatstat.explore_3.2-6  jsonlite_1.8.8          GetoptLong_1.0.5        ellipsis_0.3.2         
# [65] tidygraph_1.3.1         progressr_0.14.0        iterators_1.0.14        ggridges_0.5.6         
# [69] survival_3.5-7          foreach_1.5.2           tools_4.3.2             treeio_1.26.0          
# [73] ica_1.0-3               Rcpp_1.0.12             glue_1.7.0              gridExtra_2.3          
# [77] qvalue_2.34.0           GenomeInfoDb_1.38.8     dplyr_1.1.4             withr_3.0.0            
# [81] fastmap_1.1.1           fansi_1.0.6             digest_0.6.34           gridGraphics_0.5-1     
# [85] R6_2.5.1                mime_0.12               colorspace_2.1-0        scattermore_1.2        
# [89] GO.db_3.18.0            tensor_1.5              spatstat.data_3.0-4     RSQLite_2.3.6          
# [93] utf8_1.2.4              tidyr_1.3.1             generics_0.1.3          data.table_1.15.0      
# [97] graphlayouts_1.1.1      httr_1.4.7              htmlwidgets_1.6.4       scatterpie_0.2.2       
# [101] uwot_0.1.16             pkgconfig_2.0.3         gtable_0.3.4            blob_1.2.4             
# [105] ComplexHeatmap_2.18.0   lmtest_0.9-40           XVector_0.42.0          shadowtext_0.1.3       
# [109] htmltools_0.5.7         dotCall64_1.1-1         fgsea_1.28.0            clue_0.3-65            
# [113] scales_1.3.0            png_0.1-8               ggfun_0.1.4             rstudioapi_0.16.0      
# [117] rjson_0.2.23            reshape2_1.4.4          nlme_3.1-163            org.Hs.eg.db_3.18.0    
# [121] GlobalOptions_0.1.2     cachem_1.0.8            zoo_1.8-12              stringr_1.5.1          
# [125] KernSmooth_2.23-22      parallel_4.3.2          miniUI_0.1.1.1          HDO.db_0.99.1          
# [129] proxyC_0.4.1            pillar_1.9.0            vctrs_0.6.5             RANN_2.6.1             
# [133] slam_0.1-50             promises_1.2.1          xtable_1.8-4            cluster_2.1.4          
# [137] packrat_0.9.2           cli_3.6.2               compiler_4.3.2          rlang_1.1.3            
# [141] crayon_1.5.2            future.apply_1.11.1     labeling_0.4.3          plyr_1.8.9             
# [145] fs_1.6.3                stringi_1.8.3           viridisLite_0.4.2       deldir_2.0-2           
# [149] BiocParallel_1.36.0     munsell_0.5.0           Biostrings_2.70.3       lazyeval_0.2.2         
# [153] spatstat.geom_3.2-8     Matrix_1.6-5            RcppHNSW_0.6.0          patchwork_1.2.0        
# [157] bit64_4.0.5             future_1.33.1           statmod_1.5.0           KEGGREST_1.42.0        
# [161] shiny_1.8.0             ROCR_1.0-11             igraph_2.1.4            memoise_2.0.1          
# [165] ggtree_3.10.1           fastmatch_1.1-4         bit_4.0.5               gson_0.1.0             
# [169] ape_5.8   
