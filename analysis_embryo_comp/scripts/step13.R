# =======================================================================================
# step 13
# prepare PatchSeq data - takes input from gene expression analysis stored in RDS_files
# details on alignment/gene expression analysis in ../MADM-CloneSeq/alignment_details
# =======================================================================================

library (Seurat)
library (openxlsx)
library (dplyr)
library (ggplot2)
library (pheatmap)
library (ggbeeswarm)

# this is based on the scripts from Cheung et al. Neuron 2024

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "/MADM-CloneSeq/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

set.seed(2401)

###############
##
# prepare environment hash table as lookup for ENSMUSG - Symbol and ENTREZID
##
###############

ensmusg_symbol <- new.env()
  
#read conversion table
ensmusg_symbol_chr <- read.table(paste(base, "annotation/exchange_table_M27.tsv", sep="/"), fill = T, header = F, stringsAsFactors = F)
colnames(ensmusg_symbol_chr) <- c("ENSMUSG", "SYMBOL", "chr")
  
# sanity check for unique ENSMUSG IDs
table(table(ensmusg_symbol_chr$ENSMUSG))
  
# 1 
# 55359 
  
#fill conversion table - handled just like a list but apparently much faster
for (x in 1:nrow(ensmusg_symbol_chr)) {
  #deal here with multi mappers - one SYMBOL mapped to many ENSMUSG
  #a comma in the SYMBOL name indicates these multi-mappers
  ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]] <- ensmusg_symbol_chr$SYMBOL[x]
}

###############
##
# prepare PatchSeq dataset
##
###############

# read exon and intron counts from step 1
exon_counts <- readRDS( paste(base, "/RDS_files/raw_exon_counts.RDS", sep="") )
intron_counts <- readRDS( paste(base, "/RDS_files/raw_intron_counts.RDS", sep="") )

# sanity check
identical( rownames(intron_counts), rownames(exon_counts) )
# [1] TRUE
identical( colnames(intron_counts), colnames(exon_counts) )
# [1] TRUE 
  
#give unique symbol IDs to TPM count matrix genes
counts.matrix <- exon_counts + intron_counts
  
gn <- sapply(rownames(counts.matrix), function (x) ensmusg_symbol[[x]])
  
length(gn)
# 55359
  
#remove duplicated symbols
gn <- gn[!(duplicated(gn) | duplicated(gn, fromLast = TRUE))]
length(gn)
# 52266
  
counts.matrix <- counts.matrix[names(gn),]
rownames(counts.matrix) <- as.character(gn[rownames(counts.matrix)])

col_names <- colnames(counts.matrix)
col_names <- gsub(pattern = "STAR\\.", replacement = "", x = col_names)
col_names <- gsub(pattern = "\\.Aligned\\.sortedByCoord\\.out\\.bam", replacement = "", x = col_names)
colnames(counts.matrix) <- col_names

# read alignment data
meta.data <- readRDS( file = paste(base, "RDS_files/PatchSeq_meta_data.RDS", sep=""))
meta.data$total_reads <- as.numeric( meta.data$total_reads )
meta.data$uniquely_aligned <- as.numeric( meta.data$uniquely_aligned )
meta.data$STAR.perc <- 100 * ( meta.data$uniquely_aligned / meta.data$total_reads)
meta.data$sampleID <- rownames( meta.data )

# read sample meta data
# these are all patched samples
sample_meta <- read.xlsx( xlsxFile = paste(base, "other_data/organoids_metadata.xlsx", sep=""))
# merge meta data with sample list
meta.data <- merge( meta.data, sample_meta, by.x= "sampleID", by.y = "sample_id")
meta.data$total_reads[ which( is.na(meta.data$total_reads) ) ] <- 1
meta.data$uniquely_aligned[ which( is.na(meta.data$uniquely_aligned) ) ] <- 1
meta.data$STAR.perc[ which( is.na(meta.data$STAR.perc) ) ] <- 1

meta.data <- meta.data[which(!is.na(meta.data$sampleID)),]

rownames(meta.data) <- meta.data$sampleID
meta.data$clone_type <- gsub( pattern = " ", replacement = "", x = meta.data$clone_type )
# some QC plots

plot1 <- ggplot( meta.data, aes(x=as.factor(clone_type), y=total_reads )) + 
  geom_boxplot( outlier.size = 0.5 ) + geom_hline( yintercept = 2000000) + 
  # geom_beeswarm( size = 0.5 ) + scale_y_log10() + theme_classic() + 
  ggtitle("total reads per clone type") + ylab("log10 # total reads") + theme_classic() + 
  NoLegend()

plot2 <- ggplot( meta.data, aes(x=as.factor(plate), y=total_reads)) + 
  geom_boxplot( outlier.size = 0.5 ) + geom_hline( yintercept = 2000000) + 
  geom_beeswarm( size = 0.5 ) + scale_y_log10() + theme_classic() + 
  ggtitle("total reads per index plate") + ylab("log10 # total reads") + xlab( "index plate")

plot3 <- ggplot( meta.data, aes(x=clone_type, y=STAR.perc)) + 
  geom_boxplot( outlier.size = 0.5 ) + geom_beeswarm( size = 0.5 ) + theme_classic() + 
  ggtitle("% aligned reads per clone type") + ylab("% aligned reads") + 
  ylim(0,100)

plot4 <- ggplot( meta.data, aes(x=as.factor(plate), y=STAR.perc)) + 
  geom_boxplot( outlier.size = 0.5 ) + geom_beeswarm( size = 0.5 ) + theme_classic() + 
  ggtitle("% aligned reads per index plate") + 
  ylab("% aligned reads") + ylim(0,100) + xlab( "index plate")

comb_plot <- (plot1 + plot2) / (plot3 + plot4)
ggsave( plot = comb_plot, filename = paste(base, "QC/alignment_quality_plot01_", nrow(meta.data), "_cells.pdf", sep=""), width = 8, height=8)

# test for significance of differences
# clone_type seems to have no influence on total read number
test_aov <- aov(meta.data$total_reads ~ factor(meta.data$clone_type))
summary(test_aov)

#                               Df    Sum Sq   Mean Sq F value Pr(>F)
# factor(meta.data$clone_type)   2 9.205e+13 4.602e+13   1.944  0.145
# Residuals                    284 6.722e+15 2.367e+13  

# plate seems to have influence on total read number
test_aov <- aov(meta.data$total_reads ~ factor(meta.data$plate))
summary(test_aov)

#                          Df    Sum Sq   Mean Sq F value Pr(>F)  
# factor(meta.data$plate)   2 2.023e+14 1.012e+14   4.345 0.0138 *
#   Residuals               284 6.612e+15 2.328e+13                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(test_aov)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = meta.data$total_reads ~ factor(meta.data$plate))
# 
# $`factor(meta.data$plate)`
# diff        lwr       upr     p adj
# 2-1   968671.1  -672187.7 2609530.0 0.3470168
# 3-1 -1088741.3 -2733912.6  556429.9 0.2651516
# 3-2 -2057412.5 -3702583.7 -412241.2 0.0097183

meta.data %>%
  group_by( plate) %>%
  summarise( n=n(), mean=mean(total_reads), sd=sd(total_reads) )

# # A tibble: 3 × 4
#  plate     n     mean       sd
#  <dbl> <int>    <dbl>    <dbl>
#      1    96 4350031. 4535000.
#      2    96 5318703. 5670466.
#      3    95 3261290. 4130144.

# test for difference in % uniquely aligned reads
# arcsine transformation performed for % values
test_aov <- aov( asin(sqrt(meta.data$uniquely_aligned / meta.data$total_reads)) ~ factor(meta.data$clone_type))
summary(test_aov)

#                               Df Sum Sq Mean Sq F value Pr(>F)
# factor(meta.data$clone_type)   2  0.025 0.01274   0.183  0.833
# Residuals                    284 19.816 0.06977

# alignment rate is strongly associated with plate (batch of preparation) 
test_aov <- aov( asin(sqrt(meta.data$uniquely_aligned / meta.data$total_reads)) ~ factor(meta.data$plate))
summary(test_aov)

#                          Df Sum Sq Mean Sq F value Pr(>F)    
# factor(meta.data$plate)   2  1.313  0.6565   10.06  6e-05 ***
#   Residuals               284 18.528  0.0652                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(test_aov)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = asin(sqrt(meta.data$uniquely_aligned/meta.data$total_reads)) ~ factor(meta.data$plate))
# 
# $`factor(meta.data$plate)`
#            diff        lwr         upr     p adj
# 2-1 -0.03620582 -0.1230670  0.05065535 0.5888005
# 3-1 -0.15834763 -0.2454371 -0.07125817 0.0000745
# 3-2 -0.12214181 -0.2092313 -0.03505235 0.0030749

meta.data %>%
  group_by( plate) %>%
  summarise( n=n(), mean=mean(STAR.perc), sd=sd(STAR.perc) )

# # A tibble: 3 × 4
#   plate     n  mean    sd
#   <dbl> <int> <dbl> <dbl>
# 1     1    96  78.5  14.9
# 2     2    96  75.2  20.5
# 3     3    95  64.5  29.5

nrow(meta.data)
# 287

PatchSeq_Seurat <- CreateSeuratObject(counts = counts.matrix, min.cells = 3, min.features = 300)

# An object of class Seurat 
# 35678 features across 279 samples within 1 assay 
# Active assay: RNA (35678 features, 0 variable features)
# 1 layer present: counts

# 8 cells did not meet quality threshold for being added to Seurat object

# do scaling for all genes
all.genes <- rownames(PatchSeq_Seurat)

PatchSeq_Seurat <- NormalizeData(PatchSeq_Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
PatchSeq_Seurat <- ScaleData(PatchSeq_Seurat, features = all.genes)
PatchSeq_Seurat <- AddMetaData(object = PatchSeq_Seurat, metadata = meta.data)

# add perc mitochondrial reads - another measure of cell quality
PatchSeq_Seurat[["percent.mt"]] <- PercentageFeatureSet( PatchSeq_Seurat, pattern = "^mt-")

VlnPlot( object = PatchSeq_Seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "plate" )
ggsave(paste(base, "QC/cDNA_quality_plot02_", nrow(PatchSeq_Seurat@meta.data), "_cells.pdf", sep=""), width = 5, height=4)

# test for statistical significance of differences between plates for gene based features
test_aov <- aov(PatchSeq_Seurat@meta.data$nCount_RNA ~ factor(PatchSeq_Seurat@meta.data$plate))
summary(test_aov)

#                                          Df    Sum Sq   Mean Sq F value  Pr(>F)   
# factor(PatchSeq_Seurat@meta.data$plate)   2 1.475e+14 7.374e+13   4.734 0.00952 **
# Residuals                               276 4.299e+15 1.558e+13                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD( test_aov )

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = PatchSeq_Seurat@meta.data$nCount_RNA ~ factor(PatchSeq_Seurat@meta.data$plate))
# 
# $`factor(PatchSeq_Seurat@meta.data$plate)`
# diff        lwr       upr     p adj
# 2-1   862163.6  -498045.6 2222372.8 0.2955478
# 3-1  -923316.3 -2287236.4  440603.9 0.2494094
# 3-2 -1785479.9 -3153022.2 -417937.6 0.0064983

test_aov <- aov(PatchSeq_Seurat@meta.data$nFeature_RNA ~ factor(PatchSeq_Seurat@meta.data$plate))
summary(test_aov)
#                                          Df    Sum Sq  Mean Sq F value Pr(>F)  
# factor(PatchSeq_Seurat@meta.data$plate)   2 1.945e+08 97225016   3.191 0.0427 *
# Residuals                               276 8.410e+09 30472363                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD( test_aov )

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = PatchSeq_Seurat@meta.data$nFeature_RNA ~ factor(PatchSeq_Seurat@meta.data$plate))
# 
# $`factor(PatchSeq_Seurat@meta.data$plate)`
#           diff       lwr        upr     p adj
# 2-1   124.6807 -1777.811 2027.17231 0.9869390
# 3-1 -1710.4607 -3618.143  197.22127 0.0890874
# 3-2 -1835.1414 -3747.890   77.60684 0.0631600

###############
##
# prepare reference dataset
##
################

  seurat_EB <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.RDS", sep=""))
  Idents(seurat_EB) <- "group"
  seurat_EB <- subset( seurat_EB, idents = "EB" )
  Idents(seurat_EB) <- "age"
  seurat_EB <- subset( seurat_EB, idents = c("D20", "D25") )
  
  ###############
  ##
  # NMS calculation
  ##
  ################
  
  Idents(seurat_EB) <- "cellType"
  all.markers <- FindAllMarkers( seurat_EB, only.pos = T )
  
  all.markers %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC) -> top_markers
  
  cell_types <- table( seurat_EB@meta.data$cellType )
  cell_types <- cell_types / sum( cell_types )
  
  #         aNSC        astro           CR           iN           IP         OBNB        oligo          RGP 
  # 0.1609831029 0.0748484114 0.0002748808 0.5208990217 0.0312555583 0.0922790848 0.1144797478 0.0049801924 
  
  # calculate the mean expression for each gene in each cell type
  means <- sapply( names(cell_types), function (cell_type) {
    apply(seurat_EB[["RNA"]]$data[unique(top_markers$gene),
                                   rownames(seurat_EB@meta.data[which(seurat_EB@meta.data$cellType == cell_type),])], 1, mean)
  })
  
  length(unique(top_markers$gene))
  # 1270 genes
  
  # remove genes where the median is 0
  # focus on cell types that are > 1% abundance
  med2plot <- means
  med2plot <- med2plot[which(apply(med2plot[ , names( which( cell_types > 0.01 ) ) ], 1, max) > 0), names( which( cell_types > 0.01 ) ) ]
  
  #check min number of markers per tissue
  min_n_gl <- min(table(apply(med2plot, 1, function (x) which(x == max(x)))))
  # [1] 154
  
  #determine cell type with highest expression
  marker_list <- apply(med2plot, 1, function (x) which(x == max(x)))
  # number of marker genes per cell type

  table(marker_list)
  # marker_list
  # 1   2   3   4   5   6 
  # 235 282 205 154 192 202 
  
  marker_list_df <- data.frame(gene = names(marker_list), idx = as.numeric(marker_list))
  marker_list_df$cell_type <- colnames(med2plot)[marker_list_df$idx]
  
  #cut the gene list by taking most expressed genes
  cut_gl <- lapply(unique(marker_list_df$cell_type), function (cell_type) {
    sort(med2plot[marker_list_df[which(marker_list_df$cell_type == cell_type), "gene"], cell_type], decreasing = T)[1:min_n_gl]
  })
  
  # prepare a character vector
  cut_gl <- names(unlist(cut_gl))
  
  # filter marker list 
  marker_list_df <- marker_list_df[which(marker_list_df$gene %in% cut_gl),]
  
  #determine mean marker expression in reference cell type
  ref_med_expr <- sapply( unique(marker_list_df$cell_type), function (cell_types) {
    tmp <- apply(seurat_EB[["RNA"]]$data[ marker_list_df$gene[which(marker_list_df$cell_type == cell_types)],
                                           rownames(seurat_EB@meta.data[which(seurat_EB@meta.data$cellType == cell_types),]) ], 1, mean)
    return(mean(tmp))
  })
  
  # clean up marker list to only contain genes also informative in PatchSeq
  marker_list_df <- marker_list_df[which(marker_list_df$gene %in% rownames(PatchSeq_Seurat[["RNA"]]$data)),]
  
  # expression matrix from PatchSeq with marker genes
  lPS <- as.matrix( PatchSeq_Seurat[["RNA"]]$data[marker_list_df$gene,] )
  
  # calculate mean marker expression for each cell
  mean_marker <- sapply(colnames(lPS), function (cell){
    
    sapply ( unique(marker_list_df$cell_type), function (x) {
      tmp_counts <- lPS[marker_list_df[which(marker_list_df$cell_type == x), "gene"], cell]
      # tmp_counts <- tmp_counts[which(tmp_counts > 0)]
      mean(tmp_counts)
    })
  })
  
  #now the NMS calculation
  for (cell_type in unique(marker_list_df$cell_type) ) {
    mean_marker[cell_type,] <- mean_marker[cell_type,] / ref_med_expr[cell_type]
  }
  
  # renaming vector - will be necessary for final reporting
  ren_vec <- c("NMS_Neurons", "NMS_aNSC", "NMS_Astro", "NMS_IP", "NMS_OBNB", "NMS_oligo")
  names(ren_vec) <- c("iN", "aNSC", "astro", "IP", "OBNB", "oligo")
  
  #prepare a df with all NMS scores for each cell type
  NMS_df <- data.frame(sample_id = colnames( mean_marker ), NMS = mean_marker[1,] )
  colnames(NMS_df) <- c("sample_id", rownames(mean_marker)[1])
  
  for (i in 2: nrow(mean_marker)) {
    tmp <- data.frame(sample_id = colnames( mean_marker ), NMS =  mean_marker[i,] )
    colnames(tmp) <- c("sample_id", rownames(mean_marker)[i])
    NMS_df <- merge(NMS_df, tmp, by="sample_id")
  }
  
  NMS_mat <- NMS_df[,c(2:ncol(NMS_df))]
  rownames(NMS_mat) <- NMS_df$sample_id
  colnames(NMS_mat) <- ren_vec[ colnames(NMS_mat) ]
  
  # identify the maximum NMS score from non-neuronal (aka contaminating) t-types
  NMS_cont <- apply( NMS_mat[, c("NMS_aNSC","NMS_Astro","NMS_IP","NMS_OBNB","NMS_oligo")], 1, max)
  
  ##############################################################
  ##
  # determine a cutoff for cells which are likely not neurons
  ##
  ##############################################################
  
  # first have a look at the plot of NMS scores neurons vs contaminating
  df2plot <- data.frame ( NMS_neuron = NMS_mat$NMS_Neurons, NMS_cont = as.numeric(NMS_cont), sample_id = rownames( NMS_mat ) )
  
  ggplot( df2plot, aes ( x= NMS_neuron, y = NMS_cont)) + 
    geom_point() + theme_classic() + ggtitle( "NMS plot" )
  ggsave(filename = paste(base, "/QC/NMS_neuron_vs_NMS_cont_plot.pdf", sep=""))
  
  # looks like positive correlation of the 2 scores with outliers
  # check for correlation
  cor.test( df2plot$NMS_neuron, df2plot$NMS_cont )
  
  #   Pearson's product-moment correlation
  # 
  # data:  df2plot$NMS_neuron and df2plot$NMS_cont
  # t = 8.2116, df = 277, p-value = 8.377e-15
  # alternative hypothesis: true correlation is not equal to 0
  # 95 percent confidence interval:
  #  0.3428466 0.5322410
  # sample estimates:
  #       cor 
  # 0.4424642 
  
  # To define a cutoff for these outliers investigate the NMS scores closer
  # first have a look at the heatmap of the NMS scores
  heat_clust <- pheatmap( mat = NMS_mat, scale = "row", filename = paste(base, "/QC/NMS_heatmap_", nrow(PatchSeq_Seurat@meta.data), "cells.pdf", sep=""), show_rownames = F)
  dev.off()
  # cut the tree at height 4 to get 3 clusters
  row_tree_clusters <- cutree( heat_clust$tree_row, k = 3)
  # smaller ones are the cells high in contaminating cell type markers
  table(row_tree_clusters)
  #   1   2   3 
  # 249  17  13 
  
  row_tree_clusters_df <- data.frame( clusters = row_tree_clusters)
  row_tree_clusters_df$clusters[which(row_tree_clusters_df$clusters == 1)] <- "neuron"
  row_tree_clusters_df$clusters[which(row_tree_clusters_df$clusters == 2)] <- "cont1"
  row_tree_clusters_df$clusters[which(row_tree_clusters_df$clusters == 3)] <- "cont2"
  
  pheatmap( mat = NMS_mat, 
            scale = "row", 
            annotation_row = row_tree_clusters_df,
            filename = paste(base, "/Supplement/NMS_heatmap_", nrow(PatchSeq_Seurat@meta.data), "cells_clusters.pdf", sep=""))
  
  df2plot$cluster <- row_tree_clusters[ df2plot$sample_id ]
  
  ggplot( df2plot, aes(x=NMS_neuron, y=NMS_cont, color=as.factor(cluster))) + 
    geom_point() + scale_color_manual( values = c("1" = "grey80", "2"="black", "3"="black")) + theme_classic()
  ggsave(filename = paste(base, "/QC/NMS_cutoff_plot_heatmap_clusters.pdf", sep=""))
  
  # calculate a NMS_neuron to NMS_cont ration as a score for cutoff
  df2plot$ratio <- df2plot$NMS_neuron / df2plot$NMS_cont
  rownames(df2plot) <- df2plot$sample_id
  
  # determine the 95th percentile of that cluster as a cutoff
  quantile( df2plot[ names(row_tree_clusters)[( which ( row_tree_clusters %in% c(2,3) ))], "ratio"], probs = 0.95 )
  #      95% 
  # 1.195426
  
  # define a minimum NMS cutoff for neurons - 0.3 looks like a good compromise
  quantile( df2plot[ names(row_tree_clusters)[( which ( row_tree_clusters == 1 ))], "NMS_neuron"], probs = c(0.01, 0.05, 0.1, 0.15))
  #         1%        5%       10%       15% 
  #  0.1731167 0.2616698 0.3396072 0.3824269
  
  # define the categories
  df2plot$category <- "cont"
  df2plot$category[which(df2plot$ratio > 1.2)] <- "Neuron"
  df2plot$category[which(df2plot$NMS_neuron < 0.3)] <- "cont"
  
  # plot the NMS scores for the 2 categories
  ggplot( df2plot, aes ( x= NMS_neuron, y = NMS_cont, color = category)) + 
    geom_point() + geom_smooth(method = "glm", level = 0.99) + theme_classic() + NoLegend()
  ggsave(filename = paste(base, "/QC/NMS_cutoff_plot.pdf", sep=""))
  # `geom_smooth()` using formula = 'y ~ x'
  
  # test for correlation of NMS_neuron and NMS_cont in high quality neurons - much improved
  cor.test ( df2plot[which(df2plot$ratio > 1.2), "NMS_neuron"], df2plot[which(df2plot$ratio > 1.2), "NMS_cont"])
  
  
  #   Pearson's product-moment correlation
  # 
  # data:  df2plot[which(df2plot$ratio > 1.2), "NMS_neuron"] and df2plot[which(df2plot$ratio > 1.2), "NMS_cont"]
  # t = 14.681, df = 227, p-value < 2.2e-16
  # alternative hypothesis: true correlation is not equal to 0
  # 95 percent confidence interval:
  #  0.6247605 0.7588615
  # sample estimates:
  #       cor 
  # 0.6978773
  
  PatchSeq_Seurat <- AddMetaData( object = PatchSeq_Seurat, metadata = NMS_mat)
  PatchSeq_Seurat <- AddMetaData( object = PatchSeq_Seurat, metadata = df2plot[,c("ratio", "category")])
  
  # save the complete meta data before filtering
  wb <- createWorkbook()
  
  sheetName <- "PatchSeq metadata"
  addWorksheet(wb, sheetName)
  
  writeData( wb, sheetName, PatchSeq_Seurat@meta.data )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol( PatchSeq_Seurat@meta.data ) )
  setColWidths( wb, sheetName, cols = 1:ncol( PatchSeq_Seurat@meta.data ), widths="auto" )
  
  saveWorkbook(wb, file = paste(base, "Supplement/PatchSeq_metadata.xlsx", sep=""), overwrite = T) 
  
  # extract only Neurons
  Idents( PatchSeq_Seurat ) <- "category"
  PatchSeq_Seurat <- subset(PatchSeq_Seurat, idents = "Neuron")
  nrow(PatchSeq_Seurat@meta.data)
  # 215

  ###############
  ##
  # prepare PatchSeq data
  ##
  ################
  
  PatchSeq_Seurat <- SCTransform(PatchSeq_Seurat, vars.to.regress = c("plate", "NMS_Neurons"), verbose = FALSE)
  
  PatchSeq_Seurat <- RunPCA(PatchSeq_Seurat, assay = "SCT")
  PatchSeq_Seurat <- RunUMAP(PatchSeq_Seurat, dims = 1:12)
  
  PatchSeq_Seurat <- FindNeighbors(PatchSeq_Seurat, dims = 1:12, verbose = FALSE)
  PatchSeq_Seurat <- FindClusters(PatchSeq_Seurat, resolution = 2.5, verbose = FALSE)
  
  DimPlot( PatchSeq_Seurat, group.by = c("seurat_clusters"), pt.size = 2, label=T)
  ggsave(filename = paste(base, "/QC/PatchSeq_UMAP1.", nrow(PatchSeq_Seurat@meta.data), "cells.pdf", sep=""))
  
  DimPlot( PatchSeq_Seurat, group.by = c("plate", "clone_type"), pt.size = 2)
  ggsave(filename = paste(base, "/QC/PatchSeq_UMAP1.", nrow(PatchSeq_Seurat@meta.data), "plate_cloneTypecells.pdf", sep=""), width = 8)
  
  saveRDS( object = PatchSeq_Seurat, file = paste(nobackup_base, "RDS_files/Seurat_PatchSeq.RDS", sep="/") )
  
  
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
  #   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8   
  # [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
  # [11] LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
  # 
  # time zone: Europe/Vienna
  # tzcode source: system (glibc)
  # 
  # attached base packages:
  #   [1] stats     graphics  grDevices utils     datasets  methods   base     
  # 
  # other attached packages:
  #   [1] ggbeeswarm_0.7.2   pheatmap_1.0.12    ggplot2_3.5.0      dplyr_1.1.4        openxlsx_4.2.5.2   Seurat_5.0.1       SeuratObject_5.0.1
  # [8] sp_2.1-3          
  # 
  # loaded via a namespace (and not attached):
  #   [1] RColorBrewer_1.1-3          rstudioapi_0.16.0           jsonlite_1.8.8              magrittr_2.0.3              spatstat.utils_3.1-0       
  # [6] farver_2.1.1                zlibbioc_1.48.2             ragg_1.3.0                  vctrs_0.6.5                 ROCR_1.0-11                
  # [11] DelayedMatrixStats_1.24.0   spatstat.explore_3.2-6      RCurl_1.98-1.14             S4Arrays_1.2.1              htmltools_0.5.7            
  # [16] SparseArray_1.2.4           sctransform_0.4.1           parallelly_1.37.0           KernSmooth_2.23-22          htmlwidgets_1.6.4          
  # [21] ica_1.0-3                   plyr_1.8.9                  plotly_4.10.4               zoo_1.8-12                  igraph_2.1.4               
  # [26] mime_0.12                   lifecycle_1.0.4             pkgconfig_2.0.3             Matrix_1.6-5                R6_2.5.1                   
  # [31] fastmap_1.1.1               GenomeInfoDbData_1.2.11     MatrixGenerics_1.14.0       fitdistrplus_1.1-11         future_1.33.1              
  # [36] shiny_1.8.0                 digest_0.6.34               colorspace_2.1-0            S4Vectors_0.40.2            patchwork_1.2.0            
  # [41] tensor_1.5                  RSpectra_0.16-1             irlba_2.3.5.1               GenomicRanges_1.54.1        textshaping_0.3.7          
  # [46] labeling_0.4.3              progressr_0.14.0            fansi_1.0.6                 spatstat.sparse_3.0-3       mgcv_1.9-0                 
  # [51] httr_1.4.7                  polyclip_1.10-6             abind_1.4-5                 compiler_4.3.2              withr_3.0.0                
  # [56] fastDummies_1.7.3           MASS_7.3-60                 DelayedArray_0.28.0         tools_4.3.2                 vipor_0.4.7                
  # [61] lmtest_0.9-40               beeswarm_0.4.0              zip_2.3.1                   httpuv_1.6.14               future.apply_1.11.1        
  # [66] goftest_1.2-3               glmGamPoi_1.14.3            glue_1.7.0                  nlme_3.1-163                promises_1.2.1             
  # [71] grid_4.3.2                  Rtsne_0.17                  cluster_2.1.4               reshape2_1.4.4              generics_0.1.3             
  # [76] gtable_0.3.4                spatstat.data_3.0-4         tidyr_1.3.1                 data.table_1.15.0           XVector_0.42.0             
  # [81] utf8_1.2.4                  BiocGenerics_0.48.1         spatstat.geom_3.2-8         RcppAnnoy_0.0.22            ggrepel_0.9.5              
  # [86] RANN_2.6.1                  pillar_1.9.0                stringr_1.5.1               spam_2.10-0                 RcppHNSW_0.6.0             
  # [91] limma_3.58.1                later_1.3.2                 splines_4.3.2               lattice_0.21-9              survival_3.5-7             
  # [96] deldir_2.0-2                tidyselect_1.2.1            miniUI_0.1.1.1              pbapply_1.7-2               gridExtra_2.3              
  # [101] IRanges_2.36.0              SummarizedExperiment_1.32.0 scattermore_1.2             stats4_4.3.2                Biobase_2.62.0             
  # [106] statmod_1.5.0               matrixStats_1.2.0           stringi_1.8.3               lazyeval_0.2.2              codetools_0.2-19           
  # [111] tibble_3.2.1                cli_3.6.2                   uwot_0.1.16                 xtable_1.8-4                reticulate_1.44.1          
  # [116] systemfonts_1.0.6           munsell_0.5.0               GenomeInfoDb_1.38.8         Rcpp_1.0.12                 globals_0.16.2             
  # [121] spatstat.random_3.2-2       png_0.1-8                   ggrastr_1.0.2               parallel_4.3.2              ellipsis_0.3.2             
  # [126] presto_1.0.0                dotCall64_1.1-1             sparseMatrixStats_1.14.0    bitops_1.0-7                listenv_0.9.1              
  # [131] viridisLite_0.4.2           scales_1.3.0                ggridges_0.5.6              crayon_1.5.2                leiden_0.4.3.1             
  # [136] packrat_0.9.2               purrr_1.0.2                 rlang_1.1.3                 cowplot_1.1.3
