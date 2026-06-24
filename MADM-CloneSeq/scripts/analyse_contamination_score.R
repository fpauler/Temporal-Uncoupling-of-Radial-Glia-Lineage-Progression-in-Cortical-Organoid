# ==============================================================
#
# analyse cells removed due to non excitatory neuron signature
# this analysis gives a hint as to how many cells of a 
# MADM clone are actually neurons
# only used for response to reviewers
#
# ==============================================================

library (openxlsx)
library (ggplot2)
library (dplyr)
library (ggbeeswarm)
library (ggupset)
library (patchwork)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "/MADM-CloneSeq/", sep="" )

meta_data <- read.xlsx(xlsxFile = paste( base, "/Supplement/PatchSeq_metadata.xlsx", sep="") )

# remove clones with unclear clone association
#meta_data <- meta_data[ which( !meta_data$clone == "na" ), ]

nrow( meta_data )
# 279

# ----------------------------------
# overview figure of broad category
# ----------------------------------
table( meta_data$category )
# cont Neuron 
#   64    215 

cat2plot <- as.data.frame( table( meta_data$category ) / nrow( meta_data ) )

ggplot( cat2plot, aes(x=Var1, y=Freq)) + geom_bar( stat="identity", position="dodge" ) +
  ylim(0,0.8) + theme_classic()
ggsave( filename = paste(base, "/plots/PatchSeq_category_plot.pdf", sep=""))

# -------------------------------------
# focus on non excitatory neuron cells
# -------------------------------------

# extract cells labeled as contamination before
cont_meta <- meta_data[ which( meta_data$category == "cont" ), ]

nrow(cont_meta)
# 64

# identify the most likely cell type of this contamination cell
cont_max_id <- sapply( 1:nrow(cont_meta), function (x) {
  
  colOi <- c("NMS_aNSC","NMS_Astro","NMS_IP","NMS_OBNB","NMS_oligo")
  tmp <- cont_meta[x, colOi]
  idx <- which( tmp == max(tmp) )
  return( colOi[idx])
})
cont_meta$contID <- cont_max_id

cont_max <- sapply( 1:nrow(cont_meta), function (x) {
  
  colOi <- c("NMS_aNSC","NMS_Astro","NMS_IP","NMS_OBNB","NMS_oligo")
  tmp <- cont_meta[x, colOi]

  return( max(tmp) )
})
cont_meta$NMS_cont <- cont_max

neuron_meta <- meta_data[ which( meta_data$category == "Neuron" ), ]

# ----------------------------------------------------------------------------------
# plot the score of non-excitatory Neurons / Neurons
# some non-excitatoy Neurons show scores similar to neurons
# these are likely real cells and can help to indicate the % Neurons in each clone
# ----------------------------------------------------------------------------------

df2plot <- data.frame( score = c(neuron_meta$NMS_Neurons, cont_meta$NMS_cont), group = c(neuron_meta$category, cont_meta$category))

ggplot( df2plot, aes( score, after_stat(density), group = group, fill = group)) + 
  geom_histogram(linewidth = 0.5, position = "dodge") + theme_classic()
ggsave( filename = paste(base, "/plots/PatchSeq_NMSscore_plot.pdf", sep=""))

# true contamination
length( which(cont_meta$NMS_cont <= min( neuron_meta$NMS_Neurons ) ) ) / nrow(cont_meta) # 31 cells
# 0.4761905

# probably true cell type different than neuron
length( which(cont_meta$NMS_cont > min( neuron_meta$NMS_Neurons ) ) ) / nrow(cont_meta) # 33 cells
# 0.5238095

# ----------------------------------------------------------------------
#
# plot the cell type of the cell that show a high contamination score
# most likely these cells are 'real', non neuronal cells within a clone
#
# ----------------------------------------------------------------------

freq_df <- as.data.frame(table( cont_max_id[ which(cont_meta$NMS_cont > min( neuron_meta$NMS_Neurons ) ) ] ) )
freq_df <- freq_df[order(freq_df$Freq, decreasing = T),]
freq_df$rel <- freq_df$Freq / sum( freq_df$Freq )
freq_df$Var1 <- factor(freq_df$Var1, levels=freq_df$Var1)

# actual data:
# Var1 Freq        rel
# 1  NMS_aNSC   12 0.36363636
# 3    NMS_IP    9 0.27272727
# 4  NMS_OBNB    6 0.18181818
# 5 NMS_oligo    4 0.12121212
# 2 NMS_Astro    2 0.06060606

ggplot( freq_df, aes(x=Var1, y=rel)) + geom_bar( stat="identity", position="dodge") + ylim(0, 0.4) + 
  theme_classic()
ggsave( filename = paste(base, "/plots/PatchSeq_hq_cont_cells_identity.pdf", sep=""))

# -----------------------------------------
# how many high quality cells are neurons
# -----------------------------------------

high_qual_cont_scores <- cont_meta[ which(cont_meta$NMS_cont > min( neuron_meta$NMS_Neurons )), ]
nrow( neuron_meta ) / (nrow( neuron_meta ) + nrow( high_qual_cont_scores ) )
# [1] 0.8669355

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
#   [1] patchwork_1.2.0  ggbeeswarm_0.7.2 ggupset_0.4.1    dplyr_1.1.4      ggplot2_3.5.0    openxlsx_4.2.5.2
# 
# loaded via a namespace (and not attached):
#   [1] vctrs_0.6.5       zip_2.3.1         cli_3.6.2         rlang_1.1.3       packrat_0.9.2     stringi_1.8.3     generics_0.1.3   
# [8] textshaping_0.3.7 glue_1.7.0        labeling_0.4.3    colorspace_2.1-0  ragg_1.3.0        scales_1.3.0      fansi_1.0.6      
# [15] grid_4.3.2        munsell_0.5.0     tibble_3.2.1      lifecycle_1.0.4   vipor_0.4.7       compiler_4.3.2    Rcpp_1.0.12      
# [22] pkgconfig_2.0.3   rstudioapi_0.16.0 beeswarm_0.4.0    systemfonts_1.0.6 farver_2.1.1      R6_2.5.1          tidyselect_1.2.1 
# [29] utf8_1.2.4        pillar_1.9.0      magrittr_2.0.3    tools_4.3.2       withr_3.0.0       gtable_0.3.4     
