# ====================================================================================================
# step 10
# script to analyse age matching of EBs based on reference data
# done by plotting abundance of key cell types that have defined abundance changes during development
# ====================================================================================================

library (dplyr)
library (ggplot2)
library (loomR)
library (data.table)

set.seed( 2401 )

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

# load loom libraries
lfile <- connect(filename = paste( nobackup_base, "/dev_all.loom", sep=""), mode = "r+", skip.validate = T)

all_abundances_df <- data.frame()
# get cells from Forebrain
regionIDX <- which( lfile$col.attrs$Region[] == "Forebrain" )
regionIDX <- regionIDX[which(!(duplicated(lfile$col.attrs$CellID[regionIDX]) | duplicated(lfile$col.attrs$CellID[regionIDX], fromLast = T)))]

# identify meta data from these cells
all_abundances_df <- data.frame( age = lfile$col.attrs$Age[ regionIDX ], 
                                 class = lfile$col.attrs$Class[ regionIDX ], 
                                 clusterName = lfile$col.attrs$ClusterName[ regionIDX ],
                                 subclass = lfile$col.attrs$Subclass[ regionIDX ],
                                 location = lfile$col.attrs$Location_E9_E11[ regionIDX ],
                                 label = lfile$col.attrs$Label[ regionIDX ])

# remove low abundant cell types
all_abundances_df <- all_abundances_df[ which( !all_abundances_df$subclass %in% names( which( table( all_abundances_df$subclass ) < 500 ) ) ), ]

# include only key cell types that change in abundance with developmental time
cellTypes_to_include <- c("Cortical hem", "Cortical or hippocampal glutamatergic", "Dorsal forebrain", 
                          "Forebrain", "Forebrain astrocyte", "Forebrain glutamatergic", 
                          "Neuronal intermediate progenitor",
                          "Committed oligodendrocyte precursor", "Oligodendrocyte precursor cell", "Oligodendrocyte")

all_abundances_df <- all_abundances_df[ which(all_abundances_df$subclass %in% cellTypes_to_include), ]

# remove cell types with ambiguous origin
all_abundances_df <- all_abundances_df[ which(!all_abundances_df$class %in% c("Ependymal", "Glioblast", "Cajal-Retzius")), ]

# truncate developmental age
all_abundances_df$simple_age <- sapply( all_abundances_df$age, function (x) strsplit( x = x, split = "\\.")[[1]][1] )
all_abundances_df$simple_age[ which(all_abundances_df$simple_age == "e9") ] <- "e09"
# combine class labels
all_abundances_df$class[ which(grepl(pattern = "astrocyte", x = all_abundances_df$subclass )) ] <- "Glia"
all_abundances_df$class[ which(grepl(pattern = "Oligodendrocyte", x = all_abundances_df$subclass )) ] <- "Glia"

unique( all_abundances_df$class )
# "Radial glia" "Glia"        "Neuroblast"  "Neuron" 

# calculate abundance and bring into shape useful for plotting
all_abundances_df %>% 
  group_by(simple_age, class) %>% 
  summarise(n = n()) -> cellType_summary

all_rel_abundance_lM <- data.frame()
for (age in unique( cellType_summary$simple_age)) {
  tmp_df <- cellType_summary[ which( cellType_summary$simple_age == age ), ]
  tmp_df$mean <- tmp_df$n / sum( tmp_df$n )
  tmp_df$sd <- NA
  tmp_df$orig <- "LaManno"
  all_rel_abundance_lM <- rbind( all_rel_abundance_lM, tmp_df )
}

table( all_abundances_df$simple_age )
#  e09   e10   e11   e12   e13   e14   e15   e16   e17   e18 
# 1133  3434  3032 11776  7176  6490 10683  9145  8647  6292 

# diBella reference
meta_data <- fread(paste(nobackup_base, "DiBella/metaData_scDevSC.txt", sep=""))
meta_data <- as.data.frame(meta_data)
meta_data <- meta_data[2:nrow(meta_data),]

# remove spaces from ttype IDs - makes problems with saving intermediate files
meta_data$New_cellType <- gsub( pattern = " ", replacement = "_", x = meta_data$New_cellType)

cell_types_to_include <- c("Apical_progenitors", "Intermediate_progenitors", "Immature_neurons",
                           "Migrating_neurons", "SCPN", "CThPN", "DL_CPN", "UL_CPN", "Layer_4", "NP", "Layer_6b",
                           "Astrocytes", "Oligodendrocytes")

meta_data <- meta_data[ which(meta_data$New_cellType %in% cell_types_to_include),]

meta_data$class <- meta_data$New_cellType

neuron_class_vec <- c("Immature_neurons", "Migrating_neurons", "SCPN", "CThPN", "DL_CPN", "UL_CPN", "Layer_4", "NP", "Layer_6b")
glia_class_vec <- c("Astrocytes", "Oligodendrocytes")

meta_data$class[ which( meta_data$class %in% neuron_class_vec)] <- "Neuron"
meta_data$class[ which( meta_data$class %in% glia_class_vec)] <- "Glia"

meta_data$orig_ident[ which( meta_data$orig_ident %in% c("E18_S1", "E18_S3"))] <- "E18"
meta_data$orig_ident[ which( meta_data$orig_ident %in% c("P1_S1", "P1"))] <- "P1"

table( meta_data$orig_ident )
#  E10   E11   E12   E13   E14   E15   E16   E17   E18    P1    P4 
# 2395  3710  7753  8403  4002  9961  4864  7542 14297  9422  4806 

meta_data %>% 
  group_by(orig_ident, class) %>% 
  summarise(n = n()) -> cellType_summary_dB

all_rel_abundance_dB <- data.frame()
for (age in unique( cellType_summary_dB$orig_ident)) {
  tmp_df <- cellType_summary_dB[ which( cellType_summary_dB$orig_ident == age ), ]
  colnames(tmp_df) <- c("simple_age", "class", "n")
  tmp_df$mean <- tmp_df$n / sum( tmp_df$n )
  tmp_df$sd <- NA
  tmp_df$orig <- "diBella"
  all_rel_abundance_dB <- rbind( all_rel_abundance_dB, tmp_df )
}

# combine the plots
age2idx <- 1:12
names( age2idx) <- c("E09", "E10", "E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "P1", "P4")

# make cell descriptions similar
unique( all_rel_abundance_dB$class )
all_rel_abundance_dB$class[ which( all_rel_abundance_dB$class == "Apical_progenitors") ] <- "Radial glia"
all_rel_abundance_lM$class[ which( all_rel_abundance_lM$class == "Neuroblast") ] <- "Intermediate_progenitors"

df2plot <- rbind(all_rel_abundance_dB, all_rel_abundance_lM)
df2plot$idx <- age2idx[ toupper(df2plot$simple_age) ]
ggplot( df2plot, aes(x=idx, y=mean, group=orig, color=orig)) + geom_point() + facet_grid(~class)
ggsave(filename = paste(base, "plots/QC/reference_cellType_abundance.pdf", sep=""), width=8)

# our data
meta_data <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.meta_data.RDS", sep=""))
# remove our embryo data from this analysis
meta_data <- meta_data[ which(meta_data$group %in% c("EB")), ]

meta_data <- meta_data[ which(meta_data$cellType %in% c("RGP", "IP", "iN", "astro", "oligo")), ]
meta_data$cellType[which(meta_data$cellType == "RGP")] <- "Radial glia"
meta_data$cellType[which(meta_data$cellType == "IP")] <- "Intermediate_progenitors"
meta_data$cellType[which(meta_data$cellType == "iN")] <- "Neuron"
meta_data$cellType[which(meta_data$cellType %in% c("astro", "oligo"))] <- "Glia"

cellTypeAbundance <- data.frame()
for (orig_ident in unique(meta_data$orig.ident)) {
  for (cell_type in unique(meta_data$cellType) ) {
    abundance <- length ( which( meta_data$cellType == cell_type & meta_data$orig.ident == orig_ident ) )
    tmp_df <- data.frame( abundance = abundance, cell_type = cell_type, group = orig_ident)
    cellTypeAbundance <- rbind( cellTypeAbundance, tmp_df )
  }
}

abundance_vec <- table( meta_data$orig.ident )
cellTypeAbundance$tot_cells <- abundance_vec[ cellTypeAbundance$group ]
cellTypeAbundance$rel <- cellTypeAbundance$abundance / cellTypeAbundance$tot_cells
cellTypeAbundance$age <- sapply( cellTypeAbundance$group, function (x) strsplit(x, "-")[[1]][1] )

cellTypeAbundance %>%
  group_by( age, cell_type ) %>%
  summarize(mean = mean(rel), sd = sd(rel),
            n = n()) -> cellTypeAbundanceSummary

colnames( cellTypeAbundanceSummary ) <- c( "simple_age", "class", "mean", "sd", "n")
cellTypeAbundanceSummary$orig <- "EB"
cellTypeAbundance$orig <- "EB"
colnames(cellTypeAbundance) <- c("abundance", "class",  "group", "tot_cells", "rel", "age", "orig")

age2idx_this <- c(2,5,8,11)
names(age2idx_this) <- c("D8", "D13", "D20", "D25")
cellTypeAbundanceSummary$idx <- age2idx_this[ cellTypeAbundanceSummary$simple_age ]
cellTypeAbundance$idx <- age2idx_this[ cellTypeAbundance$age ] 

df2plot <- rbind( df2plot, cellTypeAbundanceSummary[, colnames(df2plot)] )

# do the plotting 
# indicate only relevant developmental ages
# show reference only as smooth
ggplot() + 
  geom_point( data = cellTypeAbundance, aes(x=idx, y=rel, color=orig)) + 
  geom_smooth(data = df2plot[which(df2plot$orig %in% c("diBella", "LaManno")),], aes(x=idx, y=mean), linewidth = 1, fill = "grey70", colour = "grey50", level=0.9) + 
  geom_line( data = cellTypeAbundanceSummary, aes(x=idx, y=mean, color=orig, group=orig)) +
  geom_errorbar( data = cellTypeAbundanceSummary, aes(ymin=mean-sd, ymax=mean+sd, x=idx)) +
  scale_x_continuous( breaks = c(2,5,8,11,12), labels = c("E10", "E13", "E16", "P0", "P4")) +
  facet_grid(~class) + theme_classic()

ggsave( paste(base, "plots/Sup_Fig_4/S4f_cellType_abundance.pdf", sep=""), width=8, height=3 )

ggplot() + 
  geom_point( data = df2plot[which(df2plot$orig %in% c("diBella", "LaManno")),], aes(x=idx, y=mean, color=orig)) + 
  geom_smooth(data = df2plot, aes(x=idx, y=mean), linewidth = 1, fill = "grey70", colour = "grey50", level=0.9) + 
  facet_grid(~class) + theme_classic()
#`geom_smooth()` using method = 'loess' and formula = 'y ~ x'
ggsave(filename = paste(base, "plots/QC/reference_cellType_abundance_with _smooth.pdf", sep=""), width=8)

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
#   [1] data.table_1.15.0 loomR_0.2.1.9000  hdf5r_1.3.10      R6_2.5.1          ggplot2_3.5.0     dplyr_1.1.4      
# 
# loaded via a namespace (and not attached):
#   [1] Matrix_1.6-5      bit_4.0.5         gtable_0.3.4      compiler_4.3.2    tidyselect_1.2.1  stringr_1.5.1    
# [7] parallel_4.3.2    textshaping_0.3.7 splines_4.3.2     systemfonts_1.0.6 scales_1.3.0      lattice_0.21-9   
# [13] labeling_0.4.3    generics_0.1.3    tibble_3.2.1      munsell_0.5.0     pillar_1.9.0      rlang_1.1.3      
# [19] utf8_1.2.4        stringi_1.8.3     bit64_4.0.5       cli_3.6.2         withr_3.0.0       magrittr_2.0.3   
# [25] mgcv_1.9-0        grid_4.3.2        rstudioapi_0.16.0 pbapply_1.7-2     packrat_0.9.2     lifecycle_1.0.4  
# [31] nlme_3.1-163      vctrs_0.6.5       glue_1.7.0        farver_2.1.1      ragg_1.3.0        fansi_1.0.6      
# [37] colorspace_2.1-0  tools_4.3.2       pkgconfig_2.0.3  
