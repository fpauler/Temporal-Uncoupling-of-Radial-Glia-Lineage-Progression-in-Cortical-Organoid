# =========================================
# step 08
# plot basic QC features of scRNA-Seq data
# =========================================

library (dplyr)
library (ggplot2)
library (Seurat)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

  seurat_merged <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.RDS", sep=""))
  
  # plot cell number per sample
  cell_number_df <- as.data.frame( table( seurat_merged@meta.data$orig.ident) )
  cell_number_df$age <- sapply( cell_number_df$Var1, function (x) strsplit( x = as.character(x), split = "-")[[1]][1] )
  cell_number_df$age <- gsub( pattern = "D8", replacement = "D08", x = cell_number_df$age )
  cell_number_df$line <- sapply( cell_number_df$Var1, function (x) strsplit( x = as.character(x), split = "-")[[1]][2] )
  cell_number_df$line <- gsub( pattern = "line", replacement = "", x = cell_number_df$line)
  cell_number_df$age <- paste( cell_number_df$age, cell_number_df$line, sep="_")
  
  cell_number_plot <- ggplot( cell_number_df, aes(x = age, y=Freq)) + geom_bar(stat="identity", position = "dodge") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,10000) +
    ggtitle("# cells per sample") + ylab("")
  
  seurat_merged@meta.data$line <- gsub( pattern = "line", replacement = "", x = seurat_merged@meta.data$orig.ident)
  seurat_merged@meta.data$line <- gsub( pattern = "D8", replacement = "D08", x = seurat_merged@meta.data$line)
  
  seurat_merged@meta.data$line<- factor( seurat_merged@meta.data$line, levels = sort( unique( seurat_merged@meta.data$line )))
  
  ggplot( seurat_merged@meta.data, aes(x=line, y=nCount_RNA)) + 
    geom_boxplot(outlier.size = 1, outlier.shape = 4) + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("# transcripts detected") + ylab("") + xlab("") -> plot_nCount
  
  ggplot( seurat_merged@meta.data, aes(x=line, y=nFeature_RNA)) + 
    geom_boxplot(outlier.size = 1, outlier.shape = 4) + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("# genes detected") + ylab("") + xlab("") -> plot_nFeatures
  
  ggplot( seurat_merged@meta.data, aes(x=line, y=percent.mt)) + 
    geom_boxplot(outlier.size = 1, outlier.shape = 4) + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("% mitochondria reads") + ylab("") + xlab("") -> plot_percMt
  

  plot_nCount / plot_nFeatures / plot_percMt / cell_number_plot
  
  ggsave( filename = paste(base, "plots/Sup_Fig_4/S4abcd_feature_plots_rev.pdf", sep=""), width=4, height=8)
  
  # write out raw data
  write.csv( x = seurat_merged@meta.data, file = paste(base, "/plots/Sup_Fig_4/S4abcd_plot_raw_data.csv", sep=""), row.names = T) 
  
