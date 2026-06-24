# ========================================================
# step 16
# prepare waterfall diagram of top reference markers
# ========================================================

library (openxlsx)
library (dplyr)
library (ggplot2)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

layer_markers <- read.xlsx( xlsxFile = paste(base, "/plots/Sup_Fig_22/S22a_EB_ref_layer_marker_overlap.xlsx", sep=""), sheet = 2 )
layer_markers <- layer_markers[ which(layer_markers$p_val_adj == 0),]
layer_markers <- layer_markers[ order(layer_markers$avg_log2FC), ]
layer_markers$label <- layer_markers$gene
layer_markers$label[ which(!layer_markers$gene %in% c("Foxp2", "Sox5", "Bcl11b", "Pou3f3", "Satb2", "Cux2"))] <- ""
layer_markers$gene <- factor(layer_markers$gene, levels = layer_markers$gene)
layer_markers$color <- ifelse( layer_markers$avg_log2FC < 0, "cornflowerblue", "goldenrod")

ggplot( layer_markers, aes(x=gene, y=avg_log2FC, fill=color)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_identity() +
  scale_x_discrete( labels = layer_markers$label) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle( paste("top", nrow( layer_markers ), "UL/DL DEG") )

ggsave( paste(base, "/plots/Sup_Fig_22/S22a_reference_layer_marker_waterfall.pdf", sep=""), width=3)
