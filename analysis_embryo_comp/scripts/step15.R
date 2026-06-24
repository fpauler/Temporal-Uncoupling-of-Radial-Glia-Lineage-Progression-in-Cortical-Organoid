# =======================================
# step 15
# create an overview UMAP to explain which 
# cells were extracted for downstream analyses
# =======================================

library ( Seurat )
library ( ggplot2 )

set.seed(2401)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

EB_neurons <-readRDS( file = paste(nobackup_base, "/RDS_files/EB_neurons.RDS", sep="") )
Seurat_neuron <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.RDS", sep="") )

umap_df <- as.data.frame( Embeddings( Seurat_neuron, reduction = "umap") )

cellIDs <- rownames( Embeddings( EB_neurons, reduction = "umap") )

umap_df$group <- "other"
umap_df$group[ which(rownames(umap_df) %in% cellIDs)] <- "EB Neurons"

ggplot(  ) + 
  geom_point( data = umap_df[ which( umap_df$group == "other" ), ], aes(x=umap_1, y=umap_2, color=group ), size=2 ) + 
  geom_point( data = umap_df[ which( umap_df$group == "EB Neurons" ), ], aes(x=umap_1, y=umap_2, color=group ), size=0.5 ) + 
  scale_color_manual ( values = c("EB Neurons" = "black", "other" = "grey50") ) +
  theme_classic()

ggsave( filename = paste(base, "/plots/Fig_4/Fig4i_Neuron_umap.png", sep=""), width = 15, height = 12 )
