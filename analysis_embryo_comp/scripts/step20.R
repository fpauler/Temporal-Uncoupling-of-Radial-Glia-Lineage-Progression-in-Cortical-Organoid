# ===================================
# step 20
# plot Emx1 expression in pseudobulk
# ===================================

library (Seurat)
library (ggplot2)
library (ggrepel)

set.seed(2401)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

seurat_sketch <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.downsample.RDS", sep=""))
Idents( seurat_sketch ) <- "toremove"
seurat_sketch <- subset( seurat_sketch, idents = "Yes", invert = T)

Idents( seurat_sketch ) <- "group"
embryo <- subset( seurat_sketch, idents = "embryo")
EB <- subset( seurat_sketch, idents = "EB")

# embryo@meta.data$cellType[ which( embryo@meta.data$cellType == "RGP" & embryo@meta.data$age == "E10") ] <- "RGP_1"
# embryo@meta.data$cellType[ which( embryo@meta.data$cellType == "RGP" & embryo@meta.data$age == "E13") ] <- "RGP_2"

embryo@meta.data$cellType[ which( embryo@meta.data$cellType %in% c("RGP_1" ,"RGP_2") ) ] <- "RGP"

# EB@meta.data$cellType[ which( EB@meta.data$cellType == "RGP" & EB@meta.data$age == "D8") ] <- "RGP_1"
# EB@meta.data$cellType[ which( EB@meta.data$cellType == "RGP" & EB@meta.data$age == "D13") ] <- "RGP_2"

EB@meta.data$cellType[ which( EB@meta.data$cellType %in% c("RGP_1" ,"RGP_2") ) ] <- "RGP"

embryo_aggrExpr <- AggregateExpression(
  embryo,
  group.by = "cellType", return.seurat = T)

EB_aggrExpr <- AggregateExpression(
  EB,
  group.by = "cellType", return.seurat = T)

eb_expr_df <- as.data.frame( EB_aggrExpr@assays$RNA$data["Emx1",] )
eb_expr_df$cellType <- rownames( eb_expr_df )
colnames( eb_expr_df ) <- c("in_vitro", "cellType")

embryo_expr_df <- as.data.frame( embryo_aggrExpr@assays$RNA$data["Emx1",] )
embryo_expr_df$cellType <- rownames( embryo_expr_df )
colnames( embryo_expr_df ) <- c("in_vivo", "cellType")

#sanity check
identical( eb_expr_df$cellType, embryo_expr_df$cellType)

# combine data for plotting
df2plot <- merge(eb_expr_df, embryo_expr_df, by="cellType")

ggplot( df2plot, aes(x=in_vitro, y=in_vivo, label=cellType)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point( size = 2 ) + geom_label_repel() + theme_classic() +
  ylim(0,2) + xlim(0,2) + ggtitle("Emx1 pseudobulk expression")
ggsave( paste(base, "plots/Sup_Fig_7/S7c_Emx1_pseudobulk_expression.pdf", sep=""), width = 5, height = 5)

