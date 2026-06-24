# =====================================================
# step 18
# this analysis uses embryo and EB cells together to map dev trajectories
# UMAPs and clusters are determined in the context of all cells
# this analysis produces images for a video that made a reviewer's figure
# =====================================================

library (Seurat)
library (slingshot)
library (Matrix)
library (ggplot2)
library(rgl)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# prepare figures with varying parameters for slingshot
seurat_merged <- readRDS (file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.RDS", sep=""))
DefaultAssay( seurat_merged) <- "RNA"

# plot 3d UMAP for the main figure
n_neighbors <- 10
min_dist <- 0.2 
res <- 0.8

     # prepare main UMAP as 3D
     seurat_merged <- RunUMAP(seurat_merged, reduction = "css", dims = 1:ncol(Embeddings(seurat_merged, "css")), 
                           n.neighbors = 10, metric = "cosine", n.components = 3,
                           min.dist = 0.2, return.model = F, seed.use = 2401, reduction.name = "umap_3d")
  
     seurat_merged <- FindClusters(seurat_merged, resolution = 0.8 )
     
     # identify surat clusters mainly belonging to earliest age
     # not integrated well and confuses the lineage analysis
     tot_freq_age <- table( seurat_merged@meta.data$age )
     
     age_assoc <- sapply( unique( as.character(seurat_merged@meta.data$seurat_clusters) ), function (cluster){
       freq <- table( seurat_merged@meta.data$age[ which(seurat_merged@meta.data$seurat_clusters == cluster)] )
       rel_freq <- freq[ c("E10", "E13", "E16", "P0") ] / tot_freq_age[ c("E10", "E13", "E16", "P0") ]
       if (any(!is.na(rel_freq))) {
         return ( names( which( rel_freq == max( rel_freq, na.rm = T ) ) ) )
       } else {
         return (NA)
       }
     })
     
     # restrict the analysis to only the neuronal cells in the embryo
     # otherwise the trajectories become way too complex
     meta_filter <- seurat_merged@meta.data
     
     # create colors for cell types as in fig 1
     col_list <- gg_color_hue ( length(unique(meta_filter$cellType)) )
     names(col_list) <- as.character( sort(unique(meta_filter$cellType)) )
     
     meta_filter$seurat_clusters <- as.character( meta_filter$seurat_clusters )
     meta_filter <- meta_filter[ which(meta_filter$seurat_clusters %in% names( age_assoc [ which(!age_assoc == "E10") ] )), ]
     #meta_filter <- meta_filter[ which(meta_filter$group == "EB"), ]
     
     lineages <- getLineages(
       data           = seurat_merged@reductions$umap_3d@cell.embeddings[ rownames(meta_filter), ],
       clusterLabels  = meta_filter$seurat_clusters,
       dist.method    = "simple"
     ) 
     
     # show the 3d trajectory
     plot3d( seurat_merged@reductions$umap_3d@cell.embeddings[ rownames(meta_filter), ], 
            col = col_list[ meta_filter$cellType ], size = 2, alpha = 0.2, decorate = F)
     plot3d.SlingshotDataSet( SlingshotDataSet(lineages), type='lineages', add=T, col = "black", linInd = 3 )
     plot3d.SlingshotDataSet( SlingshotDataSet(lineages), type='lineages', add=T, col = "grey", linInd = c(1:2, 4:8))
     
     # save viewpoint set manually
     if (!file.exists(paste(base, "/plots/Neuron_context_traj/view_matrix.RDS", sep=""))) {
       um <- par3d()$userMatrix
       saveRDS( object = um, file = paste(base, "/plots/Neuron_context_traj/view_matrix.RDS", sep=""))
     } else {
       um <- readRDS( file = paste(base, "/plots/Neuron_context_traj/view_matrix.RDS", sep="") )
     }
     
     # set viewpoint
     view3d(userMatrix = um)
     
     # rotate / make movie
     play3d(spin3d(axis = c(0, 0, 1), rpm = 2), duration = 100)
     movie3d(spin3d(axis = c(0, 0, 1), rpm = 2), duration = 100, dev = cur3d(), fps = 10, movie = "movie", 
             dir = paste(base, "plots/Neuron_context_traj/", sep=""), 
             convert = NULL, clean = F, verbose = TRUE, top = TRUE, 
             type = "gif", startTime = 0)

     nrow(meta_filter)
     # 116115
     
     # script to convert the images to animated GIF is in also in the folder
     
