library (dplyr)
library (ggplot2)
library (loomR)
library (data.table)
library (Seurat)

set.seed( 2401 )
base <- "./"
nobackup_base <- "./"

# =======================
# define functions here
# =======================
get_pseudobulk <- function(sobj, by1 = "celltype", by2 = "stage"){
  expr <- GetAssayData(sobj, layer = "data")   # normalized data
  meta <- sobj@meta.data[, c(by1, by2)]
  
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
                                 label = lfile$col.attrs$Label[ regionIDX ],
                                 cellID = lfile$col.attrs$CellID[ regionIDX ])
# record the loom index for downstream analysis
all_abundances_df$loom_idx <- regionIDX

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

rand_idx <- sample( 1:nrow(all_abundances_df), 10000)
cell_idx <- all_abundances_df$loom_idx[ rand_idx ]

data.subset <- lfile[["matrix"]][ cell_idx, ]

#give cellID as row
rownames(data.subset) <- all_abundances_df$cellID[rand_idx]
#give Gene name as col - remove duplicates genes
colnames(data.subset) <- lfile$row.attrs$Gene[]
data.subset <- data.subset[ , which(!(duplicated(colnames(data.subset)) | duplicated(colnames(data.subset), fromLast = T))) ]

rownames(all_abundances_df) <- all_abundances_df$cellID

# create a Seurat object with this data
seurat.lm <- CreateSeuratObject(counts = t(data.subset), 
                                project = "LaManno_ref",
                                meta.data = all_abundances_df,
                                min.cells = 3, min.features = 200)

seurat.lm <- NormalizeData( seurat.lm )
seurat.lm <- FindVariableFeatures(seurat.lm, selection.method = "vst", nfeatures = 2000)

seurat.lm <- ScaleData(seurat.lm)
seurat.lm <- RunPCA(seurat.lm, features = VariableFeatures(object = seurat.lm))

seurat.lm <- RunUMAP(seurat.lm, dims = 1:25)

DimPlot( seurat.lm, group.by = "class") -> umap_class
FeaturePlot( seurat.lm, 
             features = c("Slc1a3", "Hes5", "Neurod1", "Neurog1", "Fabp7", "Rbfox3"), 
             order = T, min.cutoff = "q15") -> umap_marker_expression

ggsave( plot = umap_class, filename = "./umap.class.png", width=7, height=5)
ggsave( plot = umap_marker_expression, filename = "./featPlot.png", width=10, height=15)

peudobulk <- get_pseudobulk( seurat.lm, by1 = "class", by2 = "age")

saveRDS( object = peudobulk, file = "./LaManno_pseudobulk.rds")

sessionInfo()
