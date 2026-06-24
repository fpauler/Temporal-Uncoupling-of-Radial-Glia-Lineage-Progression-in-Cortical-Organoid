# ====================================================================================
# step 09
# prepare a plot to compare dynamics and relative abundance of cell type production 
# between organoid and embryo data (including reference)
# this demands to annotate LaManno et al./diBella et al. reference based on our data
# use E11.5, E13.5, E16, P1 from diBella
# use e11.0, e13.5, e16.5, e18.0 from LaManno
# ====================================================================================

# finally a plot is prepared to display the data

library (Seurat)
library (dplyr)
library (ggplot2)
library (clusterProfiler)
library (org.Mm.eg.db)

set.seed( 2401 )

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

# read our processed data with annotation
seurat_sketch <- readRDS ( file = paste(nobackup_base, "RDS_files/all_seurat_merged.sketched.downsample.RDS", sep=""))

# cell cycle gene information
# read cell cycle genes to regress out cell cycle effects
cell_cycle_phase_gene <- read.table(file = paste(base, "other_data/Mus_musculus.csv", sep=""), sep = ",", header = T)
gene_conversion <- bitr(geneID = cell_cycle_phase_gene$geneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
cell_cycle_phase_gene <- merge(cell_cycle_phase_gene, gene_conversion, by.x="geneID", by.y="ENSEMBL")

# first do LaManno data

if (!file.exists( paste(nobackup_base, "RDS_files/LaManno_thisCTAnno.RDS", sep="") )) {
  
  # here the samples for annotation are defined
  # note that tissue specificity comes from these sample choices
  sampleID_list <- list()
  sampleID_list[[ "e11.0" ]] <- c( "10X40_3", "10X40_4", "10X40_5", "10X40_6" )
  sampleID_list[[ "e13.5" ]] <- c( "10X14_4" )
  sampleID_list[[ "e16.5" ]] <- c( "10X17_4" )
  sampleID_list[[ "e18.0" ]] <- c( "10X32_2" )
  
  age_list <- c("E10", "E13", "E16", "P0")
  names( age_list ) <- c("e11.0","e13.5", "e16.5", "e18.0")

  # load loom libraries
  library (loomR)
  # process loom file
  lfile <- connect(filename = paste( nobackup_base, "/dev_all.loom", sep=""), mode = "r+", skip.validate = T)
  
  meta_df_list <- list()
  UMAP_plot_list <- list()
  # loop through developmental ages
  for ( age in names( sampleID_list ) ) {
    
    # extract reference cells
    cellIDX <- which( lfile$col.attrs$Age[] == age )
    cellIDX <- cellIDX[which(!(duplicated(lfile$col.attrs$CellID[cellIDX]) | duplicated(lfile$col.attrs$CellID[cellIDX], fromLast = T)))]
    
    tissue_df <- data.frame()
    for ( sampleID in unique( lfile$col.attrs$SampleID[cellIDX] ) ) {
      SampleIDidx <- which( lfile$col.attrs$SampleID[] == sampleID )
      tmpIdx <- intersect(cellIDX, SampleIDidx)
      tmp_df <- as.data.frame( table( lfile$col.attrs$Tissue[ tmpIdx ] ) )
      tmp_df$sampleIdx <- sampleID
      tissue_df <- rbind( tissue_df, tmp_df)
    }
    
    sampleIDX <- which( lfile$col.attrs$SampleID[] %in% sampleID_list[[ age ]] )
    
    cellTypes_to_include <- c("Cajal-Retzius", "Cortical hem", "Cortical or hippocampal glutamatergic", "Dorsal forebrain", 
                              "Forebrain", "Forebrain astrocyte", "Forebrain glutamatergic", 
                              "Neuronal intermediate progenitor", "Glioblast",
                              "Mixed region astrocytes", 
                              "Committed oligodendrocyte precursor", "Oligodendrocyte precursor cell", "Oligodendrocyte" )

    cellTypeIdx <- which( lfile$col.attrs$Subclass[] %in% cellTypes_to_include )
    
    cellIDX <- intersect( cellIDX, sampleIDX)
    cellIDX <- intersect( cellIDX, cellTypeIdx)
    
    data.subset <- lfile[["matrix"]][ cellIDX, ]
    
    #give cellID as row
    rownames(data.subset) <- lfile$col.attrs$CellID[cellIDX]
    #give Gene name as col - remove duplicates genes
    colnames(data.subset) <- lfile$row.attrs$Gene[]
    
    data.subset <- data.subset[ , which(!(duplicated(colnames(data.subset)) | duplicated(colnames(data.subset), fromLast = T))) ]
    
    meta_data <- list()
    meta_data[["cell_type"]] <- lfile$col.attrs$Class[cellIDX]
    names(meta_data[["cell_type"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["Age"]] <- lfile$col.attrs$Age[cellIDX]
    names(meta_data[["Age"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["Tissue"]] <- lfile$col.attrs$Tissue[cellIDX]
    names(meta_data[["Tissue"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["Location"]] <- lfile$col.attrs$Location_E9_E11[cellIDX]
    names(meta_data[["Location"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["clusterID"]] <- lfile$col.attrs$ClusterName[cellIDX]
    names(meta_data[["clusterID"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    meta_data[["subclass"]] <- lfile$col.attrs$Subclass[ cellIDX ]
    names(meta_data[["subclass"]]) <- lfile$col.attrs$CellID[cellIDX]
    
    Idents( seurat_sketch ) <- "orig.ident"
    tmp_embryo_obj <- subset(seurat_sketch, idents = age_list[ age ] )
    
    # get raw counts to prepare new objects with matching gene sets
    embryo_counts <- GetAssayData(object = tmp_embryo_obj, assay = "RNA", layer = "counts")
    
    # determin common genes - use of identical gene lists is necessary for later integration
    common_genes <- intersect(rownames(embryo_counts), colnames(data.subset))
    
    thisStudy_Seurat <- CreateSeuratObject(counts = embryo_counts[common_genes, ], project = paste(age, "_mTmG", sep=""), min.cells = 3, min.features = 200)
    LaManno_Seurat <- CreateSeuratObject(counts = as.matrix(t(data.subset[, common_genes ])), project = paste(age, "_LaManno", sep=""), min.cells = 3, min.features = 200)
    
    thisStudy_Seurat <- AddMetaData( thisStudy_Seurat, tmp_embryo_obj@meta.data)
    
    # Process LaManno reference data
    
    # CellCycleScoring 
    LaManno_Seurat <- NormalizeData( LaManno_Seurat )
    LaManno_Seurat <- CellCycleScoring (
      object = LaManno_Seurat,
      g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
      s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
    )
    
    # regress out cell cycle effects
    LaManno_Seurat$CC.Difference <- LaManno_Seurat$S.Score - LaManno_Seurat$G2M.Score
    LaManno_Seurat <- SCTransform(LaManno_Seurat, vars.to.regress = c("CC.Difference"), verbose = FALSE)
    LaManno_Seurat <- RunPCA(LaManno_Seurat, assay = "SCT")
    
    LaManno_Seurat <- RunUMAP(LaManno_Seurat, reduction = "pca", dims = 1:15, reduction.name = "umap")
    
    # process data from this study for integration
    
    # CellCycleScoring 
    thisStudy_Seurat <- NormalizeData(thisStudy_Seurat)
    thisStudy_Seurat <- CellCycleScoring (
      object = thisStudy_Seurat, 
      g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
      s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
    )
    
    # regress out cell cycle effects
    thisStudy_Seurat$CC.Difference <- thisStudy_Seurat$S.Score - thisStudy_Seurat$G2M.Score
    thisStudy_Seurat <- SCTransform(thisStudy_Seurat, vars.to.regress = c("CC.Difference"), verbose = FALSE)
    thisStudy_Seurat <- RunPCA( thisStudy_Seurat, assay = "SCT", npcs = 30)
    
    # do the integration / label transfer
    anchors <- FindTransferAnchors(reference = thisStudy_Seurat, query = LaManno_Seurat, 
                                   normalization.method = "SCT",
                                   reference.reduction = "pca", dims = 1:15)
    
    predictions_cellType <- TransferData(anchorset = anchors, refdata = thisStudy_Seurat$cellType, dims = 1:15)
    
    # add the labels to the original UMAP - sanity check not saved
    LaManno_Seurat <- AddMetaData(LaManno_Seurat, metadata = predictions_cellType)
    UMAP_plot_list[[age]] <- DimPlot( LaManno_Seurat, group.by = c("predicted.id") )
    
    meta_df_list[[ age ]] <- LaManno_Seurat@meta.data
  }
  
  saveRDS ( object = meta_df_list, file = paste(nobackup_base, "RDS_files/LaManno_thisCTAnno.RDS", sep=""))
  
}

# diBella
if (!file.exists( paste(nobackup_base, "RDS_files/diBella_thisCTAnno.RDS", sep="") )) {
  
  library( data.table)
  
  # metaData_scDevSC.txt downloaded from single cell portal
  # at the time of download (03/2022):
  # https://singlecell.broadinstitute.org/single_cell/study/SCP1290/molecular-logic-of-cellular-diversification-in-the-mammalian-cerebral-cortex
  meta_data <- fread(paste(nobackup_base, "DiBella/metaData_scDevSC.txt", sep=""))
  meta_data <- as.data.frame(meta_data)
  meta_data <- meta_data[2:nrow(meta_data),]
  
  # remove spaces from ttype IDs - makes problems with saving intermediate files
  meta_data$New_cellType <- gsub( pattern = " ", replacement = "_", x = meta_data$New_cellType)
  
  # files come from GSE153162_RAW.tar downloaded from GEO
  file_list <- c("GSM4635072_E11_5_filtered_gene_bc_matrices_h5.h5", "GSM4635073_E12_5_filtered_gene_bc_matrices_h5.h5",
                 "GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5", "GSM4635075_E14_5_filtered_gene_bc_matrices_h5.h5",
                 "GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5", "GSM4635077_E16_filtered_gene_bc_matrices_h5.h5",
                 "GSM5277844_E17_5_filtered_feature_bc_matrix.h5", "GSM4635078_E18_5_S1_filtered_gene_bc_matrices_h5.h5",
                 "GSM4635079_E18_S3_filtered_gene_bc_matrices_h5.h5", "GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5",
                 "GSM4635081_P1_S2_filtered_gene_bc_matrices_h5.h5", "GSM5277843_E10_v1_filtered_feature_bc_matrix.h5" )
  
  age_list <- c("E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18_S1", "E18_S3", "P1_S1", "P1", "E10")
  ct_prefix <- c("E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E18", "P1", "P1", "E10")
  
  cell_types_to_include <- c("Apical_progenitors", "Intermediate_progenitors", "Immature_neurons", "Cajal_Retzius_cells",
                             "Migrating_neurons", "SCPN", "CThPN", "DL_CPN", "UL_CPN", "Layer_4", "NP", "Layer_6b",
                             "Cycling_glial_cells", "Astrocytes", "Oligodendrocytes")
  
  meta_df_list <- list()
  UMAP_plot_list <- list()
  
  for (i in c(1, 3, 6, 10)) {
    
    fileName <- file_list[i]
    
    sc.data <- Read10X_h5(filename = paste(nobackup_base, "DiBella/", fileName, sep=""))
    
    sub_meta_data <- meta_data[ which( meta_data$biosample_id ==  age_list[i] & meta_data$New_cellType %in% cell_types_to_include ), ]
    
    cellIDs <- sub_meta_data$NAME
    
    cellIDs <- gsub( pattern = "-1", replacement = "", x = cellIDs )
    
    cellIDs <- sapply(cellIDs, function (x) {
      tmp <- strsplit(x = x, split = "_")[[1]]
      return ( tmp[length(tmp)] )
    })
    
    cellIDs <- paste( cellIDs, 1, sep="-" )
    
    rownames(sub_meta_data) <- cellIDs
    
    colIDX <- which( colnames(sc.data) %in% cellIDs)
    
    Idents( seurat_sketch ) <- "orig.ident"
    
    # matching developmental ages
    if ( age_list[i] == "E11") {
      tmp_embryo_obj <- subset(seurat_sketch, idents = "E10" )
    } else if ( age_list[i] == "P1_S1") {
      tmp_embryo_obj <- subset(seurat_sketch, idents = "P0" )
    } else {
      tmp_embryo_obj <- subset(seurat_sketch, idents = age_list[i] )
    }
    
    # get raw counts to prepare new objects with matching gene sets
    embryo_counts <- GetAssayData(object = tmp_embryo_obj, assay = "RNA", layer = "counts")
    
    # determine common genes - use of identical gene lists is necessary for later integration
    common_genes <- intersect(rownames(embryo_counts), rownames(sc.data))
    
    thisStudy_Seurat <- CreateSeuratObject(counts = embryo_counts[common_genes, ], project = paste(age_list[i], "_mTmG", sep=""), min.cells = 3, min.features = 200)
    diBella_Seurat <- CreateSeuratObject(counts = sc.data[common_genes,colIDX], project = paste(age_list[i], "_diBella", sep=""), min.cells = 3, min.features = 200)
    
    thisStudy_Seurat <- AddMetaData( thisStudy_Seurat, tmp_embryo_obj@meta.data)
    
    # Process diBella reference data 
    # CellCycleScoring 
    diBella_Seurat <- NormalizeData( diBella_Seurat )
    diBella_Seurat <- CellCycleScoring (
      object = diBella_Seurat,
      g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
      s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
    )
    
    # regress out cell cycle effects
    diBella_Seurat$CC.Difference <- diBella_Seurat$S.Score - diBella_Seurat$G2M.Score
    diBella_Seurat <- SCTransform(diBella_Seurat, vars.to.regress = c("CC.Difference"), verbose = FALSE)
    diBella_Seurat <- RunPCA(diBella_Seurat, assay = "SCT")
    
    diBella_Seurat <- RunUMAP(diBella_Seurat, reduction = "pca", dims = 1:15, reduction.name = "umap")
    
    # process data from this study for integration
    # CellCycleScoring 
    thisStudy_Seurat <- NormalizeData(thisStudy_Seurat)
    thisStudy_Seurat <- CellCycleScoring (
      object = thisStudy_Seurat, 
      g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
      s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
    )
    
    # regress out cell cycle effects
    thisStudy_Seurat$CC.Difference <- thisStudy_Seurat$S.Score - thisStudy_Seurat$G2M.Score
    thisStudy_Seurat <- SCTransform(thisStudy_Seurat, vars.to.regress = c("CC.Difference"), verbose = FALSE)
    thisStudy_Seurat <- RunPCA( thisStudy_Seurat, assay = "SCT", npcs = 30)
    
    # do the integration / label transfer
    anchors <- FindTransferAnchors(reference = thisStudy_Seurat, query = diBella_Seurat, 
                                   normalization.method = "SCT",
                                   reference.reduction = "pca", dims = 1:15)
    
    predictions_cellType <- TransferData(anchorset = anchors, refdata = thisStudy_Seurat$cellType, dims = 1:15)
    
    # add the labels to the original UMAP
    diBella_Seurat <- AddMetaData(diBella_Seurat, metadata = predictions_cellType)
    
    # add the labels to the original UMAP - sanity check not saved
    UMAP_plot_list[[ age_list[i] ]] <- DimPlot( diBella_Seurat, group.by = c("predicted.id") )
    
    meta_df_list[[ age_list[i] ]] <- diBella_Seurat@meta.data
    
  }
  
  saveRDS ( object = meta_df_list, file = paste(nobackup_base, "RDS_files/diBella_thisCTAnno.RDS", sep=""))
}

#####################
##
# finally plotting
##
#####################

# read the references
meta_df_list_dB <- readRDS ( file = paste(nobackup_base, "RDS_files/diBella_thisCTAnno.RDS", sep=""))
meta_df_list_LaManno <- readRDS ( file = paste(nobackup_base, "RDS_files/LaManno_thisCTAnno.RDS", sep=""))

# prepare count tables for each dataset
df2plot <- data.frame()
for (idx in names(meta_df_list_dB)) {
  tmp_df <- as.data.frame( table( meta_df_list_dB[[ idx ]]$predicted.id ) )
  tmp_df$age <- idx
  tmp_df$rel <- tmp_df$Freq / sum(tmp_df$Freq)
  tmp_df$study <- "diBella"
  df2plot <- rbind(df2plot, tmp_df)
}

for (idx in names(meta_df_list_LaManno)) {
  tmp_df <- as.data.frame( table( meta_df_list_LaManno[[ idx ]]$predicted.id ) )
  tmp_df$age <- idx
  tmp_df$rel <- tmp_df$Freq / sum(tmp_df$Freq)
  tmp_df$study <- "LaManno"
  df2plot <- rbind(df2plot, tmp_df)
}

for (idx in c("E10", "E13", "E16", "P0")) {
  tmp_df <- as.data.frame( table( seurat_sketch$cellType[which(seurat_sketch$orig.ident == idx)] ) )
  tmp_df$age <- idx
  tmp_df$rel <- tmp_df$Freq / sum(tmp_df$Freq)
  tmp_df$study <- "this"
  df2plot <- rbind(df2plot, tmp_df)
}

age2idx <- c(1,2,3,4,1,2,3,4,1,2,3,4)
names(age2idx) <- c("E11", "E13", "E16", "P1_S1", "e11.0", "e13.5", "e16.5", "e18.0", "E10", "E13", "E16", "P0")

df2plot$idx <- age2idx[ df2plot$age] 

# calculate an average per developmental age to plot
df2plot %>%
  group_by( idx, Var1 ) %>%
  summarize(mean = mean(rel), sd = sd(rel),
            n = n()) -> embryo_cellTypeAbundanceSummary

# fill missing observations with zeros
for (cellType in unique( embryo_cellTypeAbundanceSummary$Var1 )) {
  for (idx in unique(embryo_cellTypeAbundanceSummary$idx )) {
    tmp <- embryo_cellTypeAbundanceSummary[ which(embryo_cellTypeAbundanceSummary$idx == idx & embryo_cellTypeAbundanceSummary$Var1 == cellType),]
    if (nrow(tmp) == 0) {
      embryo_cellTypeAbundanceSummary <- rbind( embryo_cellTypeAbundanceSummary, data.frame( idx = idx, Var1 = cellType, mean = 0, sd = 0, n=0))
    }
  }
}
colnames( embryo_cellTypeAbundanceSummary ) <- c( "idx", "cell_type", "mean", "sd", "n")

embryo_cellTypeAbundanceSummary$sd[ which( is.na(embryo_cellTypeAbundanceSummary$sd)) ] <- 0

# determine cellType abundances in organoid data at replicate level
cellTypeAbundance <- data.frame()
for (orig_ident in unique(seurat_sketch@meta.data$orig.ident)) {
  for (cell_type in unique(seurat_sketch@meta.data$cellType) ) {
    abundance <- length ( which( seurat_sketch@meta.data$cellType == cell_type & seurat_sketch@meta.data$orig.ident == orig_ident ) )
    tmp_df <- data.frame( abundance = abundance, cell_type = cell_type, group = orig_ident)
    cellTypeAbundance <- rbind( cellTypeAbundance, tmp_df )
  }
}

abundance_vec <- table( seurat_sketch@meta.data$orig.ident )
cellTypeAbundance$tot_cells <- abundance_vec[ cellTypeAbundance$group ]
cellTypeAbundance$rel <- cellTypeAbundance$abundance / cellTypeAbundance$tot_cells
cellTypeAbundance$age <- sapply( cellTypeAbundance$group, function (x) strsplit(x, "-")[[1]][1] )

cellTypeAbundance %>%
  group_by( age, cell_type ) %>%
  summarize(mean = mean(rel), sd = sd(rel),
            n = n()) -> cellTypeAbundanceSummary

age2idx <- c(1,2,3,4,1,2,3,4)
names(age2idx) <- unique(cellTypeAbundance$age)

age2group <- c( "embryo", "embryo", "embryo", "embryo", "EB", "EB", "EB", "EB" )
names(age2group) <- unique(cellTypeAbundance$age)

cellTypeAbundance$idx <- age2idx[ cellTypeAbundance$age ]
cellTypeAbundance$orig <- age2group[ cellTypeAbundance$age ]

cellTypeAbundanceSummary$idx <- age2idx[ cellTypeAbundanceSummary$age ]
cellTypeAbundanceSummary$orig <- age2group[ cellTypeAbundanceSummary$age ]

tmp1 <- cellTypeAbundanceSummary[ which(cellTypeAbundanceSummary$orig == "EB") ,c("idx", "cell_type", "mean", "sd", "orig")]
tmp2 <- embryo_cellTypeAbundanceSummary[ ,c("idx", "cell_type", "mean", "sd")]
tmp2$orig <- "embryo"
comb <- cellTypeAbundanceSummary <- rbind( tmp1, tmp2 )

cellType_order <- c("RGP", "CR", "IP", "iN", "aNSC", "astro", "oligo", "OBNB")
embryo_cellTypeAbundanceSummary$cell_type <- factor( embryo_cellTypeAbundanceSummary$cell_type, levels = cellType_order)
comb$cell_type <- factor( comb$cell_type, levels = cellType_order)
cellTypeAbundance$cell_type <- factor( cellTypeAbundance$cell_type, levels = cellType_order)

# informative label for x-axis
age2idx <- c(1,2,3,4,1,2,3,4,1,2,3,4)
names(age2idx) <- c("E11", "E13", "E16", "P1_S1", "e11.0", "e13.5", "e16.5", "e18.0", "E10", "E13", "E16", "P0")
x_label <- sapply( 1:4, function (i) paste( names(age2idx)[ which(age2idx == i) ], collapse = "," ) )
names(x_label) <- 1:4

ggplot(  ) + 
  geom_ribbon( data = embryo_cellTypeAbundanceSummary, aes(x=idx, ymin=mean-sd, ymax=mean+sd, group = cell_type), fill = "grey80") + 
  geom_line( data = comb, aes(x=idx, y=mean, group = orig, color=orig)) +
  geom_point( data = cellTypeAbundance[ which(cellTypeAbundance$orig == "EB"),], aes(x=idx, y=rel, color=orig)) +
  scale_color_manual( values = c("EB" = "black", "embryo" = "grey50")) +
  scale_x_continuous( label=x_label) +
  facet_wrap( ~cell_type, nrow = 2) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave ( filename = paste(base, "plots/Fig_1/F1i_cellAbundance_main.pdf", sep=""), width=7, height=5 )

colnames(df2plot) <- c("cell_type", "Freq", "age", "rel", "study", "idx")
ggplot(  ) + 
  geom_ribbon( data = embryo_cellTypeAbundanceSummary, aes(x=idx, ymin=mean-sd, ymax=mean+sd, group = cell_type), fill = "grey80") + 
  geom_line( data = comb[which(comb$orig == "embryo"),], aes(x=idx, y=mean, group = orig)) +
  #geom_errorbar( aes(ymin=mean-sd, ymax=mean+sd)) + 
  geom_point( data = df2plot, aes(x=idx, y=rel, color=study)) +
  #scale_color_manual( values = c("EB" = "black", "embryo" = "grey50")) +
  facet_wrap( ~cell_type, nrow = 2) +
  theme_classic()

ggsave ( filename = paste(base, "plots/Fig_1/cellAbundance_embryo_only.pdf", sep=""), width=7, height=5 )


