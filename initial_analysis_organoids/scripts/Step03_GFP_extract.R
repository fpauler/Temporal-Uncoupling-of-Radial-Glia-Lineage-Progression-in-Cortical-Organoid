# ========================================================================
# Step 3 
# extract GFP cells, remove stressed cells, add cell cycle score
# ========================================================================

library (Seurat)
library (pheatmap)
library (ggplot2)
library (SeuratObject)
library (gruffi)
library (GO.db)
library (org.Mm.eg.db)
library (clusterProfiler)
library (SeuratWrappers)

set.seed(2401)

base <- "~/initial_analysis_organoids/"
nobackup_base <- "~/noSave/"

if ( !file.exists(paste( nobackup_base, "RDS_files/seurat_obj.GFP.merged.RDS", sep = "")) ) {
  
  GFPAbundance <- data.frame()
  
  # define clusters to include
  # comes from the p-value heatmap prepared below
  clusters2include <- list()
  clusters2include[["D8"]] <- c("2", "3", "9") 
  clusters2include[["D13"]] <- c("0", "1", "2", "3", "4", "5")
  clusters2include[["D20"]] <- c("0", "1", "2", "3", "4", "6", "8", "9", "10", "11")
  clusters2include[["D25"]] <- c("0", "1", "2", "3", "4", "6", "7", "8", "10", "11",
                                 "12", "13", "15" )
    
  
  seurat_obj_list <- list()
  
  for ( ID in c("D8", "D13", "D20", "D25" ) ) {
    message( ID )
    sig_mat <- readRDS ( file = paste( nobackup_base, "RDS_files/", ID, ".GFP.pvalue.mar.RDS", sep = "") )
  
    sig_mat <- log10(sig_mat) *-1
    rownames(sig_mat) <- as.character( 0 : ( ncol(sig_mat) - 1 ) )
    colnames(sig_mat) <- as.character( 0 : ( ncol(sig_mat) - 1 ) )
    
    # exchange infinite value with some reasonable value to plot on the heatmap
    sig_mat[ which(is.infinite(sig_mat)) ] <- max (sig_mat[ which(!is.infinite(sig_mat)) ])
    pheatmap( sig_mat, scale="none", cluster_rows = T, cluster_cols = T, 
              filename = paste( base, "plots/QC/", ID,".eGFP.pvalue.pdf", sep="" ), main = paste( ID, "p-value plot" ))
    
    seurat.obj <- readRDS( file =  paste( nobackup_base,  "/RDS_files/seurat.", ID, ".RDS", sep="" ) )
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:25)
    seurat.obj <- FindClusters(seurat.obj, resolution = 0.3)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:25)
    
    # add information about GFP expression for supplemental figure
    seurat.obj@meta.data$GFP_cluster <- "NO"
    seurat.obj@meta.data$GFP_cluster[ which( as.character(seurat.obj@meta.data$seurat_clusters) %in% clusters2include[[ ID ]] ) ] <- "YES"
    
    # prepare a count table for fraction of GFP expressing clusters
    rep_abundance <- table (seurat.obj$orig.ident)
    tmp_df_all <- data.frame()
    for (ident in unique(seurat.obj$orig.ident) ) {
      tmp_df <- as.data.frame( table( seurat.obj$GFP_cluster[ which( seurat.obj$orig.ident == ident)]) )
      tmp_df$rep <- ident
      tmp_df$age <- ID
      tmp_df_all <- rbind( tmp_df_all, tmp_df )
    }
    
    tmp_df_all$total <- rep_abundance[ tmp_df_all$rep ]
    GFPAbundance <- rbind( GFPAbundance, tmp_df_all )
    
    seurat.obj <- subset( seurat.obj, idents = clusters2include[[ ID ]] )
    
    umap_plot <- DimPlot( seurat.obj, reduction = "umap", group.by = "seurat_clusters", label = T) + NoLegend()
    DimPlot( seurat.obj, reduction = "umap", group.by = "orig.ident", label = T) + NoLegend()
    ggsave( plot = umap_plot, file = paste ( base, "plots/QC/", ID, ".GFP.extracted.pdf", sep="" ), width = 10, height = 10 )
  
    total_number <- table(  seurat.obj@meta.data$orig.ident )
    
    all_abundance_df <- data.frame()
    for( name in unique(seurat.obj@meta.data$seurat_clusters) ) {
      tmp_df <- as.data.frame( table( seurat.obj@meta.data[ which( seurat.obj@meta.data$seurat_clusters == name ), "orig.ident"] ) )
      tmp_df$age <- sapply( as.character(tmp_df$Var1), function (x) strsplit(x, "-")[[1]][1])
      tmp_df$total <- total_number[ as.character(tmp_df$Var1) ]
      tmp_df$rel <- 100 * tmp_df$Freq / tmp_df$total
      tmp_df$group <- name
      all_abundance_df <- rbind( all_abundance_df, tmp_df)
    }
    
    # plot the number of cells in the GFP clusters
    abundance_plot <- ggplot(all_abundance_df, aes(x=group, y=rel, fill=Var1)) + 
      geom_bar( stat = "identity", position = "dodge" )+ theme_classic() +
      theme(text = element_text(size = 15), axis.text = element_text(size = 15))
    
    ggsave( plot = abundance_plot, file = paste ( base, "plots/QC/", ID, ".GFP.extracted.abundances.pdf", sep="" ) )
    
    # remove unnecessary data from object to make it smaller
    for (layerID in c("lineA1", "lineA2", "lineB1", "lineB2")) {
      LayerData( seurat.obj, paste("data.", ID, "-lineA1", sep="") ) <- NULL
    }
    LayerData( seurat.obj, "scale.data" ) <- NULL
    LayerData( seurat.obj, "pca" ) <- NULL
    
    seurat_obj_list[[ length(seurat_obj_list)+1 ]] <- JoinLayers( seurat.obj )
  }
  
  GFPAbundance$replicate <- sapply( GFPAbundance$rep, function (x) strsplit( as.character( x ), "-")[[1]][2] )
  GFPAbundance$replicate <- gsub( pattern = "line", replacement = "", x = GFPAbundance$replicate )
  
  # [[1]]
  # An object of class Seurat 
  # 16585 features across 10006 samples within 1 assay 
  # Active assay: RNA (16585 features, 2000 variable features)
  # 2 layers present: data, counts
  # 2 dimensional reductions calculated: pca, umap
  # 
  # [[2]]
  # An object of class Seurat 
  # 15787 features across 35374 samples within 1 assay 
  # Active assay: RNA (15787 features, 2000 variable features)
  # 2 layers present: data, counts
  # 2 dimensional reductions calculated: pca, umap
  # 
  # [[3]]
  # An object of class Seurat 
  # 15893 features across 33369 samples within 1 assay 
  # Active assay: RNA (15893 features, 2000 variable features)
  # 2 layers present: data, counts
  # 2 dimensional reductions calculated: pca, umap
  # 
  # [[4]]
  # An object of class Seurat 
  # 16396 features across 32084 samples within 1 assay 
  # Active assay: RNA (16396 features, 2000 variable features)
  # 2 layers present: data, counts
  # 2 dimensional reductions calculated: pca, umap
  
  seurat.obj <- merge( x = seurat_obj_list[[ 1 ]], y = c(seurat_obj_list[[ 2 ]], seurat_obj_list[[ 3 ]], seurat_obj_list[[ 4 ]]))
  
  saveRDS ( seurat.obj, file = paste( nobackup_base, "RDS_files/seurat_obj.GFP.merged.RDS", sep = "") )
  
  stop( "restart R to free memory and re-run the script" )
}

if (!file.exists(paste( nobackup_base, "RDS_files/seurat_obj.stress.removed.RDS", sep = ""))) {
  
  seurat.obj <- readRDS ( file = paste( nobackup_base, "RDS_files/seurat_obj.GFP.merged.RDS", sep = "") )
  
  # redo standard analysis
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 4000)
  
  seurat.obj <- ScaleData(seurat.obj)
  
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  
  # PC_ 1 
  # Positive:  Mapt, Crmp1, Stmn2, Dcx, Rtn1, Dpysl3, Mllt11, Stmn4, Ina, Nsg2 
  # Aplp1, Islr2, Elavl3, Gap43, Cdk5r1, Tubb3, Nsg1, L1cam, Gpm6a, Stxbp1 
  # Gria2, Ncam1, Myt1l, Ctnna2, Sez6l2, Jph4, Chl1, Tuba1a, Soga3, Tubb2b 
  # Negative:  Vim, H2afv, Sox2, Tead2, Qk, Ccna2, Notch1, Mki67, Hes5, Nt5dc2 
  # Top2a, Nr2e1, Gli3, Hist1h1e, Rrm2, Cdca3, Ccnd2, Pax6, Ctnna1, Uhrf1 
  # H2afz, Cdk1, Mcm2, Fzd2, Smc4, Tacc3, Mdk, Hist1h1b, Rrm1, Cdca8 
  # PC_ 2 
  # Positive:  Stmn1, Sox11, Jpt1, Map1b, Igfbpl1, H1fx, Top2a, Ccna2, Mki67, Nusap1 
  # Cdk1, Cdca3, Prc1, Pimreg, Ube2c, Kif11, Calm2, Tpx2, Racgap1, Ect2 
  # Ncapd2, Cdca2, Spc25, Pbk, Ncapg, H2afx, Plk1, Cdca8, Hist1h2ap, Bub1 
  # Negative:  Bcan, Cspg5, Atp1a2, Luzp2, Mfge8, Ndrg2, Ptprz1, Slc1a3, Zcchc24, Apoe 
  # Slc38a3, Ednrb, Serpine2, S100a16, Fabp7, Ncald, Smoc1, Tril, Asrgl1, Pla2g7 
  # Cmtm5, Olig1, Olig2, Abca1, Ptn, Mlc1, Mt3, Sox8, Ntrk2, Ddah1 
  # PC_ 3 
  # Positive:  Cntn1, Lhfpl3, Ncald, Prkcb, Olig1, Pdgfra, Olig2, Cdk1, Top2a, Nusap1 
  # Kif11, Serpine2, Pcdh15, Opcml, Tnr, Prc1, Mdga2, Sox10, Cspg5, Ckap2l 
  # Plppr5, Incenp, Mki67, Ccna2, Neto1, Omg, Kifc1, Nuf2, Kif18b, Gatm 
  # Negative:  Igfbpl1, Nhlh1, Eomes, Sox11, Mfng, Rcor2, Celsr1, Scrt2, Sstr2, Neurod1 
  # Rnd2, Ddit4, Mfap4, Elavl2, Nkd1, Shb, Phf21b, Abcd2, Cacna2d1, Neurog2 
  # Meis2, Sorbs2, Nxph4, Epha3, Srrm4, Kdm7a, Ppp1r14a, Gse1, Coro2b, Emx1 
  # PC_ 4 
  # Positive:  Lhfpl3, Olig1, Pdgfra, Olig2, Myt1, Sox10, Dll3, Gpr17, Pcdh15, Prkcq 
  # Neto1, Omg, Spon1, Elfn1, Enpp2, Nkx2-2, Brinp1, Dscam, Dock10, Ncald 
  # Tnr, Ets1, Nxph1, Nova1, Ppp1r14b, Cntn1, Serpine2, Pllp, Lims2, Matn4 
  # Negative:  Igfbp5, Apoe, Clu, Mlc1, Pygb, Naaa, S1pr1, Tspan18, Rgma, Tnc 
  # Sparc, Ttyh1, Jun, Pla2g7, Sparcl1, Pdpn, Itgb5, Gdpd2, Cyr61, Slc9a3r1 
  # Serpinh1, Vcam1, Aldoc, Gja1, Fgfr1, Gstm1, Emp2, Cyp46a1, Mt3, Ccdc80 
  # PC_ 5 
  # Positive:  Nfia, Scd2, GFP, Nfib, Fads2, Mmd2, Cdkn2c, Hmgcs1, Itm2b, Hmgcr 
  # Mlc1, Slc1a3, Atp1a2, Tnc, Slc15a2, Pla2g7, Scd1, Msmo1, Troap, Cyp51 
  # Nwd1, Prc1, Fdps, Naaa, Mfge8, Ckap2l, Celsr1, Cd9, Cst3, Ube2c 
  # Negative:  Crabp2, Igdcc3, Zbtb16, Prtg, Hmga2, Otx2, Fezf1, Loxl2, Rbp1, Lin28b 
  # Metrn, Lin28a, Fndc3c1, Mest, Dlk1, Scube3, Ctgf, Mecom, Slc16a3, Car14 
  # Slc25a13, Id3, Cdh22, Thbs1, Ifitm2, Slc2a3, Peg10, Slc2a1, Dtx4, Ldha 
  
  # confirm low at PC20
  elbow_plot <- ElbowPlot(seurat.obj)
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:20, reduction = "pca")
  
  # identify and remove stressed cells
  
  # first identify GO terms to prepare gene lists
  
  # which ( grepl( pattern = "glycolysis", x = Term ( GOTERM ) ) )
  Term ( GOTERM )[ which ( grepl( pattern = "glycolysis", x = Term ( GOTERM ) ) ) ]
  Term( GOTERM )[27182]
  
  # GO:0061621 
  # "canonical glycolysis" 
  
  Term ( GOTERM )[ which ( grepl( pattern = "reticulum", x = Term ( GOTERM ) ) & grepl( pattern = "stress", x = Term ( GOTERM ) ) ) ]
  
  # first 2 entries of the list
  
  # GO:0034976 
  # "response to endoplasmic reticulum stress" 
  # GO:0036483 
  # "neuron intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress" 
  
  Term ( GOTERM )[ which ( grepl( pattern = "gliogenesis", x = Term ( GOTERM ) ) ) ]
  
  # GO:0014013                           GO:0014014                           GO:0014015 
  # "regulation of gliogenesis" "negative regulation of gliogenesis" "positive regulation of gliogenesis" 
  # GO:0042063 
  # "gliogenesis" 
  
  Term ( GOTERM )[ which ( grepl( pattern = "forebrain generation of neurons", x = Term ( GOTERM ) ) ) ]
  # GO:0021872 
  # "forebrain generation of neurons"  
  
  Term ( GOTERM )[ which ( grepl( pattern = "forebrain neuroblast division", x = Term ( GOTERM ) ) ) ]
  # first entry of the list
  # GO:0021873 
  # "forebrain neuroblast division" 
  
  all_genes = keys(org.Mm.eg.db, keytype = "ENTREZID")
  tb = AnnotationDbi::select(org.Mm.eg.db, keys = all_genes, keytype = "ENTREZID", 
              columns = c("GOALL", "ONTOLOGYALL"))
  
  glycolysis_gl <- tb[ which( tb$GOALL == "GO:0061621" ), "ENTREZID"]
  er_stress_gl <- tb[ which( tb$GOALL == "GO:0034976" ), "ENTREZID"]
  gliognesis_gl <- tb[ which( tb$GOALL == "GO:0042063" ), "ENTREZID"]
  rgp_gl <- tb[ which( tb$GOALL == "GO:0021872" ), "ENTREZID"]
  rgp_gl <- unique( c( rgp_gl, tb[ which( tb$GOALL == "GO:0021873" ), "ENTREZID"] ) )
  
  glycolysis_gl <- bitr( geneID = glycolysis_gl, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db )
  er_stress_gl <- bitr( geneID = er_stress_gl, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db )
  gliognesis_gl <- bitr( geneID = gliognesis_gl, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db )
  rgp_gl <- bitr( geneID = rgp_gl, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db )
  
  # cluster identification is necessary for gruffi to run
  seurat.obj <- FindNeighbors(seurat.obj, reduction = "pca", dims = 1:30)
  seurat.obj <- FindClusters(seurat.obj, resolution = 2)
  
  # gruffi worflow from here onwards
  seurat.obj <- AutoFindGranuleResolution(obj = seurat.obj)
  
  (granule.res.4.gruffi <- GetGruffiClusteringName(seurat.obj)) # Recalled from @misc$gruffi$'optimal.granule.res.
  # character(0)
  # [1] "RNA_snn_res.61"
  # 
  seurat.obj <- ReassignSmallClusters(seurat.obj, ident = granule.res.4.gruffi) # Will be stored in meta data column as "seurat_clusters.reassigned".
  
  (granule.res.4.gruffi <- GetGruffiClusteringName(seurat.obj))
  # "RNA_snn_res.61.reassigned"
  # [1] "RNA_snn_res.61.reassigned"
  
  # add genes to relevant slots for Gruffi to use
  seurat.obj@misc$gruffi$GO <-list()
  seurat.obj@misc$gruffi$GO$GO.glycolysis <- glycolysis_gl$SYMBOL
  seurat.obj@misc$gruffi$GO$GO.er_stress <- er_stress_gl$SYMBOL
  seurat.obj@misc$gruffi$GO$GO.gliogenesis <- gliognesis_gl$SYMBOL
  seurat.obj@misc$gruffi$GO$GO.RGPs <- rgp_gl$SYMBOL
  
  seurat.obj <- JoinLayers( seurat.obj )
  
  seurat.obj <- AssignGranuleAverageScoresFromGOterm(obj = seurat.obj, GO_term = "GO:glycolysis", 
                                                     save.UMAP = F, new_GO_term_computation = T, 
                                                     clustering = granule.res.4.gruffi, plot.each.gene = F)
  
  seurat.obj <- AssignGranuleAverageScoresFromGOterm(obj = seurat.obj, GO_term = "GO:er_stress", 
                                                     save.UMAP = TRUE, new_GO_term_computation = T, 
                                                     clustering = granule.res.4.gruffi, plot.each.gene = F)
  
  seurat.obj <- AssignGranuleAverageScoresFromGOterm(obj = seurat.obj, GO_term = "GO:gliogenesis", 
                                                     save.UMAP = TRUE, new_GO_term_computation = T, 
                                                     clustering = granule.res.4.gruffi, plot.each.gene = F)
  
  seurat.obj <- AssignGranuleAverageScoresFromGOterm(obj = seurat.obj, GO_term = "GO:RGPs", 
                                                     save.UMAP = TRUE, new_GO_term_computation = T, 
                                                     clustering = granule.res.4.gruffi, plot.each.gene = F)
  
  # Create score names:
  (i1 <- ParseGruffiGranuleScoreName(obj = seurat.obj, goID = "GO:glycolysis"))
  # "RNA_snn_res.61.reassigned"
  # [1] "RNA_snn_res.61.reassigned_cl.av_GO.glycolysis"

  (i2 <- ParseGruffiGranuleScoreName(obj = seurat.obj, goID = "GO:er_stress"))
  # "RNA_snn_res.61.reassigned"
  # [1] "RNA_snn_res.61.reassigned_cl.av_GO.er_stress"
  
  (i3 <- ParseGruffiGranuleScoreName(obj = seurat.obj, goID = "GO:gliogenesis"))
  # "RNA_snn_res.61.reassigned"
  # [1] "RNA_snn_res.61.reassigned_cl.av_GO.gliogenesis"
  #
  (i4 <- ParseGruffiGranuleScoreName(obj = seurat.obj, goID = "GO:RGPs"))
  # "RNA_snn_res.61.reassigned"
  # [1] "RNA_snn_res.61.reassigned_cl.av_GO.RGPs"

  # Call Shiny app
  seurat.obj <- FindThresholdsShiny(obj = seurat.obj,
                                     stress.ident1 = i1,
                                     stress.ident2 = i2,
                                     notstress.ident3 = i3,
                                     notstress.ident4 = i4 )
  
  # set threshold manually
  # seurat.obj@misc$gruffi$thresh.stress.ident1 <- 3.586
  # seurat.obj@misc$gruffi$thresh.stress.ident2 <- 0.384
  # seurat.obj@misc$gruffi$thresh.notstress.ident3 <- 0.647
  # seurat.obj@misc$gruffi$thresh.notstress.ident4 <- -0.001
  
  "Dont forget to click the button in the app: Save New Thresholds"
  
  StressUMAP( obj = seurat.obj, GO.terms = c(i1, i2, i3, i4), reduction = "umap" )
  ggsave( file = paste ( base, "plots/QC/StressUMAP.png", sep="" ), width = 6, height = 5 )
  
  length( which( !seurat.obj@meta.data$is.Stressed ) )
  # [1] 103252
  
  stress_df <- as.data.frame ( table( seurat.obj@meta.data$orig.ident[ which( seurat.obj@meta.data$is.Stressed ) ] ) )
  stress_df$age <- sapply( as.character( stress_df$Var1 ), function (x) strsplit( x = x, split = "-")[[1]][1] )
  
  stress_df$age <- gsub( pattern = "D8", replacement = "D08", x = stress_df$age )
  stress_df <- stress_df[ order ( stress_df$age, stress_df$Freq ), ]
  stress_df$Var1 <- factor( stress_df$Var1, levels = stress_df$Var1 )
  
  tot_cells <- table( seurat.obj@meta.data$orig.ident )
  
  stress_df$total_cells <- tot_cells[ as.character(stress_df$Var1) ]
  stress_df$frac <- 100 * stress_df$Freq / stress_df$total_cells
  stress_df$rep <- sapply( stress_df$Var1, function (x) strsplit(x = as.character(x), split = "-")[[1]][2])
  stress_df$rep <- gsub( pattern = "line", replacement = "", x = stress_df$rep )
  
  stress_df$age <- sapply( stress_df$Var1, function (x) strsplit(x = as.character(x), split = "-")[[1]][1])
  stress_df$age <- factor(stress_df$age, levels = c("D8", "D13", "D20", "D25"))
  
  write.csv( x = stress_df, file = paste ( base, "plots/Supplement/fraction_stressed.csv", sep="" ))
  
  stressed_frac <- ggplot( stress_df, aes( x = rep, y=frac ) ) + geom_bar(stat="identity", position="dodge") + 
    ylim(0,100) + theme_classic() + facet_grid(~age)
  
  ggsave( plot = stressed_frac, file = paste ( base, "plots/Supplement/fraction_stressed.pdf", sep="" ), width = 6, height = 5 )
  
  Idents( seurat.obj ) <- "is.Stressed"
  seurat.obj <- subset( seurat.obj, ident = FALSE)
  
  # remove normalized and scaled data
  seurat.obj@assays$RNA$data <- NULL
  seurat.obj@assays$RNA$scale.data <- NULL
  # create an age column
  seurat.obj@meta.data$age <- sapply( seurat.obj@meta.data$orig.ident, function (x) strsplit(x, "-")[[1]][1])
  #split the object based on age - necessary for later integration
  seurat.obj[["RNA"]] <- split( seurat.obj[["RNA"]], f = seurat.obj$age )
  
  seurat.obj@misc$gruffi
  # $thresh.stress.ident1
  # [1] 3.586
  # 
  # $thresh.stress.ident2
  # [1] 0.384
  # 
  # $thresh.notstress.ident3
  # [1] 0.647
  # 
  # $thresh.notstress.ident4
  # [1] -0.001
  
  saveRDS( seurat.obj, paste( nobackup_base, "RDS_files/seurat_obj.stress.removed.RDS", sep = "") )
  stop ("restart R to free memory and re-run the script")
}

if ( !file.exists(paste( nobackup_base, "RDS_files/seurat_obj.stress.CC_added.metaData.RDS", sep = ""))) {
  
  # determine cell cycle phase
  seurat.obj <- readRDS( paste( nobackup_base, "RDS_files/seurat_obj.stress.removed.RDS", sep = "") )
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 4000)

  cell_cycle_phase_gene <- read.table(file = paste(base, "other_data/Mus_musculus.csv", sep=""), sep = ",", header = T)
  gene_conversion <- bitr(geneID = cell_cycle_phase_gene$geneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
  cell_cycle_phase_gene <- merge(cell_cycle_phase_gene, gene_conversion, by.x="geneID", by.y="ENSEMBL")
  
  # CellCycleScoring only works with joint layers
  seurat.obj <- JoinLayers( seurat.obj )
  
  seurat.obj <- CellCycleScoring (
    object = seurat.obj,
    g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
    s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
  )
  
  saveRDS ( object = seurat.obj@meta.data, file = paste( nobackup_base, "RDS_files/seurat_obj.stress.CC_added.metaData.RDS", sep = "") )
  stop ("restart R to free memory and re-run the script")
}

sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.4 LTS
# 
# Matrix products: default
# BLAS:   /opt/R/4.3.2/lib/R/lib/libRblas.so 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] SeuratWrappers_0.3.5   clusterProfiler_4.10.1 org.Mm.eg.db_3.18.0    GO.db_3.18.0           AnnotationDbi_1.64.1   IRanges_2.36.0         S4Vectors_0.40.2      
# [8] Biobase_2.62.0         BiocGenerics_0.48.1    gruffi_1.5.5           MarkdownReports_4.7.0  MarkdownHelpers_1.0.7  CodeAndRoll2_2.6.0     Stringendo_0.6.0      
# [15] magrittr_2.0.3         ggplot2_3.5.0          pheatmap_1.0.12        Seurat_5.0.1           SeuratObject_5.0.1     sp_2.1-3              
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.3                 matrixStats_1.2.0        spatstat.sparse_3.0-3    bitops_1.0-7             enrichplot_1.22.0        HDO.db_0.99.1           
# [7] httr_1.4.7               RColorBrewer_1.1-3       tools_4.3.2              sctransform_0.4.1        backports_1.4.1          utf8_1.2.4              
# [13] R6_2.5.1                 sm_2.2-6.0               lazyeval_0.2.2           uwot_0.1.16              withr_3.0.0              prettyunits_1.2.0       
# [19] gridExtra_2.3            tictoc_1.2.1             VennDiagram_1.7.3        progressr_0.14.0         cli_3.6.2                formatR_1.14            
# [25] spatstat.explore_3.2-6   fastDummies_1.7.3        scatterpie_0.2.2         ggExpress_0.9.0          spatstat.data_3.0-4      readr_2.1.5             
# [31] ggridges_0.5.6           pbapply_1.7-2            yulab.utils_0.1.4        gson_0.1.0               DOSE_3.28.2              colorRamps_2.3.4        
# [37] R.utils_2.12.3           vioplot_0.4.0            parallelly_1.37.0        sessioninfo_1.2.2        rstudioapi_0.16.0        RSQLite_2.3.6           
# [43] gridGraphics_0.5-1       generics_0.1.3           ggVennDiagram_1.5.2      RApiSerialize_0.1.2      gtools_3.9.5             ica_1.0-3               
# [49] spatstat.random_3.2-2    vroom_1.6.5              dplyr_1.1.4              zip_2.3.1                Matrix_1.6-5             futile.logger_1.4.3     
# [55] fansi_1.0.6              clipr_0.8.0              abind_1.4-5              R.methodsS3_1.8.2        terra_1.7-71             lifecycle_1.0.4         
# [61] SoupX_1.6.2              qvalue_2.34.0            gplots_3.1.3.1           BiocFileCache_2.10.2     Rtsne_0.17               grid_4.3.2              
# [67] blob_1.2.4               promises_1.2.1           crayon_1.5.2             miniUI_0.1.1.1           lattice_0.21-9           cowplot_1.1.3           
# [73] KEGGREST_1.42.0          pillar_1.9.0             knitr_1.45               fgsea_1.28.0             future.apply_1.11.1      codetools_0.2-19        
# [79] fastmatch_1.1-4          leiden_0.4.3.1           glue_1.7.0               packrat_0.9.2            ggfun_0.1.4              remotes_2.5.0           
# [85] data.table_1.15.0        treeio_1.26.0            vctrs_0.6.5              png_0.1-8                spam_2.10-0              gtable_0.3.4            
# [91] cachem_1.0.8             xfun_0.42                openxlsx_4.2.5.2         princurve_2.1.6          mime_0.12                tidygraph_1.3.1         
# [97] ReadWriter_1.5.3         tidyverse_2.0.0          survival_3.5-7           iterators_1.0.14         rgl_1.3.1                ellipsis_0.3.2          
# [103] fitdistrplus_1.1-11      ROCR_1.0-11              nlme_3.1-163             ggtree_3.10.1            Seurat.utils_2.7.3       bit64_4.0.5             
# [109] progress_1.2.3           filelock_1.0.3           RcppAnnoy_0.0.22         GenomeInfoDb_1.38.8      job_0.3.0                irlba_2.3.5.1           
# [115] KernSmooth_2.23-22       colorspace_2.1-0         DBI_1.2.2                raster_3.6-26            tidyselect_1.2.1         bit_4.0.5               
# [121] compiler_4.3.2           curl_5.2.0               xml2_1.3.6               plotly_4.10.4            shadowtext_0.1.3         stringfish_0.16.0       
# [127] checkmate_2.3.1          scales_1.3.0             caTools_1.18.2           lmtest_0.9-40            rappdirs_0.3.3           stringr_1.5.1           
# [133] digest_0.6.34            goftest_1.2-3            spatstat.utils_3.0-4     XVector_0.42.0           base64enc_0.1-3          htmltools_0.5.7         
# [139] pkgconfig_2.0.3          sparseMatrixStats_1.14.0 MatrixGenerics_1.14.0    dbplyr_2.5.0             fastmap_1.1.1            rlang_1.1.3             
# [145] htmlwidgets_1.6.4        shiny_1.8.0              farver_2.1.1             zoo_1.8-12               jsonlite_1.8.8           BiocParallel_1.36.0     
# [151] GOSemSim_2.28.1          R.oo_1.26.0              RCurl_1.98-1.14          ggplotify_0.1.2          GenomeInfoDbData_1.2.11  dotCall64_1.1-1         
# [157] patchwork_1.2.0          munsell_0.5.0            Rcpp_1.0.12              ape_5.8                  viridis_0.6.5            reticulate_1.35.0       
# [163] stringi_1.8.3            ggraph_2.2.1             zlibbioc_1.48.2          MASS_7.3-60              plyr_1.8.9               org.Hs.eg.db_3.18.0     
# [169] parallel_4.3.2           listenv_0.9.1            ggrepel_0.9.5            deldir_2.0-2             graphlayouts_1.1.1       Biostrings_2.70.3       
# [175] splines_4.3.2            tensor_1.5               hms_1.1.3                igraph_2.0.2             spatstat.geom_3.2-8      RcppHNSW_0.6.0          
# [181] reshape2_1.4.4           biomaRt_2.58.2           futile.options_1.0.1     XML_3.99-0.16.1          BiocManager_1.30.22      RcppParallel_5.1.7      
# [187] lambda.r_1.2.4           tweenr_2.0.3             tzdb_0.4.0               foreach_1.5.2            httpuv_1.6.14            RANN_2.6.1              
# [193] tidyr_1.3.1              purrr_1.0.2              polyclip_1.10-6          qs_0.26.1                future_1.33.1            scattermore_1.2         
# [199] ggforce_0.4.2            rsvd_1.0.5               xtable_1.8-4             tidytree_0.4.6           RSpectra_0.16-1          ggcorrplot_0.1.4.1      
# [205] later_1.3.2              viridisLite_0.4.2        tibble_3.2.1             aplot_0.2.2              memoise_2.0.1            cluster_2.1.4           
# [211] HGNChelper_0.8.1         globals_0.16.2           DatabaseLinke.R_1.7.0 

