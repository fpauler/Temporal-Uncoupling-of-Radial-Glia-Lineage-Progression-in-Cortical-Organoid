# prepare a GRanges object with exons and introns and performs read counting
# this is the first step in the PatchSeq analysis

library (GenomicFeatures)
library (GenomicAlignments)
library (scater)
library (BiocParallel)

base <- "./"
no_save_base <- "./"

get_introns <- function (x) {
  
  tmp <- gaps(x)
  # if one intron is present we have 2 gaps
  if (length(tmp) > 1) {
    # Note: first gap is from chromosome start to first exon - discard
    tmp <- tmp[2:length(tmp)]
    return (tmp)
  } else {
   # no intron: return empty granges object
    GRanges(seqnames= NULL, ranges = NULL, strand = NULL)
  }
  
}

if (!file.exists(paste(base, "annotation/gencode.vM27.annotation.GRangesListIntrons.RDS", sep="/"))) {
  
  # read transcript info from gtf file
  txdb <- makeTxDbFromGFF(file = paste(base, "annotation/gencode.vM27.annotation.gtf", sep="/"))
  # create GRanges list with exon info for each gene
  mm.exons <- GenomicFeatures::exonsBy(txdb, by="gene")
  # merge overlapping exons - can happen if many transcript versions are present for one gene
  mm.exons <- reduce(mm.exons)
  
  #create an intronic GRanges List - long wait!
  mm.introns <- endoapply(mm.exons, get_introns)
  
  saveRDS( object =  mm.introns, file=paste(base, "annotation/gencode.vM27.annotation.GRangesListIntrons.RDS", sep="/" ))
  saveRDS( object =  mm.exons, file=paste(base, "annotation/gencode.vM27.annotation.GRangesListExons.RDS", sep="/" ))
  
}
# stop and restart to free memory - only necessary if you have 16GB RAM

#read previously created GRanges objects
mm.introns <- readRDS( file=paste(base, "annotation/gencode.vM27.annotation.GRangesListIntrons.RDS", sep="/" ) )
mm.exons <- readRDS( file=paste(base, "annotation/gencode.vM27.annotation.GRangesListExons.RDS", sep="/" ) )

# read alignment stats - necessary for filtering
# finally prepare a meta data file from the STAR aligment logs
files <- c(list.files(path = no_save_base, pattern=glob2rx("*Log.final.out"), recursive = T, full.names = T))

meta.data <- data.frame()

for (file in files) {
  
  tmp_df <- read.table( file = file, header = F, sep = "\t", fill = T)
  file_name_parts <- strsplit( x = file, split = "/")[[1]]
  sampleID <- gsub(pattern = "STAR\\.", replacement = "", x = file_name_parts[length(file_name_parts)])
  sampleID <- gsub(pattern = "\\.Log\\.final\\.out", replacement = "", x = sampleID)
  
  meta.data <- rbind( meta.data,  data.frame( sampleID = sampleID, total_reads = tmp_df$V2[5], uniquely_aligned = tmp_df$V2[8] ) )
  
}

rownames( meta.data ) <- meta.data$sampleID
meta.data <- meta.data[,c(2,3)]
saveRDS( object = meta.data, file = paste(base, "RDS_files/PatchSeq_meta_data.RDS", sep="/"))

#read all BAM files
files <- c(list.files(path = no_save_base, pattern=glob2rx("*.bam"), recursive = T, full.names = T))

# limit here the number of parallel tasks! Otherwise throws an out of memory error!
#counts reads overlapping exons
message("counting exon reads")
mm.exons.counts <- summarizeOverlaps(features = mm.exons, reads = files,
                          singleEnd=TRUE, mode = "IntersectionNotEmpty", 
                          ignore.strand = T, inter.feature = T, BPPARAM = MulticoreParam( workers = 8 ) )

#count reads overlapping introns
message("counting intron reads") 
mm.intron.counts <- summarizeOverlaps(features = mm.introns, 
                                     reads = files,
                                     singleEnd=TRUE, mode = "IntersectionNotEmpty", 
                                     ignore.strand = T, inter.feature = T, BPPARAM = MulticoreParam( workers = 8 ) )

raw.exon.counts  <- assay(mm.exons.counts)
raw.intron.counts  <- assay(mm.intron.counts)

saveRDS( object = raw.exon.counts, file = paste(base, "RDS_files/raw_exon_counts.RDS", sep="/"))
saveRDS( object = raw.intron.counts, file = paste(base, "RDS_files/raw_intron_counts.RDS", sep="/"))

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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scater_1.30.1               ggplot2_3.5.0               scuttle_1.12.0              SingleCellExperiment_1.24.0 GenomicAlignments_1.38.2   
# [6] Rsamtools_2.18.0            Biostrings_2.70.3           XVector_0.42.0              SummarizedExperiment_1.32.0 MatrixGenerics_1.14.0      
# [11] matrixStats_1.2.0           GenomicFeatures_1.54.4      AnnotationDbi_1.64.1        Biobase_2.62.0              GenomicRanges_1.54.1       
# [16] GenomeInfoDb_1.38.8         IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.48.1        
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.2                 bitops_1.0-7              gridExtra_2.3             biomaRt_2.58.2            rlang_1.1.3              
# [6] magrittr_2.0.3            compiler_4.3.2            RSQLite_2.3.6             DelayedMatrixStats_1.24.0 png_0.1-8                
# [11] vctrs_0.6.5               reshape2_1.4.4            stringr_1.5.1             pkgconfig_2.0.3           crayon_1.5.2             
# [16] fastmap_1.1.1             dbplyr_2.5.0              utf8_1.2.4                ggbeeswarm_0.7.2          bit_4.0.5                
# [21] zlibbioc_1.48.2           cachem_1.0.8              beachmat_2.18.1           progress_1.2.3            blob_1.2.4               
# [26] DelayedArray_0.28.0       BiocParallel_1.36.0       irlba_2.3.5.1             parallel_4.3.2            prettyunits_1.2.0        
# [31] R6_2.5.1                  stringi_1.8.3             rtracklayer_1.62.0        Rcpp_1.0.12               Matrix_1.6-5             
# [36] tidyselect_1.2.1          rstudioapi_0.16.0         abind_1.4-5               yaml_2.3.8                viridis_0.6.5            
# [41] codetools_0.2-19          curl_5.2.0                lattice_0.21-9            tibble_3.2.1              plyr_1.8.9               
# [46] withr_3.0.0               KEGGREST_1.42.0           BiocFileCache_2.10.2      xml2_1.3.6                pillar_1.9.0             
# [51] BiocManager_1.30.22       filelock_1.0.3            packrat_0.9.2             generics_0.1.3            RCurl_1.98-1.14          
# [56] hms_1.1.3                 sparseMatrixStats_1.14.0  munsell_0.5.0             scales_1.3.0              glue_1.7.0               
# [61] tools_4.3.2               BiocIO_1.12.0             BiocNeighbors_1.20.2      ScaledMatrix_1.10.0       XML_3.99-0.16.1          
# [66] grid_4.3.2                colorspace_2.1-0          GenomeInfoDbData_1.2.11   beeswarm_0.4.0            BiocSingular_1.18.0      
# [71] restfulr_0.0.15           vipor_0.4.7               cli_3.6.2                 rsvd_1.0.5                rappdirs_0.3.3           
# [76] fansi_1.0.6               S4Arrays_1.2.1            viridisLite_0.4.2         dplyr_1.1.4               gtable_0.3.4             
# [81] digest_0.6.34             SparseArray_1.2.4         ggrepel_0.9.5             rjson_0.2.23              memoise_2.0.1            
# [86] lifecycle_1.0.4           httr_1.4.7                bit64_4.0.5   
