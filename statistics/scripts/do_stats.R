##############################
###
#    statics for extended data figures as indicated
###
##############################

library (openxlsx)
base <- "~/"

do_pairwise_p <- function(count_matrix, per_col = T) {
  # per_col: matrix is processed per row, all pairwise comparisons of columns are performed
  
  pvalues <- data.frame()
  
  # all pairwise comparisons determined here
  comparisons <- combn(colnames(count_matrix), 2)
  
  # total is for column
  for ( row_idx in 1:nrow(count_matrix)) {
    for (comp_idx in 1:ncol(comparisons)) {
      
      tmp_table <- as.table( rbind( c(count_matrix[row_idx, comparisons[1, comp_idx] ], sum(count_matrix[, comparisons[1, comp_idx] ]) - count_matrix[row_idx, comparisons[1, comp_idx] ]),
                                    c(count_matrix[row_idx, comparisons[2, comp_idx] ], sum(count_matrix[, comparisons[2, comp_idx] ]) - count_matrix[row_idx, comparisons[2, comp_idx] ]) )
      )
      chisq_out <- chisq.test( tmp_table, correct = T, rescale.p = F, simulate.p.value = F )
      fisher_out <- fisher.test( tmp_table )
      pvalues <- rbind( pvalues, data.frame(type = rownames(count_matrix)[row_idx], 
                                            comparison = paste( comparisons[1, comp_idx], comparisons[2, comp_idx], sep="_vs_"), 
                                            chisq_pvalue = chisq_out$p.value , fisher_pvalue = fisher_out$p.value))
      
    }
  }
  
  return( pvalues )
}

save_xlsx <- function( df2save, file_name, sheetName ) {
  
  wb <- createWorkbook()
  addWorksheet(wb, sheetName)
  writeData( wb, sheetName, df2save )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol(df2save) )
  setColWidths( wb, sheetName, cols = 1:ncol(df2save), widths="auto" )
  saveWorkbook(wb, file = file_name, overwrite = T) 
  
}


# extended data figure 10
data <- read.xlsx( xlsxFile = paste(base, "statistics/data/ExtFig10.xlsx", sep=""), rowNames = T, colNames = T)

# pairwise comparison of differentiation time points
test1 <- chisq.test( t(data)[,1:2], correct = T, rescale.p = F, simulate.p.value = F )
test2 <- chisq.test( t(data)[,2:3], correct = T, rescale.p = F, simulate.p.value = F )
test3 <- chisq.test( t(data)[,3:4], correct = T, rescale.p = F, simulate.p.value = F )

pvalue <- c( test1$p.value, test2$p.value, test3$p.value )
padj <- p.adjust ( c( test1$p.value, test2$p.value, test3$p.value ), method = "bonferroni" )
label <- c("D10-D13", "D13-D16", "D16-D20" )
chisq_out <- data.frame( comparison = label, pvalue = pvalue, padj = padj)
save_xlsx( chisq_out, paste(base, "statistics/output/ExtFig10_stats_diffAges.xlsx", sep=""), "chisq")

all_p <- do_pairwise_p ( t(data) )
all_p$chisq_padj <- p.adjust( all_p$chisq_pvalue, method = "bonferroni")
all_p$fisher_padj <- p.adjust( all_p$fisher_pvalue, method = "bonferroni")
all_p <- all_p[ order(all_p$fisher_padj, decreasing = F), ]
save_xlsx( all_p, paste(base, "statistics/output/ExtFig10_stats.xlsx", sep=""), "stats")

data <- read.xlsx( xlsxFile = paste(base, "statistics/data/ExtFig17.xlsx", sep=""), rowNames = F, colNames = T)

all_p <- data.frame()
for (ind in unique(data$induction)) {
  tmp <- data[ which(data$induction == ind), ]
  rownames( tmp ) <- paste( tmp$cloneType, tmp$induction )
  tmp <- tmp[, 3:4]
  all_p <- rbind( all_p, do_pairwise_p ( tmp ) )
}


all_p$chisq_padj <- p.adjust( all_p$chisq_pvalue, method = "bonferroni")
all_p$fisher_padj <- p.adjust( all_p$fisher_pvalue, method = "bonferroni")

save_xlsx( all_p, paste(base, "statistics/output/ExtFig17_stats.xlsx", sep=""), "stats")

data <- read.xlsx( xlsxFile = paste(base, "statistics/data/ExtFig20.xlsx", sep=""), rowNames = T, colNames = T, sheet = 1)
all_p <- do_pairwise_p ( t(data) )
all_p$chisq_padj <- p.adjust( all_p$chisq_pvalue, method = "bonferroni")
all_p$fisher_padj <- p.adjust( all_p$fisher_pvalue, method = "bonferroni")
save_xlsx( all_p, paste(base, "statistics/output/ExtFig20d_stats.xlsx", sep=""), "stats")


data <- read.xlsx( xlsxFile = paste(base, "statistics/data/ExtFig20.xlsx", sep=""), rowNames = T, colNames = T, sheet = 2)
all_p <- do_pairwise_p ( t(data) )
all_p$chisq_padj <- p.adjust( all_p$chisq_pvalue, method = "bonferroni")
all_p$fisher_padj <- p.adjust( all_p$fisher_pvalue, method = "bonferroni")
save_xlsx( all_p, paste(base, "statistics/output/ExtFig20e_stats.xlsx", sep=""), "stats")

# not really related to stats, but part of Supplemental table
# marker gene comparison of DL/UL in organoid and reference

data <- read.xlsx( xlsxFile = paste(base, "/analysis_embryo_comp/plots/Sup_Fig_22/EB_ref_layer_marker_overlap.xlsx", sep=""), sheet = 3) 

tmp_up <- data[ which(data$group == "up_up"), ]
tmp_up <- tmp_up[ order(tmp_up$Ref_FC, tmp_up$Ref_padj, decreasing = T)[1:100], ]
tmp_down <- data[ which(data$group == "down_down"), ]
tmp_down <- tmp_down[ order(tmp_down$Ref_FC, tmp_down$Ref_padj, decreasing = F)[1:100], ]

save_xlsx( tmp_up, paste(base, "statistics/output/UL_markers.xlsx", sep=""), "stats")
save_xlsx( tmp_down, paste(base, "statistics/output/DL_markers.xlsx", sep=""), "stats")

