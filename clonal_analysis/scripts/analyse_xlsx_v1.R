# ========================================================================
# analysis of clone count data - figures that show the data are indicated
# ========================================================================

library ( tidyxl )
library ( openxlsx )
library ( ggbeeswarm )
library ( dplyr )
library ( cowplot )

# ================================================
# collect and plot batch resolved data
# ================================================

# ---------------------------------------
# define a base folder and read the data
# ---------------------------------------
base <- "~/"
sheets <- getSheetNames( file = paste( base, "/clonal_analysis/data/Clones_Cell_Lines_OAM16.01.2026updated.xlsx", sep="" ) )

# ---------------------------------------------------------------------------------
# collect the data from different excel sheets
# batch data is color coded, which is converted to a suitable form for analysis here
# ---------------------------------------------------------------------------------

combined_clone_df <- data.frame()
for (sheet in sheets) {
  
  clone_df <- read.xlsx( paste( base, "/clonal_analysis/data/Clones_Cell_Lines_OAM16.01.2026updated.xlsx", sep="" ), 
                         sheet = sheet, skipEmptyRows = F, colNames = F)
  
  clone_df <- clone_df[1:24,]
  lines <- unique( clone_df[1, which(!is.na(clone_df[1,]))] )
  lines <- lines[ 1:(length(lines)-1) ]
  
  for (i in 0:(length(lines)-1)) {
    for (x in 1:3) {
      col_idx <- x+(3*i)
      cloneType <- clone_df[2, col_idx]
      counts <- clone_df[3:24, col_idx]
      counts <- counts[ which(!is.na(counts))]
      if (length(counts) > 0) {
        row_id <- 3:(length(counts)+2)
        col_id <- LETTERS[ col_idx ]
        cellID <- paste(col_id, row_id, sep="")
        df_row <- data.frame(age = sheet, cloneType = cloneType, counts = counts, cellID = cellID)
        combined_clone_df <- rbind(combined_clone_df, df_row)
      }
    }
  }
  
}

xlsx_cell_info <- xlsx_cells( paste( base, "/clonal_analysis/data/Clones_Cell_Lines_OAM16.01.2026updated.xlsx", sep="" ) )
xlsx_cell_format <- xlsx_formats( paste( base, "/clonal_analysis/data/Clones_Cell_Lines_OAM16.01.2026updated.xlsx", sep="" ) )

# prepare a conversion table
conv_vec <- c()
form_id_vec <- c()
for (batch in c("Batch 1", "Batch 2", "Batch 3", "Batch 4")) {
  tmp_form_id <- unique( xlsx_cell_info$local_format_id[ which( xlsx_cell_info$character == batch) ] )
  colors <- unique( xlsx_cell_format$local$fill$patternFill$fgColor$rgb[ tmp_form_id ] )
  conv_vec <- c( conv_vec, rep(batch, length(colors)) )
  form_id_vec <- c(form_id_vec, colors)
}
names(conv_vec) <- form_id_vec

batch_vec <- c()
for(x in 1:nrow(combined_clone_df)) {
  idx <- which( xlsx_cell_info$sheet == combined_clone_df$age[x] & xlsx_cell_info$address == combined_clone_df$cellID[x])
  form_id <- xlsx_cell_info$local_format_id [ idx ]
  color <- xlsx_cell_format$local$fill$patternFill$fgColor$rgb[ form_id ]
  if (is.na(color)) {
    batch_vec <- c( batch_vec, "x" )
  } else {
    batch_vec <- c( batch_vec, conv_vec[ color ] )
  }
  
}
combined_clone_df$batch <- batch_vec

combined_clone_df <- combined_clone_df[ which(grepl(pattern = "Batch", x = combined_clone_df$batch)),]

combined_clone_df$counts <- as.numeric( combined_clone_df$counts )

# ------------------------------------------
# create a more concise table for reporting
# ------------------------------------------

merged_df <- combined_clone_df[ which(combined_clone_df$age == "D8"), c("age", "batch", "cloneType", "counts") ]
merged_df <- merged_df[ order(merged_df$cloneType, merged_df$batch), ]
merged_df$idx <- 1:nrow(merged_df)

for (age in c("D9", "D10", "D11", "D12", "D13", "D15")) {
  tmp2 <- combined_clone_df[ which(combined_clone_df$age == age), c("age", "batch", "cloneType", "counts") ]
  tmp2 <- tmp2[ order(tmp2$cloneType, tmp2$batch), ]
  tmp2$idx <- 1:nrow(tmp2)
  
  merged_df <- merge( merged_df, tmp2, by="idx", all=T)
  
}

write.csv( x = merged_df, file = paste( base, "/clonal_analysis/output/Ext_Data_Fig16_clone_counts_comb.csv", sep="" ), row.names = T )
# Note: this file was the basis for manually curated Ext_Data_Fig16_clone_counts_comb.xlsx


# ------------------------------------------------------------
# summarize and prepare the clone data for plot and trendline
# ------------------------------------------------------------

combined_clone_df %>%
  group_by( cloneType, age ) %>%
  summarise( mean = mean(counts), sd =  sd(counts), n = n(), se = sd(counts)/sqrt(n())) -> summary_clones
  
order_age <- c("D8", "D9", "D10", "D11", "D12", "D13", "D15")
summary_clones$age <- factor(summary_clones$age, levels = order_age)
combined_clone_df$age <- factor(combined_clone_df$age, levels = order_age)

# numeric age
combined_clone_df$num_age <- as.numeric( gsub( pattern = "D", replacement = "", x = combined_clone_df$age ) )

# ----------------------------------------------------------------------------
# plot batch data as boxplots with trendlines
# ----------------------------------------------------------------------------

ggplot( ) +
  geom_boxplot( data = combined_clone_df, aes(x=num_age, y=counts, group = num_age), color = "grey40" ) +
  geom_smooth( data = combined_clone_df, aes(x=num_age, y=counts, colour = 'linear'), method = "glm", level = 0.99, se= FALSE) +
  geom_smooth( data = combined_clone_df, aes(x=num_age, y=counts, color = "polynomial"), method = 'glm', formula = y ~ poly(x,2), se= FALSE) +
  facet_wrap(~cloneType, scales = "free_y") + theme_classic()

ggsave( filename = paste( base, "/clonal_analysis/output/trendline_box_plot.pdf", sep=""), width = 8, height = 4 )

# ----------------------------------------------------------------------------
# plot batch data as dots with trendlines
# ----------------------------------------------------------------------------

ggplot( ) +
  geom_point( data = combined_clone_df, aes(x=num_age, y=counts, group = num_age), color = "grey40" ) +
  geom_smooth( data = combined_clone_df, aes(x=num_age, y=counts, colour = 'linear'), method = "glm", level = 0.99, se= FALSE) +
  geom_smooth( data = combined_clone_df, aes(x=num_age, y=counts, color = "polynomial"), method = 'glm', formula = y ~ poly(x,2), se= FALSE) +
  facet_wrap(~cloneType, scales = "free_y") + theme_classic()

ggsave( filename = paste( base, "/clonal_analysis/output/trendline_dot_plot.pdf", sep=""), width = 8, height = 4 )

# ----------------------------------------------------------------------------
# plot the data as bar plots with error bars and trend lines - used in figure
# ----------------------------------------------------------------------------
summary_clones$num_age <- as.numeric( gsub( pattern = "D", replacement = "", x = summary_clones$age ) )

ggplot( ) +
  geom_bar( data = summary_clones, aes(x=num_age, y=mean), stat="identity", position="dodge", fill = "grey60") +
  geom_errorbar( data = summary_clones, aes(x=num_age, ymax=mean+se, ymin=mean-se) ) +
  geom_smooth( data = combined_clone_df, aes(x=num_age, y=counts, colour = 'linear'), method = "glm", level = 0.99, se= FALSE) +
  geom_smooth( data = combined_clone_df, aes(x=num_age, y=counts, color = "polynomial"), method = 'glm', formula = y ~ poly(x,2), se= FALSE) +
  facet_wrap(~cloneType, scales = "free_y") + theme_classic()

ggsave( filename = paste( base, "/clonal_analysis/output/Ext_Data_Fig16def_trendline_bar_plot.pdf", sep=""), width = 8, height = 4 )

# --------------------------------
# plotting of batch resolved data
# --------------------------------

list_of_plots <- list()
for (cloneType in unique(combined_clone_df$cloneType)) {
  df2plot <- combined_clone_df[which(combined_clone_df$cloneType == cloneType),]
  list_of_plots[[cloneType]] <- ggplot( df2plot, aes(x=batch, y=counts, color = batch)) + 
    geom_beeswarm(cex=5, size = 1) + ylim(0, (max( df2plot$counts )+1)) +
    facet_wrap(~age, nrow=1, scales = "fixed") + theme_classic() + ggtitle( cloneType ) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

plot_grid( plotlist = list_of_plots, nrow = 3)
ggsave( filename = paste( base, "/clonal_analysis/output/Ext_Data_Fig16abc_batch_plot.pdf", width = 8, height = 6 ) )

# ================================================
# perform ANOVA to test for age and batch effects
# ================================================

# -----------------
# asymmetric clones
# -----------------

model <- aov( counts ~ age * batch, data = combined_clone_df[ which(combined_clone_df$cloneType == "Asym"),])
sum.model <- summary(model)
sum.model <- as.data.frame( sum.model[[1]] )
write.csv( x = sum.model, file = paste( base, "/clonal_analysis/output/Ext_Data_Fig16_ANOVA_Asym_counts_vs_age_and_batch.csv", sep=""), row.names = T )

#              Df Sum Sq Mean Sq F value Pr(>F)
# age           6    453   75.57   1.801  0.105
# batch         3    172   57.38   1.368  0.256
# age:batch    12    281   23.44   0.559  0.871
# Residuals   117   4908   41.95  

# -------------------
# proliferative clones
# -------------------

model <- aov( counts ~ age * batch, data = combined_clone_df[ which(combined_clone_df$cloneType == "Prolif"),])
sum.model <- summary(model)
sum.model <- as.data.frame( sum.model[[1]] )

write.csv( x = sum.model, file = paste( base, "/clonal_analysis/output/Ext_Data_Fig16_ANOVA_Prolif_counts_vs_age_and_batch.csv", sep=""), row.names = T )

#              Df Sum Sq Mean Sq F value  Pr(>F)   
# age           5  13571  2714.1   4.225 0.00128 **
# batch         3   3128  1042.8   1.623 0.18640   
# age:batch    10  19629  1962.9   3.056 0.00149 **
# Residuals   148  95064   642.3                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ----------------------------------------------------------
# since there was a significant interaction of age and batch
# do pairwise comparisons
# ----------------------------------------------------------

all_pw <- TukeyHSD(model,  which = "age:batch" )
all_pw <- as.data.frame( all_pw$`age:batch` )
all_pw[ which(all_pw$`p adj` < 0.05), ]
all_pw <- all_pw[ order( all_pw$`p adj`, decreasing = F), ]

# seems to be mainly D8:Batch 3
#                             diff         lwr        upr      p adj
# D8:Batch 3-D8:Batch 1   42.92935    4.363072  81.495624 0.01242961
# D8:Batch 3-D12:Batch 1  58.33929    9.711176 106.967395 0.00378169
# D8:Batch 3-D13:Batch 1  60.42500    6.860452 113.989548 0.01019426
# D8:Batch 3-D9:Batch 2   50.50000    3.520783  97.479217 0.02032927
# D8:Batch 3-D10:Batch 2  41.72024    2.682906  80.757570 0.02197933
# D8:Batch 3-D13:Batch 2  58.62500    5.060452 112.189548 0.01584917
# D10:Batch 3-D8:Batch 3 -46.31731  -88.538330  -4.096285 0.01533668
# D13:Batch 3-D8:Batch 3 -66.62500 -130.235160  -3.014840 0.02859838

length(which(all_pw$`p adj` < 0.05)) / nrow(all_pw)
# [1] 0.02898551

# write out pairwise comparisons
write.csv( x = all_pw[which(!is.na(all_pw$diff)),], file = paste( base, "/clonal_analysis/output/Ext_Data_Fig16_ANOVA_Prolif_counts_vs_age_and_batch_pw_comp.csv", sep=""), row.names = T )

# --------------------
# analyse small clones
# --------------------

model <- aov( counts ~ age * batch, data = combined_clone_df[ which(combined_clone_df$cloneType == "Small"),])
sum.model <- summary(model)
sum.model <- as.data.frame( sum.model[[1]] )
write.csv( x = sum.model, file = paste( base, "/clonal_analysis/output/Ext_Data_Fig16_ANOVA_Small_counts_vs_age_and_batch.csv", sep=""), row.names = T )

#              Df Sum Sq Mean Sq F value Pr(>F)  
# age           6  16.79  2.7980   2.011 0.0681 .
# batch         3   3.04  1.0142   0.729 0.5364  
# age:batch    12  10.29  0.8579   0.617 0.8254  
# Residuals   139 193.39  1.3913                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# =============================
# fit linear regression model
# =============================

# use a mixed model to include replicates as a random effect
library (glmmTMB)
library (rlang)
library (emmeans)

# save output as xlsx

wb <- createWorkbook()

for (type in unique(combined_clone_df$cloneType)) {
  message( type )
  data_df <- combined_clone_df[ which(combined_clone_df$cloneType == type),]
  data_df$num_batch <- gsub(pattern = "Batch ", replacement = "", x = data_df$batch)
  m_age <- glmmTMB(counts ~ age + (1|num_batch),
                   family = nbinom2, data = data_df)
  
  # test whether “age affects counts”
  m_null <- update(m_age, . ~ . - age)
  aov <- anova(m_null, m_age)
  
  aov <- as.data.frame(aov)
    
  sheetName <- paste( "nbinom2 and LRT", type)
  addWorksheet(wb, sheetName)
  writeData( wb, sheetName, aov )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol(aov) )
  setColWidths( wb, sheetName, cols = 1:ncol(aov), widths="auto" )
  
}

# "Small" clone data gave a warning - 
diagnose( m_age )

# Unusually large coefficients (|x|>10):
# test for alternative models 

m_nb2 <- glmmTMB(counts ~ age + (1|num_batch),
                 family = nbinom2, data = data_df)

m_nb1 <- glmmTMB(counts ~ age + (1|num_batch),
                 family = nbinom1, data = data_df)

m_com <- glmmTMB(counts ~ age + (1|num_batch),
                 family = compois, data = data_df)

m_pois <- glmmTMB(counts ~ age + (1|num_batch),
                 family = poisson, data = data_df)

AIC(m_nb2, m_nb1, m_com, m_pois)

#        df      AIC
# m_nb2   9 564.5097
# m_nb1   9 564.5097
# m_com   9 503.7189
# m_pois  8 562.5097

# m_com looks best, low AIC

# test whether “age affects counts”
m_null <- update(m_com, . ~ . - age)
aov <- anova(m_null, m_com)

aov <- as.data.frame(aov)

sheetName <- paste( "compois and LRT", type)
addWorksheet(wb, sheetName)
writeData( wb, sheetName, aov )
addFilter( wb, sheetName, row = 1, cols = 1:ncol(aov) )
setColWidths( wb, sheetName, cols = 1:ncol(aov), widths="auto" )

saveWorkbook(wb, file =paste( base, "/clonal_analysis/output/Ext_Data_Fig16_mixedmodel_stats.xlsx", sep=""), overwrite = T) 


# ==================================
# analyse lineage restricted clones
# Extended Data Figure 21
# ==================================

# ---------------------------------------------------
# collect data and bring in a format easy to analyse
# ---------------------------------------------------

sheets <- getSheetNames( file = paste( base, "/clonal_analysis/data/Lineage_Restriction_Batch.xlsx", sep="" ) )
list_of_data <- list()
list_of_plots <- list() 
  
for( idx in which(grepl(pattern = "Summary", x = sheets)) ) {
  data <- read.xlsx( xlsxFile = paste( base, "/clonal_analysis/data/Lineage_Restriction_Batch.xlsx", sep="" ), sheet = sheets[idx])
  data <- data[, which( grepl( pattern = "Number", x = data[1,] ) )]
  data <- apply(data, 2, as.numeric)
  data <- as.data.frame( data[which(!is.na(data[,1])),] )
  data$group <- rep( c("total", "SATB2", "CTIP2"), nrow(data)/3 )
  
  data$rep <- rep(1:(nrow(data)/3), each=3)
  data_pure <- data[which(!data$group == "total"),]
  total_df <- data[which(data$group == "total"),]
  colnames(data) <- c("Prolif", "Asym", "SN", "group", "rep")
  colnames(data_pure) <- c("Prolif", "Asym", "SN", "group", "rep")
  colnames(total_df) <- c("Prolif", "Asym", "SN", "group", "rep")
  
  df2plot <- reshape2::melt( data_pure, id=c("group", "rep") )
  total_df <- reshape2::melt( total_df, id=c("group", "rep") )
  total_df$rep_group <- paste(total_df$variable, total_df$rep, sep="_" )
  df2plot$rep_group <- paste(df2plot$variable, df2plot$rep, sep="_" )
  
  df2plot <- merge(df2plot, total_df[,c("rep_group", "value")], by="rep_group")
  colnames( df2plot ) <- c("rep_group", "group", "rep", "variable", "count", "total")
  df2plot$rel <- df2plot$count / df2plot$total
  
  list_of_plots[[ length(list_of_plots)+1 ]] <- ggplot( df2plot, aes(x=rep, y=rel)) + geom_bar(stat="identity", position="dodge") +
    facet_wrap(variable ~ group, ncol = 2) + theme_classic() + ggtitle( sheets[idx] )
  
  list_of_data[[ length(list_of_data)+1 ]] <- df2plot
}

# ----------------
# plot for figure 
# ----------------

plot_grid( plotlist = list_of_plots )
ggsave( filename = paste( base, "/clonal_analysis/output/Ext_Data_Fig21_batch_plot_restriction.pdf", width = 4, height = 6 ))

av_out <- list()
raw_out <- list()
for (i in 1:2) {
  tmp_df <- list_of_data[[ i ]]
  colnames(tmp_df) <- c("rep_group", "layer", "rep", "clone_type", "count", "total", "rel")
  tmp_df$asin <- asin(sqrt( as.numeric( tmp_df$rel ) ) )
  tmp_df$clone_type <- factor(tmp_df$clone_type)
  test <- aov(data = tmp_df, formula = asin ~ layer * clone_type )
  av_out[[i]] <- as.data.frame( summary( test )[[1]] )
  av_out[[i]]$comp <- rownames( av_out[[i]] )
  raw_out[[i]] <- tmp_df
}

# -------------------------------------------------
# write out data for reporting in statistics table 
# -------------------------------------------------

wb <- createWorkbook()
for (i in 1:2) {
  sheetName <- paste( "RAW", i)
  addWorksheet(wb, sheetName)
  writeData( wb, sheetName, raw_out[[i]] )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol(raw_out[[i]]) )
  setColWidths( wb, sheetName, cols = 1:ncol(raw_out[[i]]), widths="auto" )
  
  sheetName <- paste( "ANOVA", i)
  addWorksheet(wb, sheetName)
  writeData( wb, sheetName, av_out[[i]] )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol(av_out[[i]]) )
  setColWidths( wb, sheetName, cols = 1:ncol(av_out[[i]]), widths="auto" )
}
saveWorkbook(wb, file =paste( base, "/clonal_analysis/output/Ext_Data_Fig_21_Two_way_ANOVA_lineage_restricted.xlsx", sep=""), overwrite = T) 

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
# [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=de_AT.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] emmeans_2.0.3    rlang_1.1.3      glmmTMB_1.1.14   cowplot_1.1.3    dplyr_1.1.4      ggbeeswarm_0.7.2
# [7] ggplot2_3.5.2    openxlsx_4.2.5.2 tidyxl_1.0.10   
# 
# loaded via a namespace (and not attached):
#  [1] sandwich_3.1-1      utf8_1.2.4          generics_0.1.3      stringi_1.8.3       lattice_0.21-9     
#  [6] lme4_1.1-35.3       magrittr_2.0.3      estimability_1.5.1  grid_4.3.2          RColorBrewer_1.1-3 
# [11] mvtnorm_1.4-1       plyr_1.8.9          Matrix_1.6-5        zip_2.3.1           mgcv_1.9-0         
# [16] fansi_1.0.6         scales_1.4.0        numDeriv_2016.8-1.1 reformulas_0.4.4    Rdpack_2.6.6       
# [21] cli_3.6.2           rbibutils_2.4.1     splines_4.3.2       withr_3.0.0         packrat_0.9.2      
# [26] tools_4.3.2         reshape2_1.4.4      coda_0.19-4.1       nloptr_2.0.3        minqa_1.2.6        
# [31] boot_1.3-28.1       vctrs_0.6.5         R6_2.5.1            zoo_1.8-12          lifecycle_1.0.4    
# [36] stringr_1.5.1       vipor_0.4.7         MASS_7.3-60         pkgconfig_2.0.3     beeswarm_0.4.0     
# [41] pillar_1.9.0        gtable_0.3.6        glue_1.7.0          Rcpp_1.0.12         tibble_3.2.1       
# [46] tidyselect_1.2.1    rstudioapi_0.16.0   xtable_1.8-4        farver_2.1.1        nlme_3.1-163       
# [51] TMB_1.9.21          compiler_4.3.2     
