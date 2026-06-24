##############################
###
#    statistics for for Ext Fig 11RS
###
##############################

# fill data frame with values according to Supplemental Table 2
data <- list()
data[["D13"]] <- data.frame( label = c("Center", "Periphery"), 
            MADM = c(312, 2264), 
            MADM_Casp = c(81, 56) )
data[["D16"]] <- data.frame( label = c("Center", "Periphery"), 
                             MADM = c(217, 2865), 
                             MADM_Casp = c(77, 149) )

# calculate the Chi-square statistics

pvalues <- list()
for (day in c("D13","D16")) {
  
  tmp_table <- as.table( rbind( c(data[[day]][1, 2], data[[day]][1, 3]),
                                c(data[[day]][2, 2], data[[day]][2, 3]) )
  )
  chisq_out <- chisq.test( tmp_table, correct = T, rescale.p = F, simulate.p.value = F )
  pvalues[[day]] <- chisq_out
  
}

pvalues
# $D13
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  tmp_table
# X-squared = 228.31, df = 1, p-value < 2.2e-16
# 
# 
# $D16
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  tmp_table
# X-squared = 186.66, df = 1, p-value < 2.2e-16

pvalues[[1]]$p.value
# [1] 1.391872e-51
pvalues[[2]]$p.value
# [1] 1.707597e-42
> 