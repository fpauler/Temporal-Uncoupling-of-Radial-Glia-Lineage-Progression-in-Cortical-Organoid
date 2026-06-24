##############################
###
#    statistics for for Fig 3I
###
##############################

# fill data frame with values according to Supplemental Table 2
data <- data.frame( label = c("Organoid D8", "Organoid D13", "*In vivo"), 
            n_restr = c(31, 25, 11), 
            total = c(87, 78, 386) )

# calculate the Chi-square statistics

pvalues <- list()
for (row in c(1,2)) {
  
  tmp_table <- as.table( rbind( c(data[row, "n_restr"], data[row, "total"] - data[row, "n_restr"]),
                                c(data[3, "n_restr"], data[3, "total"] - data[3, "n_restr"]) )
  )
  chisq_out <- chisq.test( tmp_table, correct = T, rescale.p = F, simulate.p.value = F )
  pvalues[[row]] <- chisq_out
  
}

# pvalues
# [[1]]
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  tmp_table
# X-squared = 90.294, df = 1, p-value < 2.2e-16
# 
# 
# [[2]]
# 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  tmp_table
# X-squared = 73.289, df = 1, p-value < 2.2e-16

pvalues[[1]]$p.value
# [1] 2.052281e-21

pvalues[[2]]$p.value
# [1] 1.120043e-17