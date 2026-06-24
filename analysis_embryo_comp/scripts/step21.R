# ===============================================================
# step 21
# analysis clone coverage for different clone restriction types
# ===============================================================

library (openxlsx)
library (dplyr)
library (ggplot2)
library (ggbeeswarm)

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

clone_metaData <- read.xlsx( xlsxFile = paste(base, "/plots/Sup_Fig_22/MADMCloneSeq_metaData.xlsx", sep="") )
clone_metaData <- clone_metaData[ which( !clone_metaData$nn_clusters_dB == "undefined" ), ]

table( clone_metaData$nn_clusters_dB )
#  DL  UL 
# 128  79 

clone_composition <- clone_metaData %>%
  group_by(clone, clone_type, clone_size_patched, clone_size_complete) %>%
  summarise(types = paste( sort(unique(nn_clusters_dB)), collapse=","), n=n() )

clone_composition <- clone_composition[ which( clone_composition$n > 1), ]

table(clone_composition$clone_type)
#  N SN 
# 21 34 

clone_metaData <- clone_metaData[ which( clone_metaData$clone %in% clone_composition$clone ), ]
table( clone_metaData$nn_clusters_dB )
#  DL  UL 
# 119  76 

table( clone_metaData$clone_type )
#   N  SN 
# 109  86 

# analyse clone coverage

clone_composition$coverage <- clone_composition$n / as.numeric( clone_composition$clone_size_complete )

clone_composition$simple_type <- "ML"
clone_composition$simple_type[ which(clone_composition$types == "DL")] <- "DL"
clone_composition$simple_type[ which(clone_composition$types == "UL")] <- "UL"
clone_composition$simple_type_restriction <- paste(clone_composition$clone_type, clone_composition$simple_type, sep="_")

ggplot( clone_composition, aes(x=simple_type_restriction, y=coverage)) + 
  geom_boxplot() + geom_beeswarm() + 
  theme_classic()+ ylim(0,max(clone_composition$coverage))
ggsave( filename = paste(base, "/plots/Sup_Fig_23/S23e_Clone_cov_lineageRestr.pdf", sep="") )

# report stats for paper
clone_composition %>%
  group_by( simple_type_restriction, "clone_type" ) %>%
  summarize( mean = mean(coverage),
            median = median(coverage),
            n = n()) -> Panel_e_stats

Panel_e_stats <- as.data.frame( Panel_e_stats )
write.csv (x = Panel_e_stats, file = paste(base, "plots/Sup_Fig_23/Panel_e_stats.csv", sep=""), row.names = T)

# test for significance of coverage differences
clone_composition$coverage_asin <- asin(sqrt(clone_composition$coverage))
test_aov <- aov(clone_composition$coverage_asin ~ factor(clone_composition$simple_type) * factor(clone_composition$clone_type))
out <- summary (test_aov)
out <- as.data.frame( out[[1]] )
write.csv (x = out, file = paste(base, "plots/Sup_Fig_23/Panel_e_ANOVA_stats.csv", sep=""), row.names = T)


