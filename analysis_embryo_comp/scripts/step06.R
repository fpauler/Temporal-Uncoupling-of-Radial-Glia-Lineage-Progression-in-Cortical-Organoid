# ==========================================================================
# step 06
# plot the abundances GFP positive cells in LaManno reference cell types 
# this produces a more concise plot necessary to display the 4 replicates
# the meta data was saved at step 05
# ==========================================================================

library (dplyr)
library (ggplot2)

set.seed( 2401 )

# define working folders
base_folder <- "~/" # personal base folder
base <- paste( base_folder, "analysis_embryo_comp/", sep="" )
nobackup_base <- paste( base_folder, "noSave/", sep="" ) # all large files not backed up

for (age in c( "D8", "D13", "D20", "D25" ) ) {
  
  meta.data <- readRDS( file = paste( nobackup_base,  "/RDS_files/MetaData.GfpDetect.", age, ".RDS", sep="" ) )
  
  # combined annotation
  tmp_counts <- table(meta.data$SubclassID_class)
  tmp_counts <- tmp_counts / sum(tmp_counts)
  group_vec <- names( which(tmp_counts > 0.03) )
  
  gfp_frac_df <- data.frame()
  for (group in group_vec ) {
    for (replicate in unique(meta.data$orig.ident)) {
      tmp_df <- as.data.frame( table( meta.data$GFP_cluster[ which (meta.data$SubclassID_class == group & meta.data$orig.ident == replicate) ] ) )
      
      if (nrow(tmp_df) > 0) {
        tmp_df$group <- group
        tmp_df$replicate <- replicate
        gfp_frac_df <- rbind( gfp_frac_df, tmp_df)
      }
    }
  }
  
  # determine the total number of cells for each cell type/replicate
  sum_freq <- gfp_frac_df %>%
    group_by(group, replicate) %>%
    summarize(
      sum = sum(Freq)
    )
  sum_freq$group_rep <- paste(sum_freq$group, sum_freq$replicate, sep="_") 
  
  # calculate the relative abundance
  total_vec <- sum_freq$sum
  names( total_vec ) <- sum_freq$group_rep
  
  gfp_frac_df$group_rep <- paste(gfp_frac_df$group, gfp_frac_df$replicate, sep="_") 
  gfp_frac_df$total <- total_vec[ gfp_frac_df$group_rep ]
  gfp_frac_df$rel <- gfp_frac_df$Freq / gfp_frac_df$total
   
  gfp_frac_plot <- ggplot( gfp_frac_df, aes(x=Var1, y=rel, group = Var1, color = Var1)) + 
    geom_boxplot( color = "grey50") + geom_beeswarm(cex=5) + theme_classic() +
    scale_color_manual( values = c("YES" = "cornflowerblue", "NO" = "red")) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    facet_wrap(~group, nrow=1) + ylim(0,1)
  ggsave ( plot = gfp_frac_plot, filename = paste(base, "/plots/Sup_Fig_2_3/organoid_", age, "_GFPfrac.SubclassClass.pdf", sep=""), width = 8, height = 5 )
}
