





#########################
######################### desc family clusters
#########################

all_fams_agg_cluster_dt$clusterLength <- all_fams_agg_cluster_dt$maxPos-all_fams_agg_cluster_dt$midPos+1
all_fams_agg_cluster_dt$clusterLength_log10 <- log10(all_fams_agg_cluster_dt$clusterLength)
all_fams_agg_cluster_dt$clusterSize_log10 <- log10(all_fams_agg_cluster_dt$clusterSize)

fams_big1_dt <- all_fams_agg_cluster_dt[all_fams_agg_cluster_dt$clusterSize > 1,]

plot_vars <- c("clusterSize", "clusterSize_log10", "clusterLength", "clusterLength_log10")

for(plot_var in plot_vars){
  
  plotTit <- "Family cluster size"
  subTit <- paste0("# family clusters = ", length(unique(file.path(all_fams_agg_cluster_dt$clusterID, all_fams_agg_cluster_dt$family))),
                   " (# fam. = ", length(unique(all_fams_agg_cluster_dt$family)), ")") 
  
  
  p <- ggdensity(all_fams_agg_cluster_dt,
                 x = paste0(plot_var),
                 y = "..density..",
                 # combine = TRUE,                  # Combine the 3 plots
                 xlab = paste0(" family ", plot_var),
                 # add = "median",                  # Add median line.
                 rug = FALSE,                      # Add marginal rug
                 # color = "dataType",
                 # fill = "dataType",
                 palette = "d3"
  )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
    theme(plot.title=element_text(face="bold"),
          plot.subtitle=element_text(face="italic")
    )
  outFile <- file.path(outFolder, paste0("family_", plot_var, "_densityplot.", plotTypeGG))
  ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  plotTit <- "Family cluster size"
  subTit <- paste0("# family clusters = ", length(unique(file.path(fams_big1_dt$clusterID, fams_big1_dt$family))),
                   " (# fam. = ", length(unique(fams_big1_dt$family)), ")") 
  
  p <- ggdensity(fams_big1_dt,
                 x = paste0(plot_var),
                 y = "..density..",
                 # combine = TRUE,                  # Combine the 3 plots
                 xlab = paste0(" family ", plot_var),
                 # add = "median",                  # Add median line.
                 rug = FALSE,                      # Add marginal rug
                 # color = "dataType",
                 # fill = "dataType",
                 palette = "d3"
  )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
    theme(plot.title=element_text(face="bold"),
          plot.subtitle=element_text(face="italic")
    )
  outFile <- file.path(outFolder, paste0("family_", plot_var, "_bigger1_densityplot.", plotTypeGG))
  ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
}




