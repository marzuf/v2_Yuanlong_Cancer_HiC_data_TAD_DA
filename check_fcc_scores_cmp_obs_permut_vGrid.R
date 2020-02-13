
options(scipen=100)

# Rscript check_fcc_scores_cmp_obs_permut_vGrid.R

script_name <- "check_fcc_scores_cmp_obs_permut_vGrid.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")


require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)


ratio_vect <- seq(from=0, to=1, by=0.1)
ratio_fract_names <- paste0("ratio > ", ratio_vect[1:(length(ratio_vect)-1)], " and ratio <= ",ratio_vect[2:length(ratio_vect)])
ratio_fract_names <- paste0("ratio \u2208 ]", ratio_vect[1:(length(ratio_vect)-1)], ", ",ratio_vect[2:length(ratio_vect)], "]")
ratio_fract_names <- paste0("]", ratio_vect[1:(length(ratio_vect)-1)], ", ",ratio_vect[2:length(ratio_vect)], "]")
ratio_fract_names[ratio_fract_names == "]0, 0.1]"] <- "[0, 0.1]"


plotType <- "svg"
myHeight <- 8
myWidth <- 12


outFolder <- file.path("CHECK_FCC_SCORES_CMP_OBS_PERMUT_VGRID")
dir.create(outFolder, recursive = TRUE)

obs_mat_dt <- get(load("CHECK_FCC_SCORES_VGRID/mat_dt.Rdata"))
obs_mat_dt[1:5,1:5]

all_sfx <- c("all", "left", "right")

for(side in all_sfx) {
  
  
  permut_mat_dt <- get(load(file.path("CHECK_FCC_SCORES_PERMUT_VGRID", paste0(side, "permut_mat_dt.Rdata"))))
  permut_mat_dt[1:5,1:5]
  
  stopifnot(colnames(obs_mat_dt) == colnames(permut_mat_dt))
  stopifnot(rownames(obs_mat_dt) == rownames(permut_mat_dt))
  
  obs_mat_dt_norm <- obs_mat_dt/max(obs_mat_dt)
  permut_mat_dt_norm <- permut_mat_dt/max(permut_mat_dt)
  
  ratio_obs_permut_raw <- permut_mat_dt/obs_mat_dt
  ratio_obs_permut_maxNorm <- permut_mat_dt_norm/obs_mat_dt_norm
  
  
  cmp=""
  for(cmp in c("_raw", "_maxNorm")) {
    
    plot_dt <- melt(get(paste0("ratio_obs_permut", cmp)))
    
    stopifnot(grepl("ratioDown", plot_dt$Var1))
    stopifnot(grepl("ratioFC", plot_dt$Var2))
    
    plot_dt$Var1 <- factor(plot_dt$Var1, levels = paste0("ratioDown_", unique(sort(as.numeric(gsub("ratioDown_", "", plot_dt$Var1))))))
    plot_dt$Var2 <- factor(plot_dt$Var2, levels = paste0("ratioFC_", unique(sort(as.numeric(gsub("ratioFC_", "", plot_dt$Var2))))))
    
    
    plot_dt$Var1_lab <- ratio_fract_names[as.numeric(gsub("ratioDown_", "", plot_dt$Var1))]
    plot_dt$Var1_lab <- factor(plot_dt$Var1_lab, levels = ratio_fract_names)
    
    plot_dt$Var2_lab <- ratio_fract_names[as.numeric(gsub("ratioFC_", "", plot_dt$Var2))]
    plot_dt$Var2_lab <- factor(plot_dt$Var2_lab, levels = ratio_fract_names)
    
    stopifnot(!is.na(plot_dt$Var1))
    stopifnot(!is.na(plot_dt$Var1_lab))
    stopifnot(!is.na(plot_dt$Var2))
    stopifnot(!is.na(plot_dt$Var2_lab))
    
    grid_plot <- ggplot(data = plot_dt, aes(x=Var1_lab, y=Var2_lab, fill = value))+
      
      
      ggtitle(paste0("ratio # TADs permut/obs"),   
              subtitle = paste0(side, " permut. data - ", gsub("_", "", cmp), " ratio"))+
      
      scale_x_discrete(name="ratioDown")  + 
      scale_y_discrete(name="ratioFC")  + 
      geom_tile() +
      labs(fill = "ratio # TADs permut/obs")+
      scale_fill_gradient( low="blue", high="red", na.value = "white" )  +
      # scale_fill_gradientn(colours = terrain.colors(10))+
      theme(
        axis.text.x = element_text(colour = "black", size=12),
        axis.text.y= element_text(colour = "black", size=12),
        axis.title.x = element_text(colour = "black", size=14, face="bold"),
        axis.title.y = element_text(colour = "black", size=14, face="bold"),
        plot.title = element_text(hjust=0.5, size=16, face="bold"),
        plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
        panel.background = element_rect(fill = "transparent")
        # legend.background =  element_rect()
      )
    
    outFile <- file.path(outFolder, paste0("ratioDown_ratioFC_", side, "obs_permut", cmp, "Ratio_gridPlot.", plotType))
    ggsave(grid_plot, file = outFile, height=myHeight, width=myWidth)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





