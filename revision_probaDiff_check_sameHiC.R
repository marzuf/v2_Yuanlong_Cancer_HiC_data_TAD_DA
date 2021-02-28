


# Rscript revision_probaDiff_check_sameHiC.R

require(ggpubr)
require(ggsci)
require(doMC)
require(foreach)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

#outFolder <- file.path("REVISION_PROBADIFF_CHECK_SAMEHIC_V2")
#dir.create(outFolder, recursive = TRUE)
#all_inter_intra_dt <- get(load("REVISION_INTER_INTRA_PROBA_V2_CORRECTED_SAMEHIC///all_inter_intra_dt.Rdata"))

outFolder <- file.path("REVISION_PROBADIFF_CHECK_SAMEHIC")
dir.create(outFolder, recursive = TRUE)
all_inter_intra_dt <- get(load("REVISION_INTER_INTRA_PROBA_CORRECTED_SAMEHIC///all_inter_intra_dt.Rdata"))

plotType <- "png"
myWidth <- 500
myWidthGG <- 7
myHeightGG <- 5
myHeight <- 400

plotCex <- 1.2

tadSignifThresh <- 0.01
plot_qt1 <- 0.05
plot_qt2 <- 0.95

mytheme <-     theme(
  # text = element_text(family=fontFamily),
  panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.background = element_rect(fill = "transparent"),
  panel.grid.major.x =  element_blank(),
  panel.grid.minor.x =  element_blank(),
  axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
  axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
  plot.title = element_text(hjust=0.5, size = 16, face="bold"),
  plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
  legend.title = element_text(face="bold")
) 

all_pairs <- c(
  file.path("LI_40kb","GSE105381_HepG2_40kb", "TCGAlihc_norm_lihc"),
  file.path("LG1_40kb" ,"ENCSR444WCZ_A549_40kb", "TCGAluad_norm_luad"),
  file.path("LG2_40kb" ,"ENCSR444WCZ_A549_40kb" ,"TCGAluad_norm_luad"),
  file.path("LG1_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("LG2_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "TCGAprad_norm_prad"),
  file.path("GSE118514_RWPE1_40kb", "GSE118514_22Rv1_40kb", "TCGAprad_norm_prad")
)



all_inter_intra_dt$mean_inter_prevnext <- 0.5*(all_inter_intra_dt$mean_inter_prev + all_inter_intra_dt$mean_inter_next)
all_inter_intra_dt$mean_inter_prevnextNorm <- 0.5*(all_inter_intra_dt$mean_inter_prevNorm + 
                                                     all_inter_intra_dt$mean_inter_nextNorm)

all_inter_intra_dt$matched_mean_inter_prevnext <- 0.5*(all_inter_intra_dt$matched_mean_inter_prev + 
                                                         all_inter_intra_dt$matched_mean_inter_next)
all_inter_intra_dt$matched_mean_inter_prevnextNorm <- 0.5*(all_inter_intra_dt$matched_mean_inter_prevNorm + 
                                                             all_inter_intra_dt$matched_mean_inter_nextNorm)




# first have a look of count matching
all_vars <- c("mean_intra", "mean_intraNorm","mean_inter_prevnext" , "mean_inter_prevnextNorm")

plot_var = "mean_intra"
for(plot_var in all_vars) {
  
  
  plotTit <- paste0(plot_var, " ds vs. matched ds")
  
  mySub <- paste0("# DS = ", length(unique(all_inter_intra_dt$hicds)), 
                  "; # cmps = ", length(unique(file.path(all_inter_intra_dt$matched_hicds, all_inter_intra_dt$hicds))), 
                  "; # TADs = ", nrow(all_inter_intra_dt))
  
  
  myx <- all_inter_intra_dt[,paste0(plot_var)]
  myy <- all_inter_intra_dt[,paste0("matched_", plot_var)]
  
  
  
  outFile  <- file.path(outFolder, paste0(plot_var, "_ds_vs_matched_densplot.", plotType))
  do.call(plotType, list(outFile, height=myWidth, width=myWidth))
  densplot(x = myx,
           xlab =paste0(plot_var, " - ds"),
           y  = myy,
           ylab = paste0(plot_var, " - matched"),
           main=plotTit,
           cex.main=plotCex,
           cex.axis =plotCex,
           cex.lab = plotCex)
  mtext(side=3, text=mySub)
  addCorr(x=myx,y=myy, legPos="topleft", bty="n", corMet="spearman")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  plotTit <- paste0(plot_var, " ds vs. matched ds (log10)")
  
  myx <- log10(all_inter_intra_dt[,paste0(plot_var)])
  myy <- log10(all_inter_intra_dt[,paste0("matched_", plot_var)])
  
  outFile  <- file.path(outFolder, paste0(plot_var, "_ds_vs_matched_log10_densplot.", plotType))
  do.call(plotType, list(outFile, height=myWidth, width=myWidth))
  densplot(x = myx,
           xlab =paste0(plot_var, " - ds [log10]"),
           y  = myy,
           ylab = paste0(plot_var, " - matched [log10]"),
           main=plotTit,
           cex.main=plotCex,
           cex.axis =plotCex,
           cex.lab = plotCex)
  mtext(side=3, text=mySub)
  addCorr(x=myx,y=myy, legPos="topleft", bty="n", corMet="spearman")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}


