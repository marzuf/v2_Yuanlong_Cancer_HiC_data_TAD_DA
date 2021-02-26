
# need to re run 26.02.21 on corrected proba data

# outFolder <- file.path("REVISION_PROBADIFF")
# dir.create(outFolder, recursive = TRUE)
# all_inter_intra1_dt <- get(load("REVISION_INTER_INTRA_PROBA/all_inter_intra_dt.Rdata"))
# all_inter_intra2_dt <- get(load("REVISION_INTER_INTRA_PROBA2/all_inter_intra_dt.Rdata"))

outFolder <- file.path("REVISION_PROBADIFF_V2_CORRECTED")
dir.create(outFolder, recursive = TRUE)
all_inter_intra1_dt <- get(load("REVISION_INTER_INTRA_PROBA_V2_CORRECTED/all_inter_intra_dt.Rdata"))
all_inter_intra2_dt <- get(load("REVISION_INTER_INTRA_PROBA2_V2_CORRECTED/all_inter_intra_dt.Rdata"))


# Rscript revision_probaDiff.R

require(ggpubr)
require(ggsci)
require(doMC)
require(foreach)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myWidth <- 500
myWidthGG <- 7
myHeightGG <- 5
myHeight <- 400

plotCex <- 1.2

tadSignifThresh <- 0.01
plot_qt1 <- 0.05
plot_qt2 <- 0.95

# 1st part => can be done for all data
# pval vs. inter/intra 
final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
final_table_DT$signif_lab <- ifelse(final_table_DT$adjPvalComb <= 0.01, "signif.", "not signif.")

stopifnot(! all_inter_intra1_dt$hicds %in% all_inter_intra2_dt$hicds)
all_inter_intra_dt <- rbind(all_inter_intra1_dt, all_inter_intra2_dt)
stopifnot(final_table_DT$hicds %in% all_inter_intra_dt$hicds)

# compute inter/intra for the next
all_inter_intra_dt$interOverIntraNorm_next <- all_inter_intra_dt$mean_inter_nextNorm/all_inter_intra_dt$mean_intraNorm
# compute inter/intra for the prev
all_inter_intra_dt$interOverIntraNorm_prev <- all_inter_intra_dt$mean_inter_prevNorm/all_inter_intra_dt$mean_intraNorm
# compute the average next and prev inter/intra
all_inter_intra_dt$interOverIntraNorm_mean <- sapply(1:nrow(all_inter_intra_dt), function(x) {
  # if there was no next TAD -> keep only prev ratio
  if(is.na(all_inter_intra_dt$mean_inter_nextNorm[x])) {
    return(all_inter_intra_dt$interOverIntraNorm_prev[x])
  }
  # if there was no prev TAD -> keep only next ratio
  if(is.na(all_inter_intra_dt$mean_inter_prevNorm[x])) {
    return(all_inter_intra_dt$interOverIntraNorm_next[x])
  }
  return(0.5*(all_inter_intra_dt$interOverIntraNorm_next[x]+all_inter_intra_dt$interOverIntraNorm_prev[x]))
})

###### for the not norm
# compute inter/intra for the next
all_inter_intra_dt$interOverIntra_next <- all_inter_intra_dt$mean_inter_next/all_inter_intra_dt$mean_intra
# compute inter/intra for the prev
all_inter_intra_dt$interOverIntra_prev <- all_inter_intra_dt$mean_inter_prev/all_inter_intra_dt$mean_intra
# compute the average next and prev inter/intra
all_inter_intra_dt$interOverIntra_mean <- sapply(1:nrow(all_inter_intra_dt), function(x) {
  # if there was no next TAD -> keep only prev ratio
  if(is.na(all_inter_intra_dt$mean_inter_next[x])) {
    return(all_inter_intra_dt$interOverIntra_prev[x])
  }
  # if there was no prev TAD -> keep only next ratio
  if(is.na(all_inter_intra_dt$mean_inter_prev[x])) {
    return(all_inter_intra_dt$interOverIntra_next[x])
  }
  return(0.5*(all_inter_intra_dt$interOverIntra_next[x]+all_inter_intra_dt$interOverIntra_prev[x]))
})




all_inter_intra_dt$hicds_regionID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)
final_table_DT$hicds_regionID <- file.path(final_table_DT$hicds, final_table_DT$region)
stopifnot(final_table_DT$hicds_regionID %in% all_inter_intra_dt$hicds_regionID)

merge_dt <- merge(final_table_DT, all_inter_intra_dt, all=FALSE, by=c("hicds", "region", "hicds_regionID"))

legTitle <- ""
plot_var = "interOverIntraNorm_mean"


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


all_plot_vars <- c("interOverIntraNorm_mean", "interOverIntra_mean", "mean_intra")

#### COMPARE RATIO INTER/INTRA + mean_intra + IN SIGNIF VS. NOT SIGNIF
# 1) density
# 2 continuous association with pvalue)
plot_var="mean_intra"

foo <- foreach(plot_var = all_plot_vars) %dopar%{
  
  
  qts <- quantile(na.omit(merge_dt[,paste0(plot_var)]), probs=c(plot_qt1, plot_qt2))
  
  
  sub_dt <- merge_dt[merge_dt[,paste0(plot_var)] >= qts[1] & merge_dt[,paste0(plot_var)] <= qts[2] ,]
  
  sub_dt <- sub_dt[!is.na(sub_dt[,paste0(plot_var)]),]
  
  all_cmps <- unique(file.path(sub_dt$hicds, sub_dt$exprds))
  
  plotTit <- paste0(plot_var, " - signif. not signif.")
  mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(sub_dt), 
                  " (signif.: ", sum(sub_dt$adjPvalComb <= tadSignifThresh), "); crop ", plot_qt1,"-", plot_qt2, "")
  
  my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(sub_dt$signif_lab))
  
  p3 <- ggdensity(sub_dt,
                  x = paste0(plot_var),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0(plot_var),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "signif_lab",
                  fill = "signif_lab",
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = mySub)+
    scale_color_manual(values=my_cols)+
    scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    mytheme
  
  outFile <- file.path(outFolder, paste0(plot_var, "_signif_notsignif_density.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  plotTit <- paste0(plot_var, " vs. adjPvalComb")
  
  outFile  <- file.path(outFolder, paste0(plot_var, "_log10_vs_adjPvalComb_densplot.", plotType))
  do.call(plotType, list(outFile, height=myWidth, width=myWidth))
  densplot(x = -log10(sub_dt$adjPvalComb),
           xlab ="adjPvalComb [-log10]",
       y  = log10(sub_dt[,paste0(plot_var)]),
       ylab = paste0(plot_var, " [log10]"),
       main=plotTit,
       cex.main=plotCex,
       cex.axis =plotCex,
       cex.lab = plotCex)
  mtext(side=3, text=mySub)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile  <- file.path(outFolder, paste0(plot_var, "_vs_adjPvalComb_densplot.", plotType))
  do.call(plotType, list(outFile, height=myWidth, width=myWidth))
  densplot(x = -log10(sub_dt$adjPvalComb),
           xlab ="adjPvalComb [-log10]",
           y  = sub_dt[,paste0(plot_var)],
           ylab = paste0(plot_var),
           main=plotTit,
           cex.main=plotCex,
           cex.axis =plotCex,
           cex.lab = plotCex)
  mtext(side=3, text=mySub)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}



  

 #################################################################



# qts <- quantile(na.omit(merge_dt$mean_intra), probs=c(0.05, 0.95))
# sub_dt <- merge_dt[merge_dt$mean_intra >= qts[1] & merge_dt$mean_intra <= qts[2] ,]
# 
# ggdensity(sub_dt,
#           x = paste0("mean_intra"),
#           y = "..density..",
#           # combine = TRUE,                  # Combine the 3 plots
#           xlab = paste0("mean_intra"),
#           # add = "median",                  # Add median line.
#           rug = FALSE,                      # Add marginal rug
#           color = "signif_lab",
#           fill = "signif_lab",
#           palette = "jco"
# ) 




# qts <- quantile(na.omit(merge_dt$interOverIntra_mean), probs=c(0.05, 0.95))
# 
# sub_dt <- merge_dt[merge_dt$interOverIntra_mean >= qts[1] & merge_dt$interOverIntra_mean <= qts[2] ,]
# ggdensity(sub_dt,
#                 x = paste0("interOverIntra_mean"),
#                 y = "..density..",
#                 # combine = TRUE,                  # Combine the 3 plots
#                 xlab = paste0("interOverIntra_mean"),
#                 # add = "median",                  # Add median line.
#                 rug = FALSE,                      # Add marginal rug
#                 color = "signif_lab",
#                 fill = "signif_lab",
#                 palette = "jco"
# ) 
# 
# 
# 
# 
# 
# # compute inter-intra for the next
# all_inter_intra_dt$interMinusIntraNorm_next <- all_inter_intra_dt$mean_inter_nextNorm-all_inter_intra_dt$mean_intraNorm
# # compute inter/intra for the prev
# all_inter_intra_dt$interMinusIntraNorm_prev <- all_inter_intra_dt$mean_inter_prevNorm-all_inter_intra_dt$mean_intraNorm
# # compute the average next and prev inter/intra
# all_inter_intra_dt$interMinusIntraNorm_mean <- sapply(1:nrow(all_inter_intra_dt), function(x) {
#   # if there was no next TAD -> keep only prev ratio
#   if(is.na(all_inter_intra_dt$mean_inter_nextNorm[x])) {
#     return(all_inter_intra_dt$interMinusIntraNorm_prev[x])
#   }
#   # if there was no prev TAD -> keep only next ratio
#   if(is.na(all_inter_intra_dt$mean_inter_prevNorm[x])) {
#     return(all_inter_intra_dt$interMinusIntraNorm_next[x])
#   }
#   return(0.5*(all_inter_intra_dt$interMinusIntraNorm_next[x]+all_inter_intra_dt$interMinusIntraNorm_prev[x]))
# })
# 
# 
# all_inter_intra_dt$hicds_regionID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)
# final_table_DT$hicds_regionID <- file.path(final_table_DT$hicds, final_table_DT$region)
# stopifnot(final_table_DT$hicds_regionID %in% all_inter_intra_dt$hicds_regionID)
# 
# merge_dt <- merge(final_table_DT, all_inter_intra_dt, all=FALSE, by=c("hicds", "region", "hicds_regionID"))
# 
# qts <- quantile(na.omit(merge_dt$interMinusIntraNorm_mean), probs=c(0.05, 0.95))
# 
# 
# sub_dt <- merge_dt[merge_dt$interMinusIntraNorm_mean >= qts[1] & merge_dt$interMinusIntraNorm_mean <= qts[2] ,]
# 
# ggdensity(sub_dt,
#                 x = paste0("interMinusIntraNorm_mean"),
#                 y = "..density..",
#                 # combine = TRUE,                  # Combine the 3 plots
#                 xlab = paste0("interMinusIntraNorm_mean"),
#                 # add = "median",                  # Add median line.
#                 rug = FALSE,                      # Add marginal rug
#                 color = "signif_lab",
#                 fill = "signif_lab",
#                 palette = "jco"
# ) 
