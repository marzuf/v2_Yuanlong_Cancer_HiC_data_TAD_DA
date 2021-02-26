

outFolder <- file.path("REVISION_PROBADIFF_PAIREDHIC")
dir.create(outFolder, recursive = TRUE)

# Rscript revision_probaDiff_pairedHiC.R

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




all_inter_intra_dt <- get(load("REVISION_INTER_INTRA_PROBA_V2_CORRECTED_PAIREDHIC//all_inter_intra_dt.Rdata"))

# compute inter over intra - reference
all_inter_intra_dt$mean_inter_prevNext_norm <- 0.5*(all_inter_intra_dt$mean_inter_prevNorm + 
                                                      all_inter_intra_dt$mean_inter_nextNorm)
all_inter_intra_dt$interOverIntra_norm <- all_inter_intra_dt$mean_inter_prevNext_norm/all_inter_intra_dt$mean_intraNorm

# compute inter over intra - matched
all_inter_intra_dt$matched_mean_inter_prevNext_norm <- 0.5*(all_inter_intra_dt$matched_mean_inter_prevNorm + 
                                                              all_inter_intra_dt$matched_mean_inter_nextNorm)
all_inter_intra_dt$matched_interOverIntra_norm <- all_inter_intra_dt$matched_mean_inter_prevNext_norm/all_inter_intra_dt$matched_mean_intraNorm


# compute diff in mean intra
all_inter_intra_dt$meanIntra_norm_diff <- all_inter_intra_dt$mean_intraNorm_matched - all_inter_intra_dt$mean_intraNorm

# compute diff in inter over intra
all_inter_intra_dt$interOverIntra_norm_diff <- all_inter_intra_dt$matched_interOverIntra_norm - all_inter_intra_dt$interOverIntra_norm


# pval vs. inter/intra 
final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
final_table_DT$signif_lab <- ifelse(final_table_DT$adjPvalComb <= 0.01, "signif.", "not signif.")

final_table_DT$hicds_regionID <- file.path(final_table_DT$hicds, final_table_DT$region)

# merge inter intra sur final_table
final_proba_dt <- merge(final_table_DT, all_inter_intra_dt[,c("meanIntra_norm_diff", "interOverIntra_norm_diff")],
                        by=c("hicds", "regionID"))




all_plot_vars <- c("meanIntra_norm_diff", "interOverIntra_norm_diff")


plot_var="mean_intra"

foo <- foreach(plot_var = all_plot_vars) %dopar%{
  
  
  
  myx <- final_proba_dt$meanIntra_norm_diff
  myx <- final_proba_dt$adjPvalComb
  
  densplot(
    
  )
  
  
  
  
  
  
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
  
}
  
  

# hicds    region binStart binEnd   start     end mean_intra
# 1 ENCSR444WCZ_A549_40kb chr1_TAD1       14     26  560001 1080000   107.6835
# 2 ENCSR444WCZ_A549_40kb chr1_TAD2       27     31 1080001 1280000   234.3729
# 3 ENCSR444WCZ_A549_40kb chr1_TAD3       32     36 1280001 1480000   197.8352
# 4 ENCSR444WCZ_A549_40kb chr1_TAD4       37     41 1480001 1680000   213.9768
# 5 ENCSR444WCZ_A549_40kb chr1_TAD5       42     45 1680001 1840000   237.3895
# 6 ENCSR444WCZ_A549_40kb chr1_TAD6       46     51 1840001 2080000   163.9324
# mean_intraNorm mean_inter_prev mean_inter_prevNorm mean_inter_next
# 1       1.915751              NA                  NA        25.16418
# 2       1.790424        25.16418            1.763138        24.38191
# 3       1.259002        24.38191            1.287185        26.91738
# 4       2.194249        26.91738            1.438696        31.28307
# 5       2.075456        31.28307            1.586302        25.87492
# 6       1.883547        25.87492            1.426336        26.17299
# mean_inter_nextNorm matched_mean_intra matched_mean_intraNorm
# 1            1.763138           107.6835               1.915751
# 2            1.287185           234.3729               1.790424
# 3            1.438696           197.8352               1.259002
# 4            1.586302           213.9768               2.194249
# 5            1.426336           237.3895               2.075456
# 6            1.528232           163.9324               1.883547
# matched_mean_inter_prev matched_mean_inter_prevNorm matched_mean_inter_next
# 1                      NA                          NA                25.16418
# 2                       0                           0                24.38191
# 3                       0                           0                26.91738
# 4                       0                           0                31.28307
# 5                       0                           0                25.87492
# 6                       0                           0                26.17299
# matched_mean_inter_nextNorm
# 1                    1.763138
# 2                    1.287185
# 3                    1.438696
# 4                    1.586302
# 
# 
