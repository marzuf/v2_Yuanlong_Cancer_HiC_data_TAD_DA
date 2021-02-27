

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

source("revision_settings.R")


all_inter_intra_dt <- get(load("REVISION_INTER_INTRA_PROBA_V2_CORRECTED_PAIREDHIC//all_inter_intra_dt.Rdata"))

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
  
  mySub <- paste0("# hicds = ", length(unique(all_inter_intra_dt$hicds)), 
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


# compute inter over intra - reference
all_inter_intra_dt$mean_inter_prevNext_norm <- 0.5*(all_inter_intra_dt$mean_inter_prevNorm + 
                                                      all_inter_intra_dt$mean_inter_nextNorm)
all_inter_intra_dt$interOverIntra_norm <- all_inter_intra_dt$mean_inter_prevNext_norm/
                                            all_inter_intra_dt$mean_intraNorm

# compute inter over intra - matched
all_inter_intra_dt$matched_mean_inter_prevNext_norm <- 0.5*(all_inter_intra_dt$matched_mean_inter_prevNorm + 
                                                              all_inter_intra_dt$matched_mean_inter_nextNorm)
all_inter_intra_dt$matched_interOverIntra_norm <- all_inter_intra_dt$matched_mean_inter_prevNext_norm/
                                                all_inter_intra_dt$matched_mean_intraNorm


# compute diff in mean intra
all_inter_intra_dt$meanIntra_norm_diff <- all_inter_intra_dt$matched_mean_intraNorm - 
                                                        all_inter_intra_dt$mean_intraNorm

# compute diff in inter over intra
all_inter_intra_dt$interOverIntra_norm_diff <- all_inter_intra_dt$matched_interOverIntra_norm -
                                                            all_inter_intra_dt$interOverIntra_norm


# pval vs. inter/intra 
final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
final_table_DT$signif_lab <- ifelse(final_table_DT$adjPvalComb <= 0.01, "signif.", "not signif.")
final_table_DT$hicds_regionID <- file.path(final_table_DT$hicds, final_table_DT$region)

## keep only the norm tumor comparisons
final_table_DT <- final_table_DT[final_table_DT$exprds %in% basename(all_pairs),]

all_inter_intra_dt$hicds_regionID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)

# merge inter intra sur final_table
final_proba_dt <- merge(final_table_DT, 
                        all_inter_intra_dt[,c("hicds_regionID", "matched_hicds", "meanIntra_norm_diff", "interOverIntra_norm_diff")],
                        by=c("hicds_regionID"))

my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(final_proba_dt$signif_lab))


# stopifnot(!duplicated(final_proba_dt$regionID))
# there are some duplicates e.g. for A549 I have 2 matched hic (LG1 and LG2)

all_plot_vars <- c("meanIntra_norm_diff", "interOverIntra_norm_diff")

stopifnot(all_plot_vars %in% colnames(final_proba_dt))

plot_var="meanIntra_norm_diff"

final_proba_dt$adjPvalComb_log10 <- -log10(final_proba_dt$adjPvalComb)
yvar <- "adjPvalComb_log10"

ncomps <- length(unique(file.path(final_proba_dt$hicds, 
                           final_proba_dt$matched_hicds, final_proba_dt$exprds)))
nhicds <- length(unique(file.path(final_proba_dt$hicds)))

foo <- foreach(plot_var = all_plot_vars) %dopar%{
  
  plotTit <- paste0(yvar, " vs. ", plot_var)
  
  mySub <- paste0("# hicds = ", nhicds, 
                  "; # cmps = ", ncomps,
                  "; # TADs = ", nrow(final_proba_dt))

  myx <- final_proba_dt[, paste0(plot_var)]
  myy <- final_proba_dt[,paste0(yvar)]
  
  outFile  <- file.path(outFolder, paste0(plot_var, "_vs_", yvar, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myWidth, width=myWidth))
  densplot(x = myx,
           xlab =paste0(plot_var),
           y  = myy,
           ylab = paste0(yvar),
           main=plotTit,
           cex.main=plotCex,
           cex.axis =plotCex,
           cex.lab = plotCex)
  mtext(side=3, text=mySub)
  addCorr(x=myx,y=myy, legPos="topleft", bty="n", corMet="spearman")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  plotTit <- paste0(plot_var, " - signif. not signif.")
  

  
  stopifnot(sum(final_proba_dt$signif_lab == "signif.") == 
              sum(final_proba_dt$adjPvalComb <= tadSignifThresh))
  
  mySub <- paste0("# hicds = ", nhicds, 
                    "; # cmps = ", ncomps,
                  "; # TADs = ", nrow(final_proba_dt),
                  " (# signif. = ", sum(final_proba_dt$signif_lab == "signif.") , ")")
  
  legTitle <- ""
  
  p3 <- ggdensity(final_proba_dt,
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
