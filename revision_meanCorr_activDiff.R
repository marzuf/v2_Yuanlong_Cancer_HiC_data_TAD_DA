require(ggsci)
require(ggpubr)
require(ggplot2)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# Rscript revision_meanCorr_activDiff.R


outFolder <- "REVISION_MEANCORR_ACTIVDIFF_V2_CORRECTED"
dir.create(outFolder, recursive = TRUE)
matching_dt <- get(load("REVISION_PROBADIFFMATCHEDPAIRS_V2_CORRECTED/mean_intraNorm_matching_withRank_dt.Rdata"))
all_inter_intra_dt <- get(load("REVISION_INTER_INTRA_PROBA_V2_CORRECTED_PAIREDHIC//all_inter_intra_dt.Rdata"))

### THIS CANNOT BE RUN BECAUSE **** REVISION_INTER_INTRA_PROBA_V2_CORRECTED_PAIREDHIC ***
# ONLY AVAILABLE FOR V2
# outFolder <- "REVISION_MEANCORR_ACTIVDIFF_CORRECTED"
# dir.create(outFolder, recursive = TRUE)
# matching_dt <- get(load("REVISION_PROBADIFFMATCHEDPAIRS_CORRECTED/mean_intraNorm_matching_withRank_dt.Rdata"))
# all_inter_intra_dt <- get(load("REVISION_INTER_INTRA_PROBA_CORRECTED_PAIREDHIC//all_inter_intra_dt.Rdata"))


buildTable <- FALSE

plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.2
myHeightGG <- 5
myWidthGG <- 6


source("revision_settings.R")


stopifnot(matching_dt$ref_exprds %in% basename(all_pairs))
stopifnot(matching_dt$matching_exprds %in% basename(all_pairs))

tadSignifThresh <- 0.01

nBreaks <- 10

###################
### PREPARE SIGNIF DATA
###################

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)
signif_tads <- final_table_DT$regionID[final_table_DT$adjPvalComb <= tadSignifThresh]

###################
### PREPARE THE GENE meanCorr DATA -> I don't have gene level info !!!
###################




# gene_tad_signif_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
# gene_tad_signif_dt$regionID <- file.path(gene_tad_signif_dt$hicds, gene_tad_signif_dt$exprds,gene_tad_signif_dt$region )
# stopifnot(setequal(gene_tad_signif_dt$regionID, final_table_DT$regionID))

# -> assign the meanCorr to gene meanCorr quantile
# -> take the mean  by TAD
# -> signif not signif
# -> for the matched data, meanIntraDiff
# gene_tad_signif_dt$dataset <- file.path(gene_tad_signif_dt$hicds, gene_tad_signif_dt$exprds)
# 
# if(buildTable) {
#   gene_tad_histMeanCorr_dt <- do.call(rbind, by(gene_tad_signif_dt, gene_tad_signif_dt$dataset, function(sub_dt) {
#     
#     all_aggMeanCorr_hist <- hist(sub_dt$meanCorr, breaks = seq(min(sub_dt$meanCorr), max(sub_dt$meanCorr), 
#                                                       length.out=nBreaks+1), plot=FALSE)# $breaks
#     # then for each value find in which break it falls
#     meanCorr_hist <- sapply(sub_dt$meanCorr, function(x) { 
#       xbreak <- which(hist(x, breaks = all_aggMeanCorr_hist$breaks, plot=FALSE)$counts == 1);
#       stopifnot(length(xbreak) == 1); 
#       xbreak})
#     check_histqt <- factor(meanCorr_hist, levels = seq_along(all_aggMeanCorr_hist$counts))
#     stopifnot(table(check_histqt) == all_aggMeanCorr_hist$counts)
#     stopifnot(max(meanCorr_hist) <= nBreaks)
#     stopifnot(min(meanCorr_hist) >= 1)
#     stopifnot(length(meanCorr_hist) == nrow(sub_dt))
#     sub_dt$meanCorr_histBreak <- meanCorr_hist
#     sub_dt
#   }))
#   
#   outFile <- file.path(outFolder, "gene_tad_histMeanCorr_dt.Rdata")
#   save(gene_tad_histMeanCorr_dt, file=outFile, version=2)
#   cat(paste0("... written: ", outFile, "\n"))
#   
# } else {
#   outFile <- file.path(outFolder, "gene_tad_histMeanCorr_dt.Rdata")
#   gene_tad_histMeanCorr_dt <- get(load(outFile))
# }
# 
# 
# stopifnot(nrow(gene_tad_histMeanCorr_dt) == nrow(gene_tad_signif_dt))

# stop("-ok\n")

# tad_agg_meanCorrHist_dt <- aggregate(meanCorr_histBreak~hicds+exprds+region+regionID, data=gene_tad_histMeanCorr_dt,
#                             FUN=mean)
# colnames(tad_agg_meanCorrHist_dt)[ colnames(tad_agg_meanCorrHist_dt) == "meanCorr_histBreak"] <- "meanCorrbreak"

cat("prep corr. data")


final_table_DT$dataset <- file.path(final_table_DT$hicds, final_table_DT$exprds)

if(buildTable) {

  tad_agg_meanCorrHist_dt <- do.call(rbind, by(final_table_DT, final_table_DT$dataset, function(sub_dt) {
    
    all_aggMeanCorr_hist <- hist(sub_dt$meanCorr, breaks = seq(min(sub_dt$meanCorr), max(sub_dt$meanCorr), 
                                                               length.out=nBreaks+1), plot=FALSE)# $breaks
    # then for each value find in which break it falls
    meanCorr_hist <- sapply(sub_dt$meanCorr, function(x) { 
      xbreak <- which(hist(x, breaks = all_aggMeanCorr_hist$breaks, plot=FALSE)$counts == 1);
      stopifnot(length(xbreak) == 1); 
      xbreak})
    check_histqt <- factor(meanCorr_hist, levels = seq_along(all_aggMeanCorr_hist$counts))
    stopifnot(table(check_histqt) == all_aggMeanCorr_hist$counts)
    stopifnot(max(meanCorr_hist) <= nBreaks)
    stopifnot(min(meanCorr_hist) >= 1)
    stopifnot(length(meanCorr_hist) == nrow(sub_dt))
    sub_dt$meanCorr_histBreak <- meanCorr_hist
    sub_dt
  }))
  
  outFile <- file.path(outFolder, "tad_agg_meanCorrHist_dt.Rdata")
  save(tad_agg_meanCorrHist_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))

} else {
  outFile <- file.path(outFolder, "tad_agg_meanCorrHist_dt.Rdata")
  tad_agg_meanCorrHist_dt <- get(load(outFile))
}
# 


all_meanCorrbreaks <- setNames(tad_agg_meanCorrHist_dt$meanCorr_histBreak, tad_agg_meanCorrHist_dt$regionID)

cat("done")

###################
### plot ggdensity the meanCorr breaks signif and not signif 
###################

stopifnot(!duplicated(tad_agg_meanCorrHist_dt$regionID))

tad_agg_meanCorrHist_dt$signif_lab <- ifelse(tad_agg_meanCorrHist_dt$regionID %in% signif_tads,
                                       "signif.", "not signif.")

stopifnot(sum(tad_agg_meanCorrHist_dt$signif_lab=="signif.") == 
            sum(final_table_DT$adjPvalComb <= tadSignifThresh))

plotTit <- paste0("")

mySub <- paste0("# DS = ", length(unique(tad_agg_meanCorrHist_dt$hicds, tad_agg_meanCorrHist_dt$exprds)), 
                "; # TADs = ", nrow(tad_agg_meanCorrHist_dt),
                " (# signif.  = ", sum(tad_agg_meanCorrHist_dt$signif_lab=="signif."), ")")


plot_var <- "meanCorr_histBreak"

legTitle <- ""

my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(tad_agg_meanCorrHist_dt$signif_lab))


p3 <- ggdensity(tad_agg_meanCorrHist_dt,
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


###################
### merge the data --- matching regions
###################
stopifnot(matching_dt$ref_region_ID %in% names(all_meanCorrbreaks))
stopifnot(matching_dt$matching_region_ID %in% names(all_meanCorrbreaks))
matching_dt$ref_meanCorrbreaks <- all_meanCorrbreaks[paste0(matching_dt$ref_region_ID)]
matching_dt$matching_meanCorrbreaks <- all_meanCorrbreaks[paste0(matching_dt$matching_region_ID)]

# selecte the signif. tads

###################
### mean meanCorr breaks and diff in mean intra norm for matching tad
###################

matching_dt$diff_mean_intraNorm <- matching_dt$ref_mean_intraNorm-matching_dt$matching_mean_intraNorm

plot_list <- list(
  c(myxvar="tumorOverNorm_mean_intraNorm", myyvar="ref_meanCorrbreaks"),
  c(myxvar = "tumorOverNorm_mean_intraNorm", myyvar = "ref_meanCorrbreaks")
  )

toplots <- c("all", "signif")
toplot="all"
i=1
for(i in seq_along(plot_list)) {
  
  myxvar <- as.character(plot_list[[i]]["myxvar"])
  myyvar <- as.character(plot_list[[i]]["myyvar"])
  
  
  for(toplot in toplots) {
    
    if(toplot == "all") {
      toplot_dt <- matching_dt
    } else if(toplot == "signif"){
      signif_matching_dt <- matching_dt[matching_dt$ref_region_pval <= tadSignifThresh,]
      stopifnot(nrow(signif_matching_dt) > 0)
      stopifnot(signif_matching_dt$ref_region_ID %in% signif_tads)
      toplot_dt <- signif_matching_dt
    } else {
      stop("--error\n")
    }
    
    
    myx <- toplot_dt[,paste0(myxvar)]
    myy <- toplot_dt[,paste0(myyvar)]
    
    plotTit <- paste0( myyvar, " vs. ", myxvar)
    
    mySub <- paste0(toplot, " TADs - ", "# hicds = ", length(unique(toplot_dt$ref_hicds)), 
                    "; # cmps = ", length(unique(file.path(toplot_dt$matching_hicds, 
                                                           toplot_dt$ref_hicds))), 
                    "; # TADs = ", nrow(toplot_dt))
    
    outFile  <- file.path(outFolder, paste0(toplot, "TADs_", myyvar, "_vs_", myxvar, "_matched_densplot.", plotType))
    do.call(plotType, list(outFile, height=myWidth, width=myWidth))
    densplot(x=myx,
             y=myy,
              xlab = paste0(myxvar),
             ylab = paste0(myyvar),
             main  = plotTit,
             plot.main=plotCex,
             plot.main=plotCex,
             plot.main=plotCex
    )
    mtext(side=3, text=mySub)
    addCorr(x=myx,y=myy, legPos="topleft", bty="n", corMet="spearman")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
}








###################
### merge the data --- paired hi-c
###################

sub_tad_agg_meanCorrHist_dt <- tad_agg_meanCorrHist_dt[(tad_agg_meanCorrHist_dt$hicds %in% all_normal_ds | 
                                              tad_agg_meanCorrHist_dt$hicds %in% all_tumor_ds) & 
                                             tad_agg_meanCorrHist_dt$exprds %in% basename(all_pairs),]

sub_tad_agg_meanCorrHist_dt$hicds_regionID <- file.path(sub_tad_agg_meanCorrHist_dt$hicds, sub_tad_agg_meanCorrHist_dt$region)
all_inter_intra_dt$hicds_regionID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)

matched_paired_dt <- merge(sub_tad_agg_meanCorrHist_dt,all_inter_intra_dt[,c("hicds_regionID", 
                                                                   "mean_intraNorm",
                                                                   "matched_mean_intraNorm")],
      by="hicds_regionID")


matched_paired_dt$diff_mean_intraNorm <- matched_paired_dt$mean_intraNorm-
                                                  matched_paired_dt$matched_mean_intraNorm


myxvar <- "diff_mean_intraNorm"
myyvar <- "meanCorr_histBreak"

for(toplot in toplots) {
  
  if(toplot == "all") {
    toplot_dt <- matched_paired_dt
  } else if(toplot == "signif"){
    signif_matched_paired_dt <- matched_paired_dt[matched_paired_dt$regionID %in% signif_tads,]
    stopifnot(nrow(signif_matched_paired_dt) > 0)
    stopifnot(signif_matched_paired_dt$ref_region_ID %in% signif_tads)
    toplot_dt <- signif_matched_paired_dt
  } else {
    stop("--error\n")
  }
  
  
  myx <- toplot_dt[,paste0(myxvar)]
  myy <- toplot_dt[,paste0(myyvar)]
  
  plotTit <- paste0( myyvar, " vs. ", myxvar)
  
  mySub <- paste0(toplot, " TADs - ",
                  # "# hicds = ", length(unique(toplot_dt$ref_hicds)), 
                  # "; # cmps = ", length(unique(file.path(toplot_dt$matching_hicds, 
                  #                                        toplot_dt$ref_hicds))), 
                  "# hicds = ", length(unique(file.path(toplot_dt$hicds, toplot_dt$exprds))), 
                  "; # TADs = ", nrow(toplot_dt))
  
  outFile  <- file.path(outFolder, paste0(toplot, "TADs_", myyvar, "_vs_", myxvar, "_matched_densplot.", plotType))
  do.call(plotType, list(outFile, height=myWidth, width=myWidth))
  densplot(x=myx,
           y=myy,
           xlab = paste0(myxvar),
           ylab = paste0(myyvar),
           main=plotTit,
           plot.main=plotCex,
           plot.main=plotCex,
           plot.main=plotCex
  )
  mtext(side=3, text=mySub)
  addCorr(x=myx,y=myy, legPos="topleft", bty="n", corMet="spearman")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  

}







# 
# 
# 
# 
# 
# pairFolder <- file.path("REVISION_RANKDIFF_ACTIVDIFF/")
# outFile <- file.path(pairFolder, "matching_data.Rdata")
# matching_data <- get(load(file=outFile))
# 
# ds1_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["norm_matching_pval_tadRank_dt"]]))
# unique(ds1_matching_dt$ref_hicds)
# # [1] "LI_40kb"              "LG1_40kb"             "LG2_40kb"             "GSE118514_RWPE1_40kb"
# ds2_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["tumor_matching_pval_tadRank_dt"]]))
# unique(ds2_matching_dt$ref_hicds)
# # [1] "GSE105381_HepG2_40kb"      "ENCSR444WCZ_A549_40kb"     "ENCSR489OCU_NCI-H460_40kb"
# # [4] "ENCSR346DCU_LNCaP_40kb"    "GSE118514_22Rv1_40kb"     
# 
# 
# matching_withRank_dt <- rbind(ds1_matching_dt, ds2_matching_dt)
# rownames(matching_withRank_dt) <- NULL
# 
# stopifnot(matching_withRank_dt$matching_exprds == matching_withRank_dt$ref_exprds )
# 
# matching_withRank_dt$ref_region_ID <- file.path(matching_withRank_dt$ref_hicds,
#                                                 matching_withRank_dt$ref_exprds,
#                                                 matching_withRank_dt$refID
# )
# matching_withRank_dt$ref_region_hicdsID <- file.path(matching_withRank_dt$ref_hicds,
#                                                      matching_withRank_dt$refID
# )
# 
# matching_withRank_dt$matching_region_ID <- file.path(matching_withRank_dt$matching_hicds,
#                                                      matching_withRank_dt$matching_exprds,
#                                                      matching_withRank_dt$matchingID_maxOverlapBp)
# 
# matching_withRank_dt$matching_region_hicdsID <- file.path(matching_withRank_dt$matching_hicds,
#                                                           matching_withRank_dt$matchingID_maxOverlapBp)
# 
# stopifnot(matching_withRank_dt$ref_region_ID %in% names(regionID_pvals))
# matching_withRank_dt$ref_region_pval <- regionID_pvals[paste0(matching_withRank_dt$ref_region_ID)]
# stopifnot(!is.na(matching_withRank_dt$ref_region_pval))
# stopifnot(round(matching_withRank_dt$ref_region_pval, 6) == round(matching_withRank_dt$adjPval, 6))
# 
# stopifnot(matching_withRank_dt$matching_region_ID %in% names(regionID_pvals))
# matching_withRank_dt$matching_region_pval <- regionID_pvals[paste0(matching_withRank_dt$matching_region_ID)]
# stopifnot(!is.na(matching_withRank_dt$matching_region_pval))
# 
# matching_withRank_dt$ref_tadSignif <- ifelse(matching_withRank_dt$adjPval <= tadSignifThresh, "signif.", "not signif.")
# my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(matching_withRank_dt$ref_tadSignif))
# 
# 
# stopifnot(matching_withRank_dt$ref_region_hicdsID %in% names(colVar_values))
# stopifnot(matching_withRank_dt$matching_region_hicdsID %in% names(colVar_values))
# 
# matching_withRank_dt[,paste0("ref_", col_var)] <- colVar_values[matching_withRank_dt$ref_region_hicdsID]
# matching_withRank_dt[,paste0("matching_", col_var)] <- colVar_values[matching_withRank_dt$matching_region_hicdsID]
# 
# matching_withRank_dt <- matching_withRank_dt[
#   !is.na(  matching_withRank_dt[,paste0("ref_", col_var)] ) & !is.na( matching_withRank_dt[,paste0("matching_", col_var)] ),
# ]
# 
# matching_withRank_dt[,paste0("norm_", col_var)] <- ifelse(matching_withRank_dt$ref_hicds %in% all_normal_ds, 
#                                                           matching_withRank_dt[,paste0("ref_", col_var)],
#                                                           ifelse(matching_withRank_dt$matching_hicds %in% all_normal_ds, 
#                                                                  matching_withRank_dt[,paste0("matching_", col_var)],NA))
# stopifnot(!is.na(matching_withRank_dt[,paste0("norm_", col_var)]))
# 
# matching_withRank_dt[,paste0("tumor_", col_var)] <- ifelse(matching_withRank_dt$ref_hicds %in% all_tumor_ds, 
#                                                            matching_withRank_dt[,paste0("ref_", col_var)],
#                                                            ifelse(matching_withRank_dt$matching_hicds %in% all_tumor_ds, 
#                                                                   matching_withRank_dt[,paste0("matching_", col_var)],NA))
# stopifnot(!is.na(matching_withRank_dt[,paste0("tumor_", col_var)]))
# 
# 
# all_cmps <- unique(file.path(matching_withRank_dt$matching_hicds, matching_withRank_dt$matching_exprds,
#                              matching_withRank_dt$ref_hicds, matching_withRank_dt$ref_exprds))
# 
# mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_withRank_dt), 
#                 " (signif.: ", sum(matching_withRank_dt$adjPval <= tadSignifThresh), ")")
# 
# plotTit <- ""
# 
# # plot
# matching_withRank_dt[,paste0("tumorOverNorm_", col_var)] <- matching_withRank_dt[,paste0("tumor_", col_var)] / matching_withRank_dt[,paste0("norm_", col_var)]  
# 
# stopifnot(matching_withRank_dt$ref_region_ID %in% names(exprVar_values))
# matching_withRank_dt[,paste0("ref_", expr_var)] <- exprVar_values[paste0(matching_withRank_dt$ref_region_ID)]
# stopifnot(!is.na(matching_withRank_dt[,paste0("ref_", expr_var)]))
# 
# outFile  <- file.path(outFolder, paste0(col_var, "_ratio_vs_ref_", expr_var, "_densplot.", plotType))
# do.call(plotType, list(outFile, height=myWidth, width=myWidth))
# 
# densplot(
#   x=matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  ,
#   xlab=paste0("tumorOverNorm_", col_var),
#   y=matching_withRank_dt[,paste0("ref_", expr_var)] ,
#   ylab=paste0("ref_", expr_var),
#   cex.main=plotCex,
#   cex.axis=plotCex,
#   cex.lab=plotCex,
#   main=plotTit
# )
# mtext(side=3, text=mySub)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# qts <- quantile(matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  , probs=c(0.05,0.95))
# crop_dt <- matching_withRank_dt[ matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  >= qts[1] &
#                                    matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  <= qts[2] ,]
# 
# mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_withRank_dt), 
#                 " (signif.: ", sum(matching_withRank_dt$adjPval <= tadSignifThresh), "); cropped 0.05-0.95")
# 
# 
# outFile  <- file.path(outFolder, paste0(col_var, "_ratio_vs_ref_", expr_var, "_cropped_densplot.", plotType))
# do.call(plotType, list(outFile, height=myWidth, width=myWidth))
# 
# densplot(
#   x=crop_dt[,paste0("tumorOverNorm_", col_var)]  ,
#   xlab=paste0("tumorOverNorm_", col_var),
#   y=crop_dt[,paste0("ref_", expr_var)] ,
#   ylab=paste0("ref_", expr_var),
#   cex.main=plotCex,
#   cex.axis=plotCex,
#   cex.lab=plotCex,
#   main=plotTit
# )
# mtext(side=3, text=mySub)
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
