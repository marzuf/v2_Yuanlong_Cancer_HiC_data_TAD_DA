require(ggsci)
require(ggpubr)
require(ggplot2)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# Rscript revision_FC_activDiff.R


outFolder <- "REVISION_FC_ACTIVDIFF_V2_CORRECTED"
dir.create(outFolder, recursive = TRUE)

matching_dt <- get(load("REVISION_PROBADIFFMATCHEDPAIRS_V2_CORRECTED/mean_intraNorm_matching_withRank_dt.Rdata"))

all_inter_intra_dt <- get(load("REVISION_INTER_INTRA_PROBA_V2_CORRECTED_PAIREDHIC//all_inter_intra_dt.Rdata"))


plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.2

all_pairs <- c(
  file.path("LI_40kb","GSE105381_HepG2_40kb", "TCGAlihc_norm_lihc"),
  file.path("LG1_40kb" ,"ENCSR444WCZ_A549_40kb", "TCGAluad_norm_luad"),
  file.path("LG2_40kb" ,"ENCSR444WCZ_A549_40kb" ,"TCGAluad_norm_luad"),
  file.path("LG1_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("LG2_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "TCGAprad_norm_prad"),
  file.path("GSE118514_RWPE1_40kb", "GSE118514_22Rv1_40kb", "TCGAprad_norm_prad")
)
all_normal_ds <- as.character(sapply(all_pairs, function(x) dirname(dirname(x))))
all_tumor_ds <-  as.character(sapply(all_pairs, function(x) basename(dirname(x))))

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
### PREPARE THE GENE FC DATA
###################
gene_tad_signif_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
gene_tad_signif_dt$regionID <- file.path(gene_tad_signif_dt$hicds, gene_tad_signif_dt$exprds,gene_tad_signif_dt$region )
stopifnot(setequal(gene_tad_signif_dt$regionID, final_table_DT$regionID))

# -> assign the fc to gene fc quantile
# -> take the mean  by TAD
# -> signif not signif
# -> for the matched data, meanIntraDiff
gene_tad_signif_dt$dataset <- file.path(gene_tad_signif_dt$hicds, gene_tad_signif_dt$exprds)

gene_tad_histfc_dt <- do.call(rbind, by(gene_tad_signif_dt, gene_tad_signif_dt$dataset, function(sub_dt) {
  
  all_aggFC_hist <- hist(sub_dt$logFC, breaks = seq(min(sub_dt$logFC), max(sub_dt$logFC), 
                                                    length.out=nBreaks+1), plot=FALSE)# $breaks
  # then for each value find in which break it falls
  logFC_hist <- sapply(sub_dt$logFC, function(x) { 
    xbreak <- which(hist(x, breaks = all_aggFC_hist$breaks, plot=FALSE)$counts == 1);
    stopifnot(length(xbreak) == 1); 
    xbreak})
  check_histqt <- factor(logFC_hist, levels = seq_along(all_aggFC_hist$counts))
  stopifnot(table(check_histqt) == all_aggFC_hist$counts)
  stopifnot(max(logFC_hist) <= nBreaks)
  stopifnot(min(logFC_hist) >= 1)
  stopifnot(length(logFC_hist) == nrow(sub_dt))
  sub_dt$logFC_histBreak <- logFC_hist
  sub_dt
}))

outFile <- file.path(outFolder, "gene_tad_histfc_dt.Rdata")
save(gene_tad_histfc_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

stopifnot(nrow(gene_tad_histfc_dt) == nrow(gene_tad_signif_dt))

# stop("-ok\n")

tad_agg_fcHist_dt <- aggregate(logFC_histBreak~hicds+exprds+region+regionID, data=gene_tad_histfc_dt,
                            FUN=mean)
colnames(tad_agg_fcHist_dt)[ colnames(tad_agg_fcHist_dt) == "logFC_histBreak"] <- "meanLogFCbreak"
all_meanLogFCbreaks <- setNames(tad_agg_fcHist_dt$meanLogFCbreak, tad_agg_fcHist_dt$regionID)

###################
### plot ggdensity the fc breaks signif and not signif 
###################

tad_agg_fcHist_dt$signif_lab <- ifelse(tad_agg_fcHist_dt$regionID %in% signif_tads, "signif.", "not signif.")

###################
### merge the data --- matching regions
###################
stopifnot(matching_dt$ref_region_ID %in% names(all_meanLogFCbreaks))
stopifnot(matching_dt$matching_region_ID %in% names(all_meanLogFCbreaks))
matching_dt$ref_meanLogFCbreaks <- all_meanLogFCbreaks[paste0(matching_dt$ref_region_ID)]
matching_dt$matching_meanLogFCbreaks <- all_meanLogFCbreaks[paste0(matching_dt$matching_region_ID)]

# selecte the signif. tads

###################
### mean FC breaks and diff in mean intra norm for matching tad
###################

signif_matching_dt <- matching_dt[matching_dt$ref_region_pval <= tadSignifThresh,]
stopifnot(nrow(signif_matching_dt) > 0)
stopifnot(signif_matching_dt$ref_region_ID %in% signif_tads)

signif_matching_dt$diff_mean_intraNorm <- signif_matching_dt$ref_mean_intraNorm-signif_matching_dt$matching_mean_intraNorm

plot(x=signif_matching_dt$tumorOverNorm_mean_intraNorm,
     y=signif_matching_dt$ref_meanLogFCbreaks)

plot(x=signif_matching_dt$diff_mean_intraNorm,
     y=signif_matching_dt$ref_meanLogFCbreaks)


###################
### merge the data --- paired hi-c
###################

sub_tad_agg_fcHist_dt <- tad_agg_fcHist_dt[(tad_agg_fcHist_dt$hicds %in% all_normal_ds | 
                                              tad_agg_fcHist_dt$hicds %in% all_tumor_ds) & 
                                             tad_agg_fcHist_dt$exprds %in% basename(all_pairs),]

sub_tad_agg_fcHist_dt$hicds_regionID <- file.path(sub_tad_agg_fcHist_dt$hicds, sub_tad_agg_fcHist_dt$region)
all_inter_intra_dt$hicds_regionID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)

matched_paired_dt <- merge(sub_tad_agg_fcHist_dt,all_inter_intra_dt[,c("hicds_regionID", 
                                                                   "mean_intraNorm",
                                                                   "matched_mean_intraNorm")],
      by="hicds_regionID")


matched_paired_dt$diff_mean_intraNorm <- matched_paired_dt$mean_intraNorm-
                                                  matched_paired_dt$matched_mean_intraNorm

plot(x=matched_paired_dt$diff_mean_intraNorm,
     y=matched_paired_dt$meanLogFCbreak)

signif_matched_paired_dt <- matched_paired_dt[matched_paired_dt$regionID %in% signif_tads,]
stopifnot(nrow(signif_matched_paired_dt) > 0)
plot(x=signif_matched_paired_dt$diff_mean_intraNorm,
     y=signif_matched_paired_dt$meanLogFCbreak)

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
