col_var <- "mean_intra"

require(ggsci)
require(ggpubr)
require(ggplot2)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# Rscript revision_expressionLevel_activDiff.R

outFolder <- "REVISION_EXPRESSIONLEVEL_ACTIVDIFF"
dir.create(outFolder, recursive = TRUE)

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

tadSignifThresh <- 0.01

###################
### PREPARE SIGNIF DATA
###################

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)

###################
### PREPARE PROBA DIFF DATA
###################
all_inter_intra1_dt <- get(load("REVISION_INTER_INTRA_PROBA/all_inter_intra_dt.Rdata"))
all_inter_intra2_dt <- get(load("REVISION_INTER_INTRA_PROBA2/all_inter_intra_dt.Rdata"))
stopifnot(! all_inter_intra1_dt$hicds %in% all_inter_intra2_dt$hicds)
all_inter_intra_dt <- rbind(all_inter_intra1_dt, all_inter_intra2_dt)
stopifnot(final_table_DT$hicds %in% all_inter_intra_dt$hicds)
stopifnot(col_var %in% colnames(all_inter_intra_dt))

all_inter_intra_dt$region_hicdsID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)
stopifnot(!duplicated(all_inter_intra_dt$region_hicdsID))

colVar_values <- setNames(all_inter_intra_dt[,paste0(col_var)], all_inter_intra_dt$region_hicdsID)


###################
### PREPARE EXPR LEVEL DATA
###################

expr_var <- "aggLog10Expr"
all_exprLevel_dt <- get(load(file.path("REVISION_EXPRESSION_LEVEL", paste0(expr_var, "_aggByTAD_mean.Rdata"))))
exprVar_values <- setNames(all_exprLevel_dt[,paste0(expr_var)], all_exprLevel_dt$regionID)

###################
### PREPARE THE pairs_data
###################

pairFolder <- file.path("REVISION_RANKDIFF_ACTIVDIFF/")
outFile <- file.path(pairFolder, "matching_data.Rdata")
matching_data <- get(load(file=outFile))

ds1_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["norm_matching_pval_tadRank_dt"]]))
unique(ds1_matching_dt$ref_hicds)
# [1] "LI_40kb"              "LG1_40kb"             "LG2_40kb"             "GSE118514_RWPE1_40kb"
ds2_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["tumor_matching_pval_tadRank_dt"]]))
unique(ds2_matching_dt$ref_hicds)
# [1] "GSE105381_HepG2_40kb"      "ENCSR444WCZ_A549_40kb"     "ENCSR489OCU_NCI-H460_40kb"
# [4] "ENCSR346DCU_LNCaP_40kb"    "GSE118514_22Rv1_40kb"     


matching_withRank_dt <- rbind(ds1_matching_dt, ds2_matching_dt)
rownames(matching_withRank_dt) <- NULL

stopifnot(matching_withRank_dt$matching_exprds == matching_withRank_dt$ref_exprds )

matching_withRank_dt$ref_region_ID <- file.path(matching_withRank_dt$ref_hicds,
                                                matching_withRank_dt$ref_exprds,
                                                matching_withRank_dt$refID
)
matching_withRank_dt$ref_region_hicdsID <- file.path(matching_withRank_dt$ref_hicds,
                                                     matching_withRank_dt$refID
)

matching_withRank_dt$matching_region_ID <- file.path(matching_withRank_dt$matching_hicds,
                                                     matching_withRank_dt$matching_exprds,
                                                     matching_withRank_dt$matchingID_maxOverlapBp)

matching_withRank_dt$matching_region_hicdsID <- file.path(matching_withRank_dt$matching_hicds,
                                                          matching_withRank_dt$matchingID_maxOverlapBp)

stopifnot(matching_withRank_dt$ref_region_ID %in% names(regionID_pvals))
matching_withRank_dt$ref_region_pval <- regionID_pvals[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_pval))
stopifnot(round(matching_withRank_dt$ref_region_pval, 6) == round(matching_withRank_dt$adjPval, 6))

stopifnot(matching_withRank_dt$matching_region_ID %in% names(regionID_pvals))
matching_withRank_dt$matching_region_pval <- regionID_pvals[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_pval))

matching_withRank_dt$ref_tadSignif <- ifelse(matching_withRank_dt$adjPval <= tadSignifThresh, "signif.", "not signif.")
my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(matching_withRank_dt$ref_tadSignif))


stopifnot(matching_withRank_dt$ref_region_hicdsID %in% names(colVar_values))
stopifnot(matching_withRank_dt$matching_region_hicdsID %in% names(colVar_values))

matching_withRank_dt[,paste0("ref_", col_var)] <- colVar_values[matching_withRank_dt$ref_region_hicdsID]
matching_withRank_dt[,paste0("matching_", col_var)] <- colVar_values[matching_withRank_dt$matching_region_hicdsID]

matching_withRank_dt <- matching_withRank_dt[
  !is.na(  matching_withRank_dt[,paste0("ref_", col_var)] ) & !is.na( matching_withRank_dt[,paste0("matching_", col_var)] ),
]

matching_withRank_dt[,paste0("norm_", col_var)] <- ifelse(matching_withRank_dt$ref_hicds %in% all_normal_ds, 
                                                          matching_withRank_dt[,paste0("ref_", col_var)],
                                                          ifelse(matching_withRank_dt$matching_hicds %in% all_normal_ds, 
                                                                 matching_withRank_dt[,paste0("matching_", col_var)],NA))
stopifnot(!is.na(matching_withRank_dt[,paste0("norm_", col_var)]))

matching_withRank_dt[,paste0("tumor_", col_var)] <- ifelse(matching_withRank_dt$ref_hicds %in% all_tumor_ds, 
                                                           matching_withRank_dt[,paste0("ref_", col_var)],
                                                           ifelse(matching_withRank_dt$matching_hicds %in% all_tumor_ds, 
                                                                  matching_withRank_dt[,paste0("matching_", col_var)],NA))
stopifnot(!is.na(matching_withRank_dt[,paste0("tumor_", col_var)]))


all_cmps <- unique(file.path(matching_withRank_dt$matching_hicds, matching_withRank_dt$matching_exprds,
                             matching_withRank_dt$ref_hicds, matching_withRank_dt$ref_exprds))

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_withRank_dt), 
                " (signif.: ", sum(matching_withRank_dt$adjPval <= tadSignifThresh), ")")

plotTit <- ""

# plot
matching_withRank_dt[,paste0("tumorOverNorm_", col_var)] <- matching_withRank_dt[,paste0("tumor_", col_var)] / matching_withRank_dt[,paste0("norm_", col_var)]  

stopifnot(matching_withRank_dt$ref_region_ID %in% names(exprVar_values))
matching_withRank_dt[,paste0("ref_", expr_var)] <- exprVar_values[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt[,paste0("ref_", expr_var)]))

outFile  <- file.path(outFolder, paste0(col_var, "_ratio_vs_ref_", expr_var, "_densplot.", plotType))
do.call(plotType, list(outFile, height=myWidth, width=myWidth))

densplot(
  x=matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  ,
  xlab=paste0("tumorOverNorm_", col_var),
  y=matching_withRank_dt[,paste0("ref_", expr_var)] ,
  ylab=paste0("ref_", expr_var),
  cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
  main=plotTit
)
mtext(side=3, text=mySub)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

qts <- quantile(matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  , probs=c(0.05,0.95))
crop_dt <- matching_withRank_dt[ matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  >= qts[1] &
                                   matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  <= qts[2] ,]

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_withRank_dt), 
                " (signif.: ", sum(matching_withRank_dt$adjPval <= tadSignifThresh), "); cropped 0.05-0.95")


outFile  <- file.path(outFolder, paste0(col_var, "_ratio_vs_ref_", expr_var, "_cropped_densplot.", plotType))
do.call(plotType, list(outFile, height=myWidth, width=myWidth))

densplot(
  x=crop_dt[,paste0("tumorOverNorm_", col_var)]  ,
  xlab=paste0("tumorOverNorm_", col_var),
  y=crop_dt[,paste0("ref_", expr_var)] ,
  ylab=paste0("ref_", expr_var),
  cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
  main=plotTit
)
mtext(side=3, text=mySub)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))






















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
# 
# require(foreach)
# require(doMC)
# registerDoMC(40)
# require(ggsci)
# require(ggpubr)
# require(ggrepel)
# require(ggplot2)
# require(patchwork)
# require(stringr)
# blank_theme <- theme_minimal()+
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.y = element_blank(),
#     axis.text.x = element_blank(),
#     panel.border = element_blank(),
#     panel.grid=element_blank(),
#     axis.ticks = element_blank(),
#     plot.title=element_text(size=14, face="bold")
#   )
# 
# 
# # Rscript revision_expressionLevel.R
# 
# setDir <- "/media/electron"
# setDir <- ""
# 
# plotType <- "png"
# myHeightGG <- 5
# myWidthGG <- 7
# myHeight <- 400
# myWidth <- 500
# 
# plotCex <- 1.2
# 
# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# 
# outFolder <- "REVISION_EXPRESSION_LEVEL"
# dir.create(outFolder)
# 
# tadSignifThresh <- 0.01
# 
# final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
# stopifnot(file.exists(final_table_file))
# final_dt <- get(load(final_table_file))
# final_table_DT <- final_dt
# final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
# stopifnot(!duplicated(final_table_DT$regionID))
# regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)
# final_table_DT$signif_lab <- ifelse(final_table_DT$adjPvalComb <= tadSignifThresh, "signif.", "not signif.")
# regionID_signif <- setNames(final_table_DT$signif_lab,final_table_DT$regionID )
# 
# ngenes1 <- str_count(final_table_DT$region_genes, pattern=",")+1
# names(ngenes1) <- final_table_DT$regionID
# 
# expr_level_dt <- get(load("REVISION_EXPRESSION_LEVELDATA/all_exprLevel_dt.Rdata"))
# expr_level_dt$regionID <- file.path(expr_level_dt$hicds, expr_level_dt$exprds, expr_level_dt$region)
# 
# ngenes2 <- setNames(as.numeric(table(expr_level_dt$regionID)), names(table(expr_level_dt$regionID)))
# 
# stopifnot(length(setdiff(names(ngenes1), names(ngenes2))) == 0)
# 
# stopifnot(ngenes1 == ngenes2[names(ngenes1)])
# 
# # > okay, now it is checked, in expr_level_dt i have only the genes used in the pipeline, i can aggregate directly this table
# 
# legTitle <- ""
# 
# mytheme <-     theme(
#   # text = element_text(family=fontFamily),
#   panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
#   panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
#   panel.background = element_rect(fill = "transparent"),
#   panel.grid.major.x =  element_blank(),
#   panel.grid.minor.x =  element_blank(),
#   axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
#   axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
#   axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
#   axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
#   plot.title = element_text(hjust=0.5, size = 16, face="bold"),
#   plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
#   legend.title = element_text(face="bold")
# ) 
# 
# all_vars <- c("aggLog10Expr", "zscoreAggExpr", "qqnormAggExpr", "histqt1AggExpr" ,"histqt2AggExpr")
# 
# expr_var <- "aggLog10Expr"
# aggFun <- "mean"
# 
# for(expr_var in all_vars) {
#   
#   agg_dt <- aggregate(as.formula(paste0(expr_var, "~regionID")), data = expr_level_dt, FUN=aggFun)
#   stopifnot(agg_dt$regionID %in% names(regionID_signif))
#   stopifnot(agg_dt$regionID %in% names(regionID_pvals))
#   agg_dt$signif_lab <- regionID_signif[paste0(agg_dt$regionID)]
#   agg_dt$pval <- regionID_pvals[paste0(agg_dt$regionID)]
#   
#   stopifnot(!is.na(agg_dt$signif_lab))
#   
#   plotTit <- paste0("Expression level and TAD signif.: ", expr_var)
#   
#   mySub <- paste0("# DS = ", length(unique(dirname(agg_dt$regionID))), "; # TADs = ", length(unique(agg_dt$regionID)),
#                   " (# signif. = ", sum(agg_dt$signif_lab == "signif."), ")")
#   
#   my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(agg_dt$signif_lab))
#   
#   
#   p3 <- ggdensity(agg_dt,
#                   x = paste0(expr_var),
#                   y = "..density..",
#                   # combine = TRUE,                  # Combine the 3 plots
#                   xlab = paste0(expr_var, " (TAD ", aggFun, ")"),
#                   # add = "median",                  # Add median line.
#                   rug = FALSE,                      # Add marginal rug
#                   color = "signif_lab",
#                   fill = "signif_lab",
#                   palette = "jco"
#   ) +
#     ggtitle(plotTit, subtitle = mySub)+
#     scale_color_manual(values=my_cols)+
#     scale_fill_manual(values=my_cols)  +
#     labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
#     guides(color=FALSE)+
#     scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
#     scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
#     mytheme
#   
#   outFile <- file.path(outFolder, paste0("tad_", aggFun, "_", expr_var, "_signif_notsignif_density.", plotType))
#   ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   plotTit <- paste0("Expression level and TAD adjPval: ", expr_var)
#   
#   mySub <- paste0("# DS = ", length(unique(dirname(agg_dt$regionID))), "; # TADs = ", length(unique(agg_dt$regionID)),
#                   " (# signif. = ", sum(agg_dt$pval <= tadSignifThresh), ")")
#   
#   myx <- agg_dt[,paste0(expr_var)]
#   myy <- -log10(agg_dt$pval)
#   outFile <- file.path(outFolder, paste0("tad_", aggFun, "_", expr_var, "_vs_adjPval.", plotType))
#   do.call(plotType, list(outFile, height=myHeight, width=myHeight))
#   densplot(
#     x= myx,
#     y=myy,
#     xlab =paste0(expr_var, " (TAD ", aggFun, ")"),
#     ylab = "adj. pval [-log10]",
#     cex.main=plotCex,
#     cex.main=plotCex,
#     cex.main=plotCex,
#     main=plotTit
#   )
#   mtext(side=3, text=mySub)
#   addCorr(x=myx,y=myy, legPos="topright", bty="n")
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))
#   
# }
# 
# all_vars <- c("histqt1AggExpr" ,"histqt2AggExpr")
# all_vars <- c("histqt1AggExpr")
# expr_var <- "histqt1AggExpr"
# 
# 
# for(expr_var in all_vars) {
#   
#   # for each TAD, count how many genes in the quantiles (take the ratio)
#   agg_tad_dt <- aggregate(as.formula(paste0("symbol~regionID + ", expr_var)), data=expr_level_dt, FUN=length)
#   colnames(agg_tad_dt)[colnames(agg_tad_dt) == "symbol"] <- "nGenes"
#   stopifnot(agg_tad_dt$regionID %in% names(regionID_signif))
#   agg_tad_dt$signif_lab <- regionID_signif[paste0(agg_tad_dt$regionID)]
#   
#   stopifnot(agg_tad_dt$regionID %in% names(ngenes2))
#   agg_tad_dt$nTotGenes <- ngenes2[paste0(agg_tad_dt$regionID)]
#   agg_tad_dt$ratioGenes <- agg_tad_dt$nGenes/agg_tad_dt$nTotGenes
#   tmpfoo <- aggregate(ratioGenes~regionID, FUN=sum, data=agg_tad_dt)
#   stopifnot(round(tmpfoo$ratioGenes,10) == 1)
#   
#   stopifnot(range(agg_tad_dt$ratioGenes) >= 0)
#   stopifnot(range(agg_tad_dt$ratioGenes) <= 1)
#         
#   # on average, how many genes by quantile ?
#   ## !!! cannot take the mean in aggregate, because I don't have the rows when 0 gene in a given hist break
#   # agg_dt <- aggregate(as.formula(paste0("ratioGenes~signif_lab + ", expr_var)), data=agg_tad_dt, FUN=mean)
#   # colnames(agg_dt)[colnames(agg_dt) == "ratioGenes"] <- "mean_ratioGenes"
#   tmp_agg_dt <- aggregate(as.formula(paste0("ratioGenes~signif_lab + ", expr_var)), data=agg_tad_dt, FUN=sum)
#   colnames(tmp_agg_dt)[colnames(tmp_agg_dt) == "ratioGenes"] <- "sum_ratioGenes"
#   
#   aggfoo <- aggregate(regionID~signif_lab, FUN=function(x)length(unique(x)), data=agg_tad_dt)
#   nTADs <- setNames(aggfoo$regionID, aggfoo$signif_lab)
#   tmp_agg_dt$nTot <- nTADs[paste0(tmp_agg_dt$signif_lab)]
#   
#   tmp_agg_dt$mean_ratioGenes <- tmp_agg_dt$sum_ratioGenes/tmp_agg_dt$nTot
#   
#   stopifnot(round(sum(tmp_agg_dt$mean_ratioGenes[tmp_agg_dt$signif_lab == "not signif."]), 6) == 1)
#   stopifnot(round(sum(tmp_agg_dt$mean_ratioGenes[tmp_agg_dt$signif_lab == "signif."]), 6) == 1)
#   
#   agg_dt <- tmp_agg_dt
#   
#   agg_dt[, paste0(expr_var)] <- factor(as.character(agg_dt[, paste0(expr_var)]), 
#                                        levels = as.character(min(agg_dt[, paste0(expr_var)]): max(agg_dt[, paste0(expr_var)])))
#   stopifnot(!is.na(  agg_dt[, paste0(expr_var)] ))
#   
#   
#   ############################################################## DRAW DOUBLE PIEPLOT
#   
#   plot_dt <- agg_dt[agg_dt$signif_lab == "signif.",]
#   # plot_dt <- plot_dt[order(as.numeric(plot_dt$histqt1AggExpr)),]
#   plot_dt$freq <- 100*plot_dt$mean_ratioGenes
#   plot_dt$freq_rd <- paste0(round(100*plot_dt$mean_ratioGenes,2), "%")
#   
#   p_signif <- ggplot(plot_dt, aes_string(x="1", y="freq", fill=expr_var)) +
#     ggtitle( "signif. TADs")+
#     geom_col() +
#     # geom_text(aes(label = freq_rd), position = position_stack(vjust = 0.5))+
#     geom_text_repel(aes(label = freq_rd), position = position_stack(vjust = 0.5))+
#     coord_polar(theta = "y") + 
#     theme_void() +    
#     labs(fill="hist. break") +
#     blank_theme
# 
#   plot_dt <- agg_dt[agg_dt$signif_lab == "not signif.",]
#   # plot_dt <- plot_dt[order(as.numeric(plot_dt$histqt1AggExpr)),]
#   plot_dt$freq <- 100*plot_dt$mean_ratioGenes
#   plot_dt$freq_rd <- paste0(round(100*plot_dt$mean_ratioGenes,2), "%")
#   
#   p_notsignif <- ggplot(plot_dt, aes_string(x="1", y="freq", fill=expr_var)) +
#     ggtitle( "not signif. TADs")+
#     geom_col() +
#     geom_text_repel(aes(label = freq_rd), position = position_stack(vjust = 0.5))+
#     coord_polar(theta = "y") + 
#     theme_void() +    
#     labs(fill="hist. break") +
#     blank_theme
#   
#   
#   p_tmp <- (p_signif + p_notsignif)
#   
#   stopifnot(nTADs["signif."] == sum(final_table_DT$adjPvalComb<=tadSignifThresh))
#            
#   
#   plot_tit <- "mean ratios"
#   subtit <- paste0("# TADs: ", paste0(names(nTADs), "=", nTADs, collapse="; "), " adjPval<=", tadSignifThresh)
# 
#     
#   p_out <- p_tmp + 
#     plot_layout(guides = 'collect')+
#     plot_annotation(
#       title = plot_tit,
#       subtitle = paste0(subtit) 
#     ) & theme(plot.title=element_text(hjust=0.5, size=14, face = "bold"), 
#               plot.subtitle=element_text(hjust=0.5, size=12, face = "italic"),
#               legend.position = 'top')
#   
#   outFile <- file.path(outFolder,paste0("mean_ratioGenes_byHistBreak_", expr_var, "_signif_notsignif_pieplot.", plotType))
#   ggsave(p_out, filename = outFile, height=myHeightGG, width=myWidthGG*1.5)
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   ############################################################## DRAW BARPLOT
#   legTitle <- paste0("hist. break")
#     
#   ggbar_p <-  ggbarplot(agg_dt, 
#                         y="mean_ratioGenes",
#                         x="signif_lab", 
#                         fill=paste0(expr_var)) + 
#     ggtitle(plotTit, subtitle=mySub)+
#     mytheme +
#     labs(x="" , y ="ratio of TADs", color=paste0(legTitle),fill=paste0(legTitle)) + 
#     theme(
#       axis.text.x = element_text(hjust=1, vjust=0.5,size=10,angle=90)
#     )
#   outFile <- file.path(outFolder,paste0("mean_ratioGenes_byHistBreak_", expr_var, "_barplot.", plotType))
#   ggsave(ggbar_p, filename = outFile, height=myHeightGG, width=myWidthGG)
#   cat(paste0("... written: ", outFile, "\n"))
# }
# 
# 
# 
# 
# 
