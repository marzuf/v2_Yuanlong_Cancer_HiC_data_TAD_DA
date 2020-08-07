#!/usr/bin/Rscript

# Rscript leDily_quantif.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(reshape2)

plotCex <- 1.2

outFolder <- file.path("LEDILY_QUANTIF")
dir.create(outFolder, recursive = TRUE)

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
gene_tad_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

# A)
# 100 top and bottom TADs by meanLogFC ("activated_TAD", "repressed_TAD" "nonresp_TAD")
# % of signif. DE genes in activated and repressed
# % of of all genes in activated and repressed

# B)
# % of signif up/down DE genes in activated/repressed TADs
# % of signif in nonresp TADs

# C)
# % of TADs with ratioDown >= 75% and <= 25% ("concordant_TAD")
# % of DE genes in concordant

final_dt$dataset <- file.path(final_dt$hicds, final_dt$exprds)

final_dt_withRank <- do.call(rbind, by(final_dt, final_dt$dataset, function(x) {
  x$meanLogFC_upRank <- rank(x$meanLogFC)
  x$meanLogFC_downRank <- rank(-x$meanLogFC)
  x }))

thresh1 <- 0.75
thresh2 <- 1-thresh1

geneSignifThresh <- 0.01
tadSignifThresh <- 0.01

gene_tad_dt$geneDElab <- ifelse(gene_tad_dt$adj.P.Val > geneSignifThresh, "notDE_gene",
                             ifelse(gene_tad_dt$logFC > 0, "upDE_gene", 
                                    ifelse(gene_tad_dt$logFC < 0, "downDE_gene", NA)))
stopifnot(!is.na(gene_tad_dt$geneDElab))

stopifnot(! (final_dt_withRank$meanLogFC_downRank <= 100 & final_dt_withRank$meanLogFC_upRank <= 100))
final_dt_withRank$concordanceLab <- ifelse(final_dt_withRank$ratioDown >= thresh1 | final_dt_withRank$ratioDown <= thresh2, 
                                           "concordant_TAD", "discordant_TAD")

final_dt_withRank$concordanceDirLab <- ifelse(final_dt_withRank$ratioDown >= thresh1,  "concordantDown_TAD",
                                              ifelse(final_dt_withRank$ratioDown <= thresh2,  "concordantUp_TAD", "discordant_TAD"))
stopifnot(
  which(final_dt_withRank$concordanceDirLab == "discordant_TAD") == which(final_dt_withRank$concordanceLab == "discordant_TAD")
)


final_dt_withRank$activationLab <- ifelse(final_dt_withRank$meanLogFC_upRank <= 100, "activated_TAD", 
                                          ifelse(final_dt_withRank$meanLogFC_downRank <= 100, "repressed_TAD", "nonresp_TAD")) 
# final_dt_withRank$dirLab <- ifelse(final_dt_withRank$meanLogFC > 0, "up", "down")
# stopifnot(!is.na(final_dt_withRank$dirLab))


final_dt_withRank$tadDAlab <- ifelse(final_dt_withRank$adjPvalComb > tadSignifThresh, "notDA_TAD",
                                ifelse(final_dt_withRank$meanLogFC > 0, "upDA_TAD", 
                                       ifelse(final_dt_withRank$meanLogFC < 0, "downDA_TAD", NA)))
stopifnot(!is.na(final_dt_withRank$tadDAlab))


all_merged_dt <- merge(gene_tad_dt[,c("hicds", "exprds", "region", "entrezID", "geneDElab")], 
                       final_dt_withRank[,c("hicds", "exprds", "region", "activationLab", "concordanceLab", "tadDAlab", "concordanceDirLab")], 
                       by=c("hicds", "exprds", "region"), all = TRUE)
stopifnot(!is.na(all_merged_dt))


all_merged_dt$dataset <- file.path(all_merged_dt$hicds, all_merged_dt$exprds)

x=all_merged_dt[all_merged_dt$dataset==all_merged_dt$dataset[1],]
save(x, file="x.Rdata", version=2)

all_stats_dt <- do.call(rbind, by(all_merged_dt, all_merged_dt$dataset, function(x) {
  
  ds <- unique(x$dataset)
  stopifnot(length(ds) == 1)
  
  # A)
  # 100 top and bottom TADs by meanLogFC ("activated_TAD", "repressed_TAD" "nonresp_TAD")
  # % of signif. DE genes in activated and repressed
  signifDE_and_resp <- mean(x$geneDElab != "notDE_gene" & x$activationLab != "nonresp_TAD")
  # % of of all genes in activated and repressed
  all_and_resp <- mean(x$activationLab != "nonresp_TAD")
  
  # % of signif. DE genes in DA
  signifDE_and_DA <- mean(x$geneDElab != "notDE_gene" & x$tadDAlab != "notDA_TAD")
  # % of of all genes in activated and repressed
  all_and_DA <- mean(x$tadDAlab != "notDA_TAD")
  
  # B)
  # % of signif up/down DE genes in activated/repressed TADs
  # => % of upDE_and_activated | downDE_and_repressed
  upDE_and_activated <- mean(x$geneDElab == "upDE_gene" & x$activationLab == "activated_TAD")
  upDE_and_upDA <- mean(x$geneDElab == "upDE_gene" & x$tadDAlab == "upDA_TAD")
  downDE_and_repressed <- mean(x$geneDElab == "downDE_gene" & x$activationLab == "repressed_TAD")
  downDE_and_downDA <- mean(x$geneDElab == "downDE_gene" & x$tadDAlab == "downDA_TAD")
  upDEdownDE_and_activrepr <- mean(
    (x$geneDElab == "upDE_gene" & x$activationLab == "activated_TAD") |
      (x$geneDElab == "downDE_gene" & x$activationLab == "repressed_TAD"))
  upDEdownDE_and_upDAdownDA <- mean(
    (x$geneDElab == "upDE_gene" & x$tadDAlab == "upDA_TAD") |
      (x$geneDElab == "downDE_gene" & x$tadDAlab == "downDA_TAD"))
  if(upDEdownDE_and_upDAdownDA != upDE_and_upDA+downDE_and_downDA) {
    cat("ds=",ds,"\n")
    cat("upDEdownDE_and_upDAdownDA=\t", upDEdownDE_and_upDAdownDA, "\n")
    cat("upDE_and_upDA=\t", upDE_and_upDA, "\n")
    cat("downDE_and_downDA=\t", downDE_and_downDA, "\n")  
  }
  
  # cat("upDEdownDE_and_activrepr=\t", upDEdownDE_and_activrepr, "\n")
  # cat("upDE_and_activated=\t", upDEdownDE_and_upDAdownDA, "\n")
  # cat("downDE_and_repressed=\t", downDE_and_repressed, "\n")
  
  stopifnot(abs(upDEdownDE_and_upDAdownDA- (upDE_and_upDA+downDE_and_downDA)) <= 10^-4)
  stopifnot(abs(upDEdownDE_and_activrepr - (upDE_and_activated+downDE_and_repressed)) <= 10^-4)
  
  # % of signif in nonresp TADs
  signifDE_and_notResp <- mean(x$geneDElab != "notDE_gene" & x$activationLab == "nonresp_TAD")
  signifDE_and_notDA <- mean(x$geneDElab != "notDE_gene" & x$tadDAlab == "notDA_TAD")
  
  # C)
  # % of TADs with ratioDown >= 75% and <= 25% ("concordant_TAD")
  # % of DE genes in concordant
  x_t <- x[,c("dataset", "region", "concordanceLab", "activationLab", "concordanceDirLab")]
  x_t <- unique(x_t)
  concordantTAD <- mean(x_t$concordanceLab == "concordant_TAD")
  signifDE_and_concord <- mean(x$geneDElab != "notDE_gene" & x$concordanceLab == "concordant_TAD")
  all_and_concord <- mean(x$concordanceLab == "concordant_TAD")
  
  # Among TAD genes:
  #   => ratio upDE in activated
  # => ratio downDE in repressed
  # => ratio upDE in nonresp.
  # => ratio downDE in nonresp.
  
  ratioGenes_inTAD_dt <- do.call(rbind, by(x, x$region, function(tad) {
    tadlab <- unique(tad$activationLab)
    stopifnot(length(tadlab) == 1)
    ratioUpDE <- mean(tad$geneDElab == "upDE_gene")
     ratioDownDE <- mean(tad$geneDElab == "downDE_gene")
    data.frame(tadlab=tadlab,
               ratioUpDE=ratioUpDE,
               ratioDownDE=ratioDownDE,
               stringsAsFactors = FALSE)
  }))
  upDE_in_activatedTAD <- mean(ratioGenes_inTAD_dt$ratioUpDE[ratioGenes_inTAD_dt$tadlab=="activated_TAD"])
  downDE_in_repressedTAD <- mean(ratioGenes_inTAD_dt$ratioDownDE[ratioGenes_inTAD_dt$tadlab=="repressed_TAD"])
  upDE_in_nonrespTAD <- mean(ratioGenes_inTAD_dt$ratioUpDE[ratioGenes_inTAD_dt$tadlab=="nonresp_TAD"])
  downDE_in_nonrespTAD <- mean(ratioGenes_inTAD_dt$ratioDownDE[ratioGenes_inTAD_dt$tadlab=="nonresp_TAD"])
  
  
  # Among all TADs:
  #   => % with ratioDown >= 0.75 in repressed and <= 0.25 in activated
 concordantDown_and_repressed <- mean(x_t$concordanceDirLab == "concordantDown_TAD" & x_t$activationLab == "repressed_TAD")
 concordantUp_and_activated <- mean(x_t$concordanceDirLab == "concordantUp_TAD" & x_t$activationLab == "activated_TAD")
  
  
  data.frame(
    dataset = ds,
    signifDE_and_resp=signifDE_and_resp,
    all_and_resp=all_and_resp,
    signifDE_and_DA=signifDE_and_DA,
    all_and_DA=all_and_DA,
    upDE_and_activated=upDE_and_activated,
    upDE_and_upDA=upDE_and_upDA,
    downDE_and_repressed=downDE_and_repressed,
    downDE_and_downDA=downDE_and_downDA,
    upDEdownDE_and_activrepr=upDEdownDE_and_activrepr,
    upDEdownDE_and_upDAdownDA=upDEdownDE_and_upDAdownDA,
    signifDE_and_notResp=signifDE_and_notResp,
    signifDE_and_notDA=signifDE_and_notDA,
    concordantTAD=concordantTAD,
    signifDE_and_concord=signifDE_and_concord,
    all_and_concord=all_and_concord,
    concordantDown_and_repressed=concordantDown_and_repressed,
    concordantUp_and_activated=concordantUp_and_activated,
    upDE_in_activatedTAD=upDE_in_activatedTAD,
    downDE_in_repressedTAD=downDE_in_repressedTAD,
    upDE_in_nonrespTAD=upDE_in_nonrespTAD,
    downDE_in_nonrespTAD=downDE_in_nonrespTAD,
    stringsAsFactors=FALSE
  )
  
  
  
  
}))

rownames(all_stats_dt) <- NULL

outFile <- file.path(outFolder, "all_stats_dt.Rdata")
save(all_stats_dt, file=outFile,version=2)
cat(paste0("... written: ", outFile,"\n"))


# boxplot(NUMS ~ GRP, data = ddf, lwd = 2, ylab = 'NUMS')
# stripchart(NUMS ~ GRP, vertical = TRUE, data = ddf, 
#            method = "jitter", add = TRUE, pch = 20, col = 'blue')

plot_dt <- all_stats_dt[,c("dataset", "signifDE_and_resp", "all_and_resp")]
plot_dt_m <- melt(plot_dt,id="dataset")
par(bty="L")
boxplot(value ~ variable, data = plot_dt_m, lwd = 2,
        cex.main=plotCex, cex.lab=plotCex,cex.axis=plotCex,
        xlab="", ylab = 'average ratio of genes')
stripchart(value ~ variable, vertical = TRUE, data = plot_dt_m, 
            method = "jitter", add = TRUE, pch = 20, col = 'blue')


plot_density_2_vars <- function(var1, var2, data, saveF=TRUE, ...) {
  stopifnot(var1 %in% colnames(data))
  stopifnot(var2 %in% colnames(data))
  plot_list <- list(data[,var1],data[,var2])
  names(plot_list) <- c(var1, var2)
  if(saveF) outFile <- file.path(outFolder, paste0(var1,"_AND_", var2, ".png"))
  if(saveF) png(outFile, height=350, width=550)
  plot_multiDens(plot_list, ...)
  if(saveF) foo <- dev.off()
}

plot_density_2_vars("signifDE_and_resp", "all_and_resp",all_stats_dt, my_xlab= "ratio of genes", legPos = "topleft")

stop("-ok\n")

outFile <- file.path(outFolder, "signifDE_and_resp_AND_all_and_resp.png")
png(outFile, height=350, width=550)
plot_multiDens(
  list(
    signifDE_and_resp = all_stats_dt$signifDE_and_resp,
    all_and_resp = all_stats_dt$all_and_resp
  ), legPos = "topleft", my_xlab="ratio of genes"
)
foo <- dev.off()

outFile <- file.path(outFolder, "signifDE_and_DA_AND_all_and_DA.png")
png(outFile, height=350, width=550)
plot_multiDens(
  list(
    signifDE_and_DA = all_stats_dt$signifDE_and_DA,
    all_and_DA = all_stats_dt$all_and_DA
  ), legPos = "topright", my_xlab="ratio of genes"
)
foo <- dev.off()

outFile <- file.path(outFolder, "upDE_and_activated_AND_upDE_and_upDA.png")
png(outFile, height=350, width=550)
plot_multiDens(
  list(
    upDE_and_activated = all_stats_dt$upDE_and_activated,
    upDE_and_upDA = all_stats_dt$upDE_and_upDA
  ), legPos = "topright", my_xlab="ratio of genes"
)
foo <- dev.off()

outFile <- file.path(outFolder, "downDE_and_repressed_AND_downDE_and_downDA.png")
png(outFile, height=350, width=550)
plot_multiDens(
  list(
    downDE_and_repressed = all_stats_dt$downDE_and_repressed,
    downDE_and_downDA = all_stats_dt$downDE_and_downDA
  ), legPos = "topright", my_xlab="ratio of genes"
)
foo <- dev.off()

outFile <- file.path(outFolder, "signifDE_and_notResp_AND_signifDE_and_notDA.png")
png(outFile, height=350, width=550)
plot_multiDens(
  list(
    signifDE_and_notResp = all_stats_dt$signifDE_and_notResp,
    signifDE_and_notDA = all_stats_dt$signifDE_and_notDA
  ), legPos = "topright", my_xlab="ratio of genes"
)
foo <- dev.off()

outFile <- file.path(outFolder, "signifDE_and_concord_AND_all_and_concord.png")
png(outFile, height=350, width=550)
plot_multiDens(
  list(
    signifDE_and_concord = all_stats_dt$signifDE_and_concord,
    all_and_concord = all_stats_dt$all_and_concord
  ), legPos = "topleft", my_xlab="ratio of genes"
)
foo <- dev.off()

outFile <- file.path(outFolder, "concordantTAD.png")
png(outFile, height=350, width=550)
plot_multiDens(
  list(
    concordantTAD = all_stats_dt$concordantTAD
  ), legPos = "topright", my_xlab="ratio of TADs"
)
foo <- dev.off()

outFile <- file.path(outFolder, "concordance_and_activation.png")
png(outFile, height=350, width=550)
plot_multiDens(
  list(
    concordantDown_and_repressed = all_stats_dt$concordantDown_and_repressed,
    concordantUp_and_activated=all_stats_dt$concordantUp_and_activated
  ), legPos = "topright", my_xlab="ratio of TADs"
)
foo <- dev.off()


outFile <- file.path(outFolder, "geneDE_within_TADs.png")
png(outFile, height=350, width=550)
plot_multiDens(
  list(
    upDE_in_activatedTAD=all_stats_dt$upDE_in_activatedTAD,
    downDE_in_repressedTAD=all_stats_dt$downDE_in_repressedTAD,
    upDE_in_nonrespTAD=all_stats_dt$upDE_in_nonrespTAD,
    downDE_in_nonrespTAD=all_stats_dt$downDE_in_nonrespTAD
  ), legPos = "topright", my_xlab="ratio of genes in TAD"
)
foo <- dev.off()

# 
# 
# 
# 
# 
# script_name <- "smilePlot.R"
# 
# options(scipen=100)
# 
# startTime <- Sys.time()
# 
# suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
# require(foreach)
# require(doMC)
# registerDoMC(40)
# 
# 
# # set files and folders
# script8_name <- "8cOnlyRatioDownFastSave_runAllDown"
# outFolder <- file.path("SMILEPLOT")
# dir.create(outFolder, recursive=TRUE)
# pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
# settingFolder <- file.path("PIPELINE", "INPUT_FILES")
# 
# # settings for plotting
# thresh1 <- 0.75
# thresh2 <- 1-thresh1
# nRandomInit <- 100000
# keepPermut <- 1000 # because converted as vector -> too big if keep 100000
# nRandom <- keepPermut
# curr_ratio_type <- "ratioDown"
# 
# ## TO SET
# ## hard-coded:
# obsCol <- "bisque"
# obsColText <- "bisque2"
# permutCol <- "mediumseagreen"
# plotType <- "svg"
# myHeight <- ifelse(plotType == "png", 1028 , 15)
# myWidth <- ifelse(plotType == "png", 686, 10)
# 
# histBreakStep <- 0.1
# stopifnot( (1/histBreakStep) %% 1 == 0 )
# #***********************************************************************************
# 
# stopifnot(dir.exists(pipFolder))
# all_hicds <- list.files(pipFolder)
# all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
# mainFolder <- "."
# stopifnot(dir.exists(file.path(mainFolder, all_hicds)))
# all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
# names(all_exprds) <- all_hicds
# 
# hicds <- "ENCSR489OCU_NCI-H460_40kb"
# exprds <- "TCGAluad_norm_luad"
# 
# setDir <-""
# 
# all_missing_data <- foreach(hicds = all_hicds) %dopar% {
#   missing_data <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
#     pipOutFold <- file.path(pipFolder, hicds, exprds)