#!/usr/bin/Rscript

# Rscript leDily_quantif.R

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

stopifnot(! (final_dt_withRank$meanLogFC_downRank <= 100 & final_dt_withRank$meanLogFC_upRank <= 100))
final_dt_withRank$concordanceLab <- ifelse(final_dt_withRank$ratioDown >= thresh1 | final_dt_withRank$ratioDown <= thresh2, "concordant_TAD", "discordant_TAD")
final_dt_withRank$activationLab <- ifelse(final_dt_withRank$meanLogFC_upRank <= 100, "activated_TAD", 
                                          ifelse(final_dt_withRank$meanLogFC_downRank <= 100, "repressed_TAD", "nonresp_TAD")) 








script_name <- "smilePlot.R"

options(scipen=100)

startTime <- Sys.time()

suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
require(foreach)
require(doMC)
registerDoMC(40)


# set files and folders
script8_name <- "8cOnlyRatioDownFastSave_runAllDown"
outFolder <- file.path("SMILEPLOT")
dir.create(outFolder, recursive=TRUE)
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

# settings for plotting
thresh1 <- 0.75
thresh2 <- 1-thresh1
nRandomInit <- 100000
keepPermut <- 1000 # because converted as vector -> too big if keep 100000
nRandom <- keepPermut
curr_ratio_type <- "ratioDown"

## TO SET
## hard-coded:
obsCol <- "bisque"
obsColText <- "bisque2"
permutCol <- "mediumseagreen"
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 1028 , 15)
myWidth <- ifelse(plotType == "png", 686, 10)

histBreakStep <- 0.1
stopifnot( (1/histBreakStep) %% 1 == 0 )
#***********************************************************************************

stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
mainFolder <- "."
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"

setDir <-""

all_missing_data <- foreach(hicds = all_hicds) %dopar% {
  missing_data <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
    pipOutFold <- file.path(pipFolder, hicds, exprds)