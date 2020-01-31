

options(scipen=100)

# Rscript check_fcc_scores.R

script_name <- "check_fcc_scores.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildTable <- TRUE

require(flux)
require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)
registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("CHECK_FCC_SCORES")
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.4
plotType <- "png"
myHeight <- 400
myWidthLeg <- 500
myWidth <- 400

all_fcc_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern = "all_obs_prodSignedRatio.Rdata", recursive = TRUE, full.names = TRUE)
all_rd_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern="all_obs_ratioDown.Rdata", recursive=TRUE, full.names = TRUE)
all_negFC_obs_files <- list.files("OBS_TAD_NEGATIVE_FC", pattern="all_obs_negFC.Rdata", recursive=TRUE, full.names = TRUE)
all_ratioFC_obs_files <- list.files("OBS_TAD_FC_RATIO", pattern="all_obs_ratioFC.Rdata", recursive=TRUE, full.names = TRUE)

stopifnot(length(all_rd_obs_files) == length(all_fcc_obs_files) )
stopifnot(length(all_rd_obs_files) == length(all_negFC_obs_files) )
stopifnot(length(all_rd_obs_files) == length(all_ratioFC_obs_files) )

pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

if(buildTable) {
  
  all_values_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    ds_values <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      fcc_value <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata")))
      
      negFC_value <- get(load(file.path("OBS_TAD_NEGATIVE_FC", hicds, exprds, "all_obs_negFC.Rdata")))
      
      rd_value <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown/all_obs_ratioDown.Rdata")))
      
      ratioFC_value <- get(load(file.path("OBS_TAD_FC_RATIO", hicds, exprds, "all_obs_ratioFC.Rdata")))
      
      stopifnot(setequal(names(fcc_value), names(negFC_value)))
      stopifnot(setequal(names(fcc_value), names(rd_value)))
      stopifnot(setequal(names(fcc_value), names(ratioFC_value)))
      stopifnot(length(fcc_value) == length(negFC_value))
      stopifnot(length(fcc_value) == length(rd_value))
      stopifnot(length(fcc_value) == length(ratioFC_value))
      
      all_tads <- names(fcc_value)
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        FCCscore = fcc_value[all_tads],
        negFC = negFC_value[all_tads],
        ratioFC = ratioFC_value[all_tads],
        rD = rd_value[all_tads],
        stringsAsFactors = FALSE
      )
    } # end foreach iterating exprds
    ds_values
  } # end foreach iterating hicds
  
  outFile <- file.path(outFolder, "all_values_dt.Rdata")
  save(all_values_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  # outFile="PIPELINE/OUTPUT_FOLDER/all_values_dt.Rdata"
  outFile <- file.path(outFolder, "all_values_dt.Rdata")
  all_values_dt <- get(load(outFile))
}

all_values_dt$abs_negFC <- abs(all_values_dt$negFC)

rbPal <- colorRampPalette(c('red','blue'))

all_values_dt$dotCols <- rev(rbPal(10))[as.numeric(cut(all_values_dt$FCCscore,breaks = 10))]

nDS <- length(unique(file.path(all_values_dt$hicds, all_values_dt$exprds)))

outFile <- file.path(outFolder, paste0("ratioDown_ratioFC_colFCC.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidthLeg))

par(xpd = T, mar = par()$mar + c(0,0,0,6))

plot(
  x = all_values_dt$ratioFC,
  y = all_values_dt$rD,
  xlab="ratioNegFC",
  ylab="ratioDown",
  col = all_values_dt$dotCols,
  pch = 16,
  cex = 0.7,
  cex.lab = plotCex,
  cex.axis = plotCex,
  main = "all datasets"
)
mtext(side=3, text = paste0("# DS = ", nDS, "; # TADs = ", nrow(all_values_dt)))

# legend("bottomright",
       legend(1.1,1,
       title="FCC score",
       legend=rev(levels(cut(all_values_dt$FCCscore,breaks = 10))),
       col =rbPal(10),
       pch=20,
       cex = 0.8,
       ncol=1,
       horiz = F,
       bty="n")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("ratioNegFC_FCCscore_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = all_values_dt$ratioFC,
  y = all_values_dt$FCCscore,
  xlab="ratioNegFC",
  ylab="FCC score",
  # col = all_values_dt$dotCols,
  pch = 16,
  cex = 0.7,
  cex.lab = plotCex,
  cex.axis = plotCex,
  main = "all datasets"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("ratioDown_FCCscore_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = all_values_dt$rD,
  y = all_values_dt$FCCscore,
  xlab="ratioDown",
  ylab="FCC score",
  # col = all_values_dt$dotCols,
  pch = 16,
  cex = 0.7,
  cex.lab = plotCex,
  cex.axis = plotCex,
  main = "all datasets"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


stopifnot( (2*all_values_dt$ratioFC - 1) * (2* all_values_dt$rD- 1) == all_values_dt$FCCscore)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))












