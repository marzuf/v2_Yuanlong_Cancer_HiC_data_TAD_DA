

options(scipen=100)

# Rscript ratioDown_ratioFC_05.R

script_name <- "ratioDown_ratioFC_05.R"

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

outFolder <- file.path("RATIODOWN_RATIOFC_05")
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.4
plotType <- "png"
myHeight <- 400
myWidthLeg <- 500
myWidth <- 400

# all_fcc_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern = "all_obs_prodSignedRatio.Rdata", recursive = TRUE, full.names = TRUE)
# all_rd_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern="all_obs_ratioDown.Rdata", recursive=TRUE, full.names = TRUE)
# all_ratioFC_obs_files <- list.files("OBS_TAD_FC_RATIO", pattern="all_obs_ratioFC.Rdata", recursive=TRUE, full.names = TRUE)
# 
# stopifnot(length(all_rd_obs_files) == length(all_fcc_obs_files) )
# stopifnot(length(all_rd_obs_files) == length(all_negFC_obs_files) )
# stopifnot(length(all_rd_obs_files) == length(all_ratioFC_obs_files) )

pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

if(buildTable) {
  
  all_values_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    ds_values <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      fcc_value <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata")))
      

      rd_value <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown/all_obs_ratioDown.Rdata")))
      
      ratioFC_value <- get(load(file.path("OBS_TAD_FC_RATIO", hicds, exprds, "all_obs_ratioFC.Rdata")))
      
      stopifnot(setequal(names(fcc_value), names(rd_value)))
      stopifnot(setequal(names(fcc_value), names(ratioFC_value)))
      stopifnot(length(fcc_value) == length(rd_value))
      stopifnot(length(fcc_value) == length(ratioFC_value))
      
      all_tads <- names(fcc_value)
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        FCCscore = fcc_value[all_tads],
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
  # outFile="RATIODOWN0_RATIOFC/all_values_dt.Rdata"
  outFile <- file.path(outFolder, "all_values_dt.Rdata")
  all_values_dt <- get(load(outFile))
}


outFile <- file.path(outFolder, paste0("ratioFC_dist_ratioDown05.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot(
  density(all_values_dt$ratioFC[all_values_dt$rD == 0.5]),
  main = paste0("ratioFC dist for ratioDown=0.5")
)
mtext(side=3, text = paste0(sum(all_values_dt$rD == 0.5), "/", nrow(all_values_dt), " TADs"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

unique(all_values_dt$rD[all_values_dt$FCCscore == 0])

outFile <- file.path(outFolder, paste0("ratioDown_dist_ratioFC05.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot(
  density(all_values_dt$rD[all_values_dt$ratioFC == 0.5]),
  main = paste0("ratioFC dist for ratioDown=0.5")
)
mtext(side=3, text = paste0(sum(all_values_dt$ratioFC == 0.5), "/", nrow(all_values_dt), " TADs"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


stopifnot( (2*all_values_dt$ratioFC - 1) * (2* all_values_dt$rD- 1) == all_values_dt$FCCscore)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))











