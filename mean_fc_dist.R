

options(scipen=100)

# Rscript mean_fc_dist.R

script_name <- "mean_fc_dist.R"

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

outFolder <- file.path("MEAN_FC_DIST")
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.4
plotType <- "png"
myHeight <- 400
myWidth <- 600



pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

if(buildTable) {
  
  all_values_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    ds_values <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      meanFC_value <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata")))
      rd_value <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown/all_obs_ratioDown.Rdata")))
      
      stopifnot(setequal(names(rd_value), names(meanFC_value)))
      stopifnot(length(rd_value) == length(meanFC_value))
      
      all_tads <- names(meanFC_value)
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        meanFC = meanFC_value[all_tads],
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



source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
# source("cancer_cols.R")
all_values_dt$cmp_type <- as.character(all_cmps[paste0(all_values_dt$exprds)])
stopifnot(!is.na(all_values_dt$cmp_type))


outFile <- file.path(outFolder, paste0("meanFC_dist_allDS_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot_multiDens(
  split(x = all_values_dt$meanFC, all_values_dt$cmp_type),
  plotTit = paste0("all DS - meanFC dist.")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("ratioDown_dist_allDS_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot_multiDens(
  split(x = all_values_dt$rD, all_values_dt$cmp_type),
  plotTit = paste0("all DS - ratioDown dist.")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
