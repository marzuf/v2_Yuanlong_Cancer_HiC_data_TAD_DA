
# Rscript sampleMeanCorr_fcc.R

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

plotType <- "png"
myHeight <- ifelse(plotType == "png", 400 , 7)
myWidth <- ifelse(plotType == "png", 400, 7)
plotCex <- 1.4

require(flux)
require(foreach)
require(doMC)
registerDoMC(40)

outFolder <- file.path("SAMPLEMEANCORR_FCC")
dir.create(outFolder, recursive=TRUE)

options(scipen=100)

startTime <- Sys.time()

buildData <- TRUE

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAlusc_norm_lusc"

runFolder <-  "."
pipFolder <- file.path(runFolder, "PIPELINE/OUTPUT_FOLDER")

script7_name <- "7sameNbr_runPermutationsMeanTADCorr"
script8_name <- "8cOnlyFCC_runAllDown"

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_obs_exprds <- sapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))



# all_obs_hicds=all_obs_hicds[1]
if(buildData){
  
  all_fcc_corr_dt <- foreach(hicds = all_obs_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_obs_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      meanCorr_samp_file <- file.path(pipFolder, hicds, exprds, script7_name, "meanCorr_sample_around_TADs_sameNbr.Rdata")
      stopifnot(file.exists(meanCorr_samp_file))
      meanCorr_data <- get(load(meanCorr_samp_file))
      
      all_meanCorr_values <- lapply(meanCorr_data, function(x) x[["meanCorr"]])
      
      meanCorr_dt <- data.frame(region=names(all_meanCorr_values), sample_meanCorr = as.numeric(all_meanCorr_values), stringsAsFactors = FALSE)
      meanCorr_dt <- na.omit(meanCorr_dt)
      
      corr_file <- file.path(pipFolder, hicds, exprds, "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(corr_file))
      corr_data <- get(load(corr_file))
      
      fcc_file <- file.path(pipFolder, hicds, exprds, script8_name, "all_obs_prodSignedRatio.Rdata")
      stopifnot(file.exists(fcc_file))
      fcc_data <- get(load(fcc_file))
      
      stopifnot(meanCorr_dt$region %in% names(fcc_data))
      stopifnot(meanCorr_dt$region %in% names(corr_data))
      
      fcc_dt <- data.frame(region=names(fcc_data), fcc = as.numeric(fcc_data), stringsAsFactors = FALSE)
      corr_dt <- data.frame(region=names(corr_data), obsMeanCorr = as.numeric(corr_data), stringsAsFactors = FALSE)
      
      out_dt <- merge(meanCorr_dt, fcc_dt, by=c("region"), all.x=TRUE, all.y=FALSE)
      stopifnot(!is.na(out_dt))
      
      out_dt <- merge(out_dt, corr_dt, by=c("region"), all.x=TRUE, all.y=FALSE)
      stopifnot(!is.na(out_dt))
      
      out_dt$hicds <- hicds
      out_dt$exprds <- exprds
      out_dt
      
    }
    hicds_dt
  }
  
  outFile <- file.path(outFolder, paste0("all_fcc_corr_dt.Rdata"))
  save(all_fcc_corr_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  inFile <- file.path(outFolder, paste0("all_fcc_corr_dt.Rdata"))
  all_dt <- get(load(inFile))
  # load("SAMPLEMEANCORR_FCC/all_fcc_corr_dt.Rdata")
}

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_fcc_corr_dt$corrRatio <- all_fcc_corr_dt$obsMeanCorr/all_fcc_corr_dt$sample_meanCorr
all_fcc_corr_dt$corrDelta <- all_fcc_corr_dt$obsMeanCorr-all_fcc_corr_dt$sample_meanCorr

all_y_vars <- c("sample_meanCorr", "corrRatio", "corrDelta")

x_var <- "fcc"

fcc_thresh <- 0.75

plotCex <- 1.4

for(y_var in all_y_vars){
  
  outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    y=  all_fcc_corr_dt[,y_var],
    x=  all_fcc_corr_dt[,x_var],
    xlab = x_var,
    ylab = y_var,
    cex.main=plotCex,
    cex.lab=plotCex,
    cex.axis=plotCex
  )
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_upper", fcc_thresh, "fcc_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    y=  all_fcc_corr_dt[all_fcc_corr_dt$fcc>fcc_thresh,y_var],
    x=  all_fcc_corr_dt[all_fcc_corr_dt$fcc>fcc_thresh,x_var],
    xlab = x_var,
    ylab = y_var,
    cex.main=plotCex,
    cex.lab=plotCex,
    cex.axis=plotCex
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}




  
