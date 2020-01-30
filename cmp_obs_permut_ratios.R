

options(scipen=100)

# Rscript cmp_obs_permut_ratios.R

script_name <- "cmp_obs_permut_ratios.R"

require(foreach)
require(doMC)
registerDoMC(40)

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildTable <- TRUE

plotCex <- 1.4
plotType <- "png"
myHeight <- 400
myWidth <- 600

outFolder <- "CMP_OBS_PERMUT_RATIOS"
dir.create(outFolder, recursive = TRUE)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))


if(buildTable){
  
  all_values <- foreach(hicds = all_hicds) %dopar%{
    ds_values <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      all_permut_FCC <- get(load(file.path("../FIGURES_V3_YUANLONG/RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT", hicds, exprds, "ds_all_permut.Rdata")))
      permut_rD <- unlist(c(lapply(all_permut_FCC, function(x) x[["ratioDown_right"]]), lapply(all_permut_FCC, function(x) x[["ratioDown_left"]])))
      permut_ratioFC <- unlist(c(lapply(all_permut_FCC, function(x) x[["ratioFC_right"]]), lapply(all_permut_FCC, function(x) x[["ratioFC_left"]])))
      
      
      obs_ratioFC <- get(load(file.path("OBS_TAD_FC_RATIO", hicds, exprds, "all_obs_ratioFC.Rdata")))
      obs_rD <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown/all_obs_ratioDown.Rdata")))
      
      stopifnot(length(obs_rD) == length(obs_ratioFC))
      stopifnot(setequal(names(obs_rD), names(obs_ratioFC)))
      
      # plot(obs_rD[names(obs_rD)]~obs_ratioFC[names(obs_rD)])
        
      
      list(
        permut_rD=permut_rD,
        permut_ratioFC=permut_ratioFC,
        obs_rD = obs_rD,
        obs_ratioFC=obs_ratioFC
      )
    }
    ds_values
  }
  outFile <- file.path(outFolder, "all_values.Rdata")
  save(all_values, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  # outFile="CMP_OBS_PERMUT_RATIOS/all_values.Rdata"
  outFile <- file.path(outFolder, "all_values.Rdata")
  all_values <- get(load(outFile))
}

all_obs_rD <- unlist(lapply(all_values, function(subl) lapply(subl, function(x) x[["obs_rD"]])))
all_obs_ratioFC <- unlist(lapply(all_values, function(subl) lapply(subl, function(x) x[["obs_ratioFC"]])))

all_permut_rD <- unlist(lapply(all_values, function(subl) lapply(subl, function(x) x[["permut_rD"]])))
all_permut_ratioFC <- unlist(lapply(all_values, function(subl) lapply(subl, function(x) x[["permut_ratioFC"]])))

outFile <- file.path(outFolder, paste0("cmp_all_permut_obs_ratioDown_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(
    obs_ratioDown = all_permut_rD,
    permut_ratioDown = all_obs_rD
  ), 
  plotTit = paste0("all datasets - ratioDown"),
  legPos = "topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("cmp_all_permut_obs_ratioFC_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(
    obs_ratioFC = all_permut_ratioFC,
    permut_ratioFC = all_obs_ratioFC
  ), 
  plotTit = paste0("all datasets - ratioFC"),
  legPos = "topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))







