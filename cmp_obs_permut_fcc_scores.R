

options(scipen=100)

# Rscript cmp_obs_permut_fcc_scores.R

script_name <- "cmp_obs_permut_fcc_scores.R"

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

outFolder <- "CMP_OBS_PERMUT_FCC_SCORES"
dir.create(outFolder, recursive = TRUE)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))


if(buildTable){
  
  all_obs_permut_FCC_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    ds_values <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      obs_FCC <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata")))
      
      all_permut_FCC <- get(load(file.path("../FIGURES_V3_YUANLONG/RANDOM_FCC_AUC_RATIO_MEANCORRPERMUT", hicds, exprds, "ds_all_permut.Rdata")))
      permut_FCC <- unlist(lapply(all_permut_FCC, function(x) x[["fcc_meanRL"]]))
      
      stopifnot(length(permut_FCC) == length(obs_FCC))
      stopifnot(setequal(names(permut_FCC), names(obs_FCC)))
      
      data.frame(
        hicds = hicds, 
        exprds = exprds,
        region = names(obs_FCC),
        obs_FCC = obs_FCC[names(obs_FCC)],
        permut_FCC= permut_FCC[names(obs_FCC)],
        stringsAsFactors = FALSE)
      
    }
    ds_values
  }
  outFile <- file.path(outFolder, "all_obs_permut_FCC_dt.Rdata")
  save(all_obs_permut_FCC_dt, file = outFile, version=2)
} else {
  outFile <- file.path(outFolder, "all_obs_permut_FCC_dt.Rdata")
  all_obs_permut_FCC_dt <- get(load(outFile))
}

nDS <- length(unique(file.path(all_obs_permut_FCC_dt$hicds, all_obs_permut_FCC_dt$exprds)))

outFile <- file.path(outFolder, paste0("obsFCC_permutFCC_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(
    obsFCC = all_obs_permut_FCC_dt$obs_FCC,
    permutFCC = all_obs_permut_FCC_dt$permut_FCC
  ), 
  plotTit = "",
  legPos = "topleft"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))







