
# Rscript check_randommidposdisc.R

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "check_randommidposdisc.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"

require(foreach)
require(doMC)
registerDoMC(40)

pipFolder<- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

discFolder <- "CHECK_DIFFTADS_RANDOMMIDPOS"

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[grepl("_RANDOMMIDPOS_", all_hicds) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


# all_hicds = all_hicds[1]
all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  
  hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    new_hicds <- gsub("_RANDOMMIDPOS_", "_RANDOMMIDPOSDISC_", hicds)
  
    logfc <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", new_hicds, exprds, "3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata")))
    reglist <- get(load(file.path("CHECK_DIFFTADS_RANDOMMIDPOS", hicds, exprds, "pipeline_regionList.Rdata")))
    stopifnot(setequal(names(logfc), reglist))
  }
}
