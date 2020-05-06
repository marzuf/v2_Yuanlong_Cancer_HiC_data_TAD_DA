# Rscript create_cond_meanFC.R



setDir <- "/media/electron"
setDir <- ""

require(foreach)
require(doMC)
registerDoMC(40)

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
pipFolder <- pipOutFolder

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds)) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


outFolder <- file.path("CREATE_COND_MEANFC")
dir.create(outFolder, recursive = TRUE)

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    settingF <- file.path("PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingF))
    source(settingF)  
    stopifnot(exists("cond1"))
    stopifnot(exists("cond2"))
    fc_file <- file.path(pipFolder,  hicds, exprds, "3_runMeanTADLogFC", "all_meanLogFC_TAD.Rdata")
    stopifnot(file.exists(fc_file))
    all_fc <- get(load(fc_file))
    
    data.frame(
      hicds = hicds,
      exprds = exprds,
      cond1 = cond1,
      cond2 = cond2,
      region = names(all_fc),
      meanFC = as.numeric(all_fc),
      stringsAsFactors = FALSE
    )
    
  }
  hicds_dt
  
}
outFile <- file.path(outFolder, "all_dt.Rdata")
save(all_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# all_dt <- get(load(file.path(outFolder, "all_dt.Rdata")))

