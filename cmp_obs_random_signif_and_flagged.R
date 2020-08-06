# Rscript cmp_obs_random_signif_and_flagged.R 

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

plotCex <- 1.2

outFolder <- file.path("CMP_OBS_RANDOM_SIGNIF_AND_FLAGGED")
dir.create(outFolder, recursive = TRUE)

all_obs_dt <- get(load(file.path("ALL_PURITYFLAGGED_FINAL/aran/CPE/log10/all_dt.Rdata")))
all_rd_dt <- get(load(file.path("ALL_PURITYFLAGGED_FINAL_RANDOMMIDPOS//aran/CPE/log10/all_dt.Rdata")))
all_rd_dt$rd_type <- gsub(".+_(.+)_40kb", "\\1", dirname(all_rd_dt$dataset))
all_rd_dt$dataset_init <- all_rd_dt$dataset
all_rd_dt$dataset <- gsub("_RANDOM.+_40kb", "_40kb", all_rd_dt$dataset_init)
stopifnot(setequal(all_rd_dt$dataset, all_obs_dt$dataset))
all_rd_dt$dataset_init <- NULL
nall_dt <- merge(all_obs_dt, all_rd_dt, by=c("dataset"), suffixes=c("_obs", "_rd"))

signif_obs_dt <- get(load(file.path("SIGNIF_PURITYFLAGGED_FINAL//aran/CPE/log10/all_dt.Rdata")))
signif_rd_dt <- get(load(file.path("SIGNIF_PURITYFLAGGED_FINAL_RANDOMMIDPOS//aran/CPE/log10/all_dt.Rdata")))
signif_rd_dt$rd_type <- gsub(".+_(.+)_40kb", "\\1", dirname(signif_rd_dt$dataset))
signif_rd_dt$dataset_init <- signif_rd_dt$dataset
signif_rd_dt$dataset <- gsub("_RANDOM.+_40kb", "_40kb", signif_rd_dt$dataset_init)
stopifnot(setequal(signif_rd_dt$dataset, signif_rd_dt$dataset))
signif_rd_dt$dataset_init <- NULL
signif_dt <- merge(signif_obs_dt, signif_rd_dt, by=c("dataset"), suffixes=c("_obs", "_rd"))

all_dt <- nall_dt # because 2nd get(load overwrite the all_dt
all_dt$cmpTypeCol <- all_cols[all_cmps[basename(all_dt$dataset)]]
stopifnot(!is.na(all_dt$cmpTypeCol))
signif_dt$cmpTypeCol <- all_cols[all_cmps[basename(signif_dt$dataset)]]
stopifnot(!is.na(signif_dt$cmpTypeCol))


rd_types <- c("RANDOMMIDPOS", "RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT")
rd=rd_types[1]
for(rd in rd_types) {
  
  curr_dt <- all_dt[all_dt$rd_type == rd,]
  stopifnot(nrow(curr_dt) > 0)
  stopifnot(!duplicated(curr_dt$dataset))
  
  outFile <- file.path(outFolder, paste0("ratioFlagged_", rd, ".svg"))
  cat(paste0(outFile), "\n")
  svg(outFile, height=6, width=6)
  plot(
    x = curr_dt$ratioFlagged_obs,
    y = curr_dt$ratioFlagged_rd,
    main=paste0("ratioFlagged"),
    pch=16,
    cex=0.7, 
    col=curr_dt$cmpTypeCol,
    xlab="observed",
    ylab=rd,
    cex.axis=plotCex,
    cex.lab=plotCex,
    cex.main=plotCex
  )
  curve(1*x, col="darkgrey", add=TRUE)
  addCorr(
    x = curr_dt$ratioFlagged_obs,
    y = curr_dt$ratioFlagged_rd,
    bty="n", legPos = "topleft"
  )
  legend(
    "bottomright",
    legend=names(all_cols),
    pch=16,
    cex=0.7,
    col=all_cols,
    bty="n"
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  curr_dt <- signif_dt[signif_dt$rd_type == rd,]
  stopifnot(nrow(curr_dt) > 0)
  stopifnot(!duplicated(curr_dt$dataset))
  
  outFile <- file.path(outFolder, paste0("ratioSignifFlagged_", rd, ".svg"))
  svg(outFile, height=6, width=6)
  plot(
    x = curr_dt$ratioSignifFlagged_obs,
    y = curr_dt$ratioSignifFlagged_rd,
    main=paste0("ratioSignifFlagged"),
    pch=16,
    cex=0.7, 
    col=curr_dt$cmpTypeCol,
    xlab="observed",
    ylab=rd,
    cex.axis=plotCex,
    cex.lab=plotCex,
    cex.main=plotCex
  )
  curve(1*x, col="darkgrey", add=TRUE)
  addCorr(
    x = curr_dt$ratioSignifFlagged_obs,
    y = curr_dt$ratioSignifFlagged_rd,
    bty="n", legPos = "topleft"
  )
  legend(
    "bottomright",
    legend=names(all_cols),
    pch=16,
    cex=0.7,
    col=all_cols,
    bty="n"
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}




