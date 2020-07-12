
# Rscript cmp_partialcorr.R

require(foreach)
require(doMC)
registerDoMC(40)

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

outFolder <- "CMP_PARTIALCORR"
dir.create(outFolder, recursive=T)

script4_name <- "4_runMeanTADCorr"
script4partial_name <- "4partial_runMeanTADCorr"

script11_name <- "11sameNbr_runEmpPvalCombined"
script11partial_name <- "11sameNbrPartial_runEmpPvalCombined"

# all_ds=all_ds[1]

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"


# all_ds = file.path(ex_hicds, ex_exprds)
  
all_full_partial_corr_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  cat(paste0("... start: ", ds, "\n"))
  hicds <- dirname(ds)
  exprds <- basename(ds)
  empPval <- get(load(file.path(pipFolder, hicds, exprds, script11_name, "emp_pval_combined.Rdata")))
  partialEmpPval <- get(load(file.path(pipFolder, hicds, exprds, script11partial_name, "emp_pval_combined.Rdata")))
  stopifnot(names(empPval) == names(partialEmpPval))
  
  meanCorr <- get(load(file.path(pipFolder, hicds, exprds,script4_name, "all_meanCorr_TAD.Rdata")))
  partialmeanCorr <- get(load(file.path(pipFolder, hicds, exprds, script4partial_name, "all_meanCorr_TAD.Rdata")))
  stopifnot(names(meanCorr) == names(partialmeanCorr))
  
  
  data.frame(
    dataset =ds,  
    empPval_full_partial_corr = cor(partialEmpPval, empPval),
    meanCorr_full_partial_corr = cor(partialmeanCorr, meanCorr),
    stringsAsFactors = FALSE
  )
   
}
outFile <- file.path(outFolder, "all_full_partial_corr_dt.Rdata")
save(all_full_partial_corr_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

all_full_partial_delta_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  cat(paste0("... start: ", ds, "\n"))
  hicds <- dirname(ds)
  exprds <- basename(ds)
  empPval <- get(load(file.path(pipFolder, hicds, exprds, script11_name, "emp_pval_combined.Rdata")))
  partialEmpPval <- get(load(file.path(pipFolder, hicds, exprds, script11partial_name, "emp_pval_combined.Rdata")))
  stopifnot(names(empPval) == names(partialEmpPval))
  
  meanCorr <- get(load(file.path(pipFolder, hicds, exprds,script4_name, "all_meanCorr_TAD.Rdata")))
  partialmeanCorr <- get(load(file.path(pipFolder, hicds, exprds, script4partial_name, "all_meanCorr_TAD.Rdata")))
  stopifnot(names(meanCorr) == names(partialmeanCorr))
  
  stopifnot(setequal(names(meanCorr) , names(empPval)))
  
  all_regs <- names(meanCorr)
  
  out_dt <- data.frame(
    dataset =ds,  
    region = all_regs,
    empPval_full_partial_delta = as.numeric(empPval[all_regs])-as.numeric(partialEmpPval[all_regs]),
    absMeanCorr_full_partial_delta = abs(as.numeric(meanCorr[all_regs]))-abs(as.numeric(partialmeanCorr[all_regs])),
    stringsAsFactors = FALSE
  )
  stopifnot(!is.na(out_dt))
  out_dt
}
outFile <- file.path(outFolder, "all_full_partial_delta_dt.Rdata")
save(all_full_partial_delta_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


load("CMP_PARTIALCORR/all_full_partial_delta_dt.Rdata")
load("CMP_PARTIALCORR/all_full_partial_corr_dt.Rdata")
load("ALLTADS_AND_PURITY/EPIC/log10/all_ds_corrPurity_dt.Rdata")

all_dt <- merge(all_full_partial_delta_dt, all_ds_corrPurity_dt, by =c("dataset", "region"), all=TRUE)
stopifnot(!is.na(all_dt))

  