require(foreach)
require(doMC)
registerDoMC(40)

purity_dt <- get(load("ALLTADS_AND_PURITY_FINAL/aran/CPE/log10/all_ds_corrPurity_dt.Rdata"))

agg_purity_dt <- aggregate(purityCorr~dataset+region, data =purity_dt, FUN=mean)

pvals_dt <- foreach(ds = unique(purity_dt$dataset), .combine='rbind')%dopar% {
  
  
  td_pvals <- get(load(file.path(
    # setDir,
    # "/mnt/etemp/marie", td_folder,
    "PIPELINE/OUTPUT_FOLDER", dirname(ds), basename(ds), "11sameNbr_runEmpPvalCombined/emp_pval_combined.Rdata")))
  
  td_pvals <- p.adjust(td_pvals, method="BH")
  
  data.frame(
    dataset=ds,
    region = names(td_pvals),
    TAD_adjPval = as.numeric(td_pvals),
    stringsAsFactors = FALSE
  )
  
  
}


merge_dt <- merge(agg_purity_dt,pvals_dt, by=c("dataset", "region"), all.x=T, all.y=F)
stopifnot(!is.na(merge_dt))

tad_pval_and_purity_dt <- merge_dt

outFolder <- "PREP_TAD_SUMMARY_DT"
dir.create(outFolder)
save(tad_pval_and_purity_dt,file= file.path(outFolder, "tad_pval_and_purity_dt.Rdata"), version=2)


require(foreach)
require(doMC)
registerDoMC(40)

setDir <- ""

# Rscript prep_tad_summary_dt.R

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

settingFolder <- file.path("PIPELINE/INPUT_FILES/")


all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

outFolder <- file.path("PREP_TAD_SUMMARY_DT")
dir.create(outFolder, recursive = TRUE)

all_ds_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  cat(paste0("... start: ", ds, "\n"))
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
  
  settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  
  samp1 <- get(load(file.path(setDir, sample1_file)))
  samp2 <- get(load(file.path(setDir, sample2_file)))
  
  
  data.frame(
    hicds=hicds,
    exprds=exprds,
    nSamp1 = length(samp1),
    nSamp2 = length(samp2),
    stringsAsFactors = FALSE
  )
  
}
save(all_ds_dt, file=file.path(outFolder, "all_ds_dt.Rdata"), version = 2)

