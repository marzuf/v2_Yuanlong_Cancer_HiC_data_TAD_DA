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