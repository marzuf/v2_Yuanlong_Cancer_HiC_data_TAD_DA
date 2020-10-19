# Rscript prep_subb_table.R

outFolder <- "PREP_SUPP_TABLE"
dir.create(outFolder, recursive = TRUE)

dt1 <- get(load("PREP_TAD_SUMMARY_DT/all_ds_dt.Rdata"))
dt1$dataset <- file.path(dt1$hicds, dt1$exprds)
dt2 <- get(load("SIGNIF_PURITYFLAGGED_FINAL/aran/CPE/log10/all_dt.Rdata"))
out1 <- merge(dt1, dt2, by="dataset")
source("../MANUSCRIPT_FIGURES/full_dataset_names.R")
out1$hicds_lab <- hicds_names[out1$hicds]
out1$exprds_lab <- exprds_names[out1$exprds]

out1 <- out1[, c("hicds_lab","exprds_lab","nSamp1","nSamp2", "nTot", "nSignif","nPurityFlagged","nSignifAndFlagged")]
stopifnot(!is.na(out1))

outFile <- file.path(outFolder, "suppTableA.txt")
write.table(out1, file = outFile, sep="\t", quote=F, col.names = T, row.names=F)
cat(paste0("... outFile ", outFile, "\n"))


dt3 <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

dt4 <- get(load("ALLTADS_AND_PURITY_FINAL/aran/CPE/log10/all_ds_corrPurity_dt.Rdata"))


purity_ds <- "aran"
pm <- "CPE"
 transfExpr <- "log10"
 runFolder <- "."
 signifThresh <- 0.01
 corrPurityQtThresh <- 0.05
purity_file <- file.path(runFolder,"ALLTADS_AND_PURITY_FINAL", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"), all=FALSE)
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))
merge_dt$purityFlagged <- merge_dt$purityCorr <= purityCorrThresh
stopifnot(!is.na(merge_dt))


out2 <- merge(resultData, merge_dt[,c("dataset", "region", "purityCorr" )], by=c("dataset", "region"), all=TRUE)
stopifnot(nrow(out2) == nrow(resultData))
out2$hicds_lab <- hicds_names[out2$hicds]
out2$exprds_lab <- exprds_names[out2$exprds]
out2$purityFlagged <- out2$purityCorr <= purityCorrThresh
out2$signif <- out2$adjPvalComb <= signifThresh

out2 <- out2[,c("hicds_lab","exprds_lab","region","start","end","region_genes",
                "meanLogFC","meanCorr", "purityCorr")]

outFile <- file.path(outFolder, "suppTableC.txt")
write.table(out2, file = outFile, sep="\t", quote=F, col.names = T, row.names=F)
cat(paste0("... outFile ", outFile, "\n"))

 
dt5 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2_WITHTHRESH_ALLBARS/8/fig4B_topSignifGO_conservedRegions_dt_allBars.Rdata"))
outFile <- file.path(outFolder, "suppTableB.txt")
write.table(dt5, file = outFile, sep="\t", quote=F, col.names = T, row.names=F)
cat(paste0("... outFile ", outFile, "\n"))
