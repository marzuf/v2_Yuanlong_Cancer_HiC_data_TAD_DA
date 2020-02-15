
require(foreach)
require(doMC)
registerDoMC(40)

# Rscript cmp_meanCorr_permut.R

plotType <- "png"
myHeight <- myWidth <- 400

outFolder <- "CMP_MEANCORR_PERMUT"
dir.create(outFolder, recursive = TRUE)
  
all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
hicds = all_hicds[1]
all_hicds <- all_hicds[grep("NCI-H460", all_hicds)]


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_ds_meanCorr <- foreach(hicds = all_hicds) %dopar% {
 get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata")))
}
names(all_ds_meanCorr) <- all_hicds

outFile <- file.path(outFolder, paste0("allDS_meanCorr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot_multiDens(
  all_ds_meanCorr, 
  plotTit = paste0("all hicds - n=", length(all_hicds), " - TAD meanCorr"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


all_ds_empPval <- foreach(hicds = all_hicds) %dopar% {
 empPval <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
 empPval <- p.adjust(empPval, method="BH")
 -log10(empPval)
}
names(all_ds_empPval) <- all_hicds

outFile <- file.path(outFolder, paste0("allDS_empPval_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot_multiDens(
  all_ds_empPval, 
  plotTit = paste0("all hicds - n=", length(all_hicds), " - adj. emp. p-val [-log10]"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))










