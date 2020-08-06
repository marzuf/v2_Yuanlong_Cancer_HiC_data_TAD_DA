# Rscript cmp_purity_init_final.R

plotCex <- 1.2

outFolder <- "CMP_PURITY_INIT_FINAL"
dir.create(outFolder, recursive = TRUE)

v0 <- get(load(file.path("../OLDER_v2_Yuanlong_Cancer_HiC_data_TAD_DA/ALLTADS_AND_PURITY/CPE/log10//all_ds_corrPurity_dt.Rdata")))
v1 <- get(load(file.path("ALLTADS_AND_PURITY_FINAL/aran/CPE/log10//all_ds_corrPurity_dt.Rdata")))

v0_samp <- v0[,c("dataset", "nSampWithPurity")]
v0_samp <- unique(v0_samp)
v1_samp <- v1[,c("dataset", "nSampWithPurity")]
v1_samp <- unique(v1_samp)

merged_samp_dt <- merge(v0_samp, v1_samp, by="dataset", suffixes=c("_v0", "_v1"))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path(outFolder, "cmp_nSamp_purity_corrected_data.svg")
svg(outFile, height=6, width=6)
plot(merged_samp_dt$nSampWithPurity_v0,merged_samp_dt$nSampWithPurity_v1,
     main = paste0("aran CPE data - # av. samples"),
     pch=16, cex=0.7,
     cex.main=plotCex,cex.axis=plotCex,cex.lab=plotCex,
     xlab="v0 data", ylab="corrected data")
curve(1*x, add=TRUE, col="darkgrey")
addCorr(x=merged_samp_dt$nSampWithPurity_v0,
        y= merged_samp_dt$nSampWithPurity_v1,
        legPos = "topleft",
        bty="n")
foo <- dev.off()


agg_dt_0 <- aggregate(purityCorr~dataset+region, data=v0[,c("dataset", "region", "purityCorr")], FUN=mean)
agg_dt_1 <- aggregate(purityCorr~dataset+region, data=v1[,c("dataset", "region", "purityCorr")], FUN=mean)

agg_merged_values_dt <- merge(agg_dt_0[,c("dataset", "region", "purityCorr")],
                              agg_dt_1[,c("dataset", "region", "purityCorr")],
                          by=c("dataset", "region"), suffixes=c("_v0", "_v1"))

outFile <- file.path(outFolder, "cmp_tadLevelCorr_purity_corrected_data.png")
png(outFile, height=400, width=400)
plot(agg_merged_values_dt$purityCorr_v0,agg_merged_values_dt$purityCorr_v1,
     main = paste0("aran CPE data - purity corr. TAD-level"),
     pch=16, cex=0.7,
     cex.main=plotCex,cex.axis=plotCex,cex.lab=plotCex,
     xlab="v0 data", ylab="corrected data")
curve(1*x, add=TRUE, col="darkgrey")
addCorr(x=agg_merged_values_dt$purityCorr_v0,
        y=agg_merged_values_dt$purityCorr_v1,
        legPos = "topleft",
        bty="n")
foo <- dev.off()




merged_values_dt <- merge(v0[,c("dataset", "region", "entrezID", "purityCorr")],
                          v1[,c("dataset", "region", "entrezID", "purityCorr")],
                          by=c("dataset", "region", "entrezID"), suffixes=c("_v0", "_v1"))

outFile <- file.path(outFolder, "cmp_geneLevelCorr_purity_corrected_data.png")
png(outFile, height=400, width=400)
plot(merged_values_dt$purityCorr_v0,merged_values_dt$purityCorr_v1,
     main = paste0("aran CPE data - purity corr. gene-level"),
     pch=16, cex=0.7,
     cex.main=plotCex,cex.axis=plotCex,cex.lab=plotCex,
     xlab="v0 data", ylab="corrected data")
curve(1*x, add=TRUE, col="darkgrey")
addCorr(x=merged_values_dt$purityCorr_v0,
        y=merged_values_dt$purityCorr_v1,
        legPos = "topleft",
        bty="n")
foo <- dev.off()





cat(paste0("# data v0\t", sum(v0$nSampWithPurity), "\n"))
cat(paste0("# data v1\t", sum(v1$nSampWithPurity), "\n"))