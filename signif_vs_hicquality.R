result_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

qt_dt<- read.delim("hicds_sparsity.csv", sep="\t", header = TRUE)

signif_thresh <- 0.01

nSignif_dt <- aggregate(adjPvalComb  ~ hicds + exprds, function(x) sum(x <= signif_thresh), data = result_dt)

nSignif_qt_dt <- merge(nSignif_dt, qt_dt, by="hicds", all.x=TRUE)
stopifnot(!is.na(nSignif_qt_dt))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path("SIGNIF_VS_HICQUALITY", "signif_vs_quality.png")
dir.create(dirname(outFile))
png(outFile, height=400, width=400)
plot(
  adjPvalComb~mean, data = nSignif_qt_dt,
  xlab = "quality (mean all chr.)",
  ylab ="# signif. TADs",
  pch=16,
  cex=0.7
)
addCorr(x=nSignif_qt_dt$mean, y=nSignif_qt_dt$adjPvalComb, bty="n")
mtext(side=3, text = paste0("all DS - n=", length(unique(file.path(nSignif_qt_dt$hicds, nSignif_qt_dt$exprds)))))
dev.off()

# Rscript signif_vs_hicquality.R