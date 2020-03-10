
require(foreach)
require(doMC)
registerDoMC(40)

# Rscript size_nbr_vs_hicquality.R

outFolder <- "SIZE_NBR_VS_HICQUALITY"
dir.create(outFolder)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <-  "svg"
myHeight <- myWidth <- 7
plotCex <- 1.4

qt_dt<- read.delim("hicds_sparsity.csv", sep="\t", header = TRUE)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
hicds = all_hicds[1]

all_sizes_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  tad_dt <- read.delim(file.path(hicds, "genes2tad", "all_assigned_regions.txt"), header=FALSE, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
  stopifnot(!duplicated(tad_dt$region))
  meanSize <- mean(tad_dt$end - tad_dt$start + 1)
  data.frame(
    hicds = hicds,
    meanSize = meanSize,
    nbrTADs = nrow(tad_dt),
    stringsAsFactors = FALSE
  )
}

merge_dt <- merge(all_sizes_dt, qt_dt, by=c("hicds"), all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(merge_dt))


outFile <- file.path(outFolder, paste0("meanSize_vs_quality.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = merge_dt$mean,
  y = log10(merge_dt$meanSize),
  ylab = "mean TAD size [log10]",
  xlab = "Hi-C data quality",
  pch=16,
  cex=0.7,
  cex.lab = plotCex,
  cex.main = plotCex,
  cex.axis= plotCex
)
addCorr(x=merge_dt$mean, y=log10(merge_dt$meanSize), bty="n")
mtext(side=3, text = paste0("all DS - n=", length(unique(file.path(merge_dt$hicds)))))
dev.off()


outFile <- file.path(outFolder, paste0("tadNbr_vs_quality.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = merge_dt$mean,
  y = merge_dt$nbrTADs,
  ylab = "# TADs",
  xlab = "Hi-C data quality",
  pch=16,
  cex=0.7,
  cex.lab = plotCex,
  cex.main = plotCex,
  cex.axis= plotCex
)
addCorr(x=merge_dt$mean, y=log10(merge_dt$meanSize), bty="n")
mtext(side=3, text = paste0("all DS - n=", length(unique(file.path(merge_dt$hicds)))))
dev.off()




