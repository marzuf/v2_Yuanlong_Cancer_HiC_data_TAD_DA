options(scipen=100)

# v2 version
# -> one single script for all and for cmpType
# -> if share more than 80% of the genes -> merge conserved regions
# -> count as conserved in one dataset at the exprds level
# -> min conserved region

SSHFS=F

# Rscript nSignif_nSignifPartial.R


startTime <- Sys.time()


require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("plot_lolliTAD_funct.R")
# source("my_heatmap.2.R")

buildTable <- TRUE

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightHeat <- myHeight * 1.8
myWidthHeat <- myWidth * 1.8

myWidthGG <- 12
myHeightGG <- 12


outFolder <- file.path("NSIGNIF_NSIGNIFPARTIAL")
dir.create(outFolder, recursive = TRUE)


result_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
resultData$adjPvalComb_full <- resultData$adjPvalComb
resultData$signif_full <- resultData$adjPvalComb <= 0.01

partial_file <- file.path("CREATE_FINAL_TABLE_PARTIAL", "all_result_dt.Rdata")
partialData <- get(load(partial_file))
partialData$dataset <- file.path(partialData$hicds, partialData$exprds)
partialData$adjPvalComb_partial <- partialData$adjPvalComb
partialData$signif_partial <- partialData$adjPvalComb <= 0.01

merge_dt <- merge(resultData[,c("dataset", "region", "signif_full")],
                  partialData[,c("dataset", "region", "signif_partial")], by=c("dataset", "region"))
merge_dt$signif_both <- merge_dt$signif_full&merge_dt$signif_partial

nFulll_merge_dt <- aggregate(signif_full~dataset, FUN=sum, data=merge_dt)
nPartial_merge_dt <- aggregate(signif_partial~dataset, FUN=sum, data=merge_dt)
nBoth_merge_dt<- aggregate(signif_both~dataset, FUN=sum, data=merge_dt)

sum_agg_dt <- merge(merge(nBoth_merge_dt, nFulll_merge_dt, by="dataset"),nPartial_merge_dt, by="dataset")
sum_agg_dt <- sum_agg_dt[order(sum_agg_dt$signif_both, decreasing = TRUE),]
  
outFile <- file.path(outFolder, "nsignif_full_partial_dt.txt")
write.table(sum_agg_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
cat(paste0("... written: ", outFile, "\n"))



