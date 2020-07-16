options(scipen=100)

##
all_dt <- get(load("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2/conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
all_dt$nDS <- sapply(all_dt$corresp_tads, function(x) length(unlist(strsplit(x, ","))))
all_dt <- all_dt[order(all_dt$nDS, decreasing = T),]
View(all_dt)


# v2 version
# -> one single script for all and for cmpType
# -> if share more than 80% of the genes -> merge conserved regions
# -> count as conserved in one dataset at the exprds level
# -> min conserved region

SSHFS=F

# Rscript nSignif_nSignifPartial.R
# Rscript nSignif_nSignifPartial.R EPIC

startTime <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
} else{
  purity_ds <- ""
}
if(purity_ds == "") {
  purity_plot_name <- "aran"
  # all the ranks are between 1 and 0
} else if(purity_ds == "EPIC") {
  purity_plot_name <- "EPIC"
} else {
  stop("--invalid purity_ds\n")
}

purity_ds <- "CPE"
purity_plot_name <- "Aran - CPE"

transfExpr <- "log10"

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


outFolder <- file.path("NSIGNIF_NSIGNIFPARTIAL", purity_ds, transfExpr)
dir.create(outFolder, recursive = TRUE)

corrPurityQtThresh <- 0.05
adjPvalThresh <- 0.01


result_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
resultData$adjPvalComb_full <- resultData$adjPvalComb
resultData$signif_full <- resultData$adjPvalComb <= adjPvalThresh

partial_file <- file.path("CREATE_FINAL_TABLE_PARTIAL", "all_result_dt.Rdata")
partialData <- get(load(partial_file))
partialData$dataset <- file.path(partialData$hicds, partialData$exprds)
partialData$adjPvalComb_partial <- partialData$adjPvalComb
partialData$signif_partial <- partialData$adjPvalComb <= adjPvalThresh

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


purity_file <- file.path("ALLTADS_AND_PURITY", purity_ds, transfExpr, "all_ds_corrPurity_dt.Rdata")
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)


pur_merge_dt <- merge(merge_dt, agg_purity, by=c("dataset", "region"))
purityCorrThresh <- as.numeric(quantile(pur_merge_dt$purityCorr[!pur_merge_dt$signif_full], probs = corrPurityQtThresh ))
pur_merge_dt$purityFlagged <- pur_merge_dt$purityCorr <= purityCorrThresh


pur_merge_dt$fullSignifAndFlagged <- pur_merge_dt$signif_full & pur_merge_dt$purityFlagged

pur_merge_dt$partialSignifAndFlagged <- pur_merge_dt$signif_partial & pur_merge_dt$purityFlagged

nFullSignifFlagged_dt <- aggregate(fullSignifAndFlagged~dataset, data=pur_merge_dt, FUN=sum)
nPartialSignifFlagged_dt <- aggregate(partialSignifAndFlagged~dataset, data=pur_merge_dt, FUN=sum)

all_dt <- merge(nPartial_merge_dt, merge(nFulll_merge_dt, merge(nFullSignifFlagged_dt, nPartialSignifFlagged_dt, by="dataset"), by="dataset"), by="dataset")


pur_merge_dt$fullSignifAndFlagged_signifPartial <- pur_merge_dt$fullSignifAndFlagged & pur_merge_dt$signif_partial

nBothSignifAndFlagged_dt <- aggregate(fullSignifAndFlagged_signifPartial~dataset, FUN=sum, data=pur_merge_dt)

all_dt <- merge(merge(nFulll_merge_dt, nFullSignifFlagged_dt, by="dataset"), nBothSignifAndFlagged_dt, by="dataset")
all_dt$ratioFlagged <- all_dt$fullSignifAndFlagged/all_dt$signif_full
all_dt <- all_dt[order(all_dt$ratioFlagged, decreasing = TRUE),]
all_dt$ratioFlagged <- round(all_dt$ratioFlagged, 4)
all_dt$flagged_alsoPartialSignif <- all_dt$fullSignifAndFlagged_signifPartial/all_dt$fullSignifAndFlagged
all_dt$flagged_alsoPartialSignif <- round(all_dt$flagged_alsoPartialSignif,4)
outFile <- file.path(outFolder, "nsignif_flagged_full_partial_dt.txt")
write.table(all_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
cat(paste0("... written: ", outFile, "\n"))


# 
# merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))
# merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
# purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))
# merge_dt$purityFlagged <- merge_dt$purityCorr <= purityCorrThresh
# merge_dt$signifFlagged <- merge_dt$signif & merge_dt$purityFlagged
# 
# aggSignif_merge_dt <- aggregate(signif~dataset, FUN=sum, data=merge_dt)
# colnames(aggSignif_merge_dt)[2] <- "nSignif"
# stopifnot(sum(aggSignif_merge_dt$nSignif) ==sum(merge_dt$signif))
# aggFlagged_merge_dt <- aggregate(purityFlagged~dataset, FUN=sum, data=merge_dt)
# colnames(aggFlagged_merge_dt)[2] <- "nPurityFlagged"
# stopifnot(sum(aggFlagged_merge_dt$nPurityFlagged) ==sum(merge_dt$purityFlagged))
# aggSignifFlagged_merge_dt <- aggregate(signifFlagged~dataset, FUN=sum, data=merge_dt)
# colnames(aggSignifFlagged_merge_dt)[2] <- "nSignifAndFlagged"
# stopifnot(sum(aggSignifFlagged_merge_dt$nSignifAndFlagged) ==sum(merge_dt$signif & merge_dt$purityFlagged ))
# 
# all_dt <- merge(merge(aggSignif_merge_dt, aggFlagged_merge_dt, by="dataset", all=TRUE ),aggSignifFlagged_merge_dt,by="dataset", all=TRUE)
# all_dt$ratioSignifFlagged <- all_dt$nSignifAndFlagged/all_dt$nSignif
# stopifnot(all_dt$ratioSignifFlagged >= 0 & all_dt$ratioSignifFlagged <= 1)
# all_dt <- all_dt[order(all_dt$ratioSignifFlagged, decreasing = TRUE),]          
# 
# all_dt$ratioSignifFlagged <- round(all_dt$ratioSignifFlagged,4)
# 
# outFile <- file.path(outFolder, "all_dt_signif_flagged.txt")
# write.table(all_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
# cat(paste0("... written: ", outFile, "\n"))
# 
# resultData$regID <- file.path(resultData$hicds, resultData$exprds, resultData$region)
# flagged_signif_dt <- resultData[resultData$regID %in% merge_dt$regID[merge_dt$signifFlagged],c("regID", "region_genes")]
# flagged_signif_dt <- flagged_signif_dt[order(flagged_signif_dt$regID),]
# 
# flagged_signif_genes_dt <- do.call(rbind, apply(flagged_signif_dt, 1, function(x) data.frame(conserved_region=unique(x["regID"]), 
#                                                                                              symbol=unlist(strsplit(x["region_genes"], ",")),
#                                                                                              stringsAsFactors = FALSE)))
# rownames(flagged_signif_genes_dt) <- NULL
# 
# nFlagged_genes_dt <- data.frame(
#   symbol = names(table(flagged_signif_genes_dt$symbol)),
#   nSignifFlagged=as.numeric(table(flagged_signif_genes_dt$symbol)),
#   stringsAsFactors = FALSE
# )
# nFlagged_genes_dt <- nFlagged_genes_dt[order(nFlagged_genes_dt$nSignifFlagged, decreasing = TRUE),]
# 
# outFile <- file.path(outFolder, "nFlagged_genes.txt")
# write.table(nFlagged_genes_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
# cat(paste0("... written: ", outFile, "\n"))
# 








