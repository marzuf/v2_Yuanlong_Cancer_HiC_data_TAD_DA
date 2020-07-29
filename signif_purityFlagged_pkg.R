options(scipen=100)

# v2 version
# -> one single script for all and for cmpType
# -> if share more than 80% of the genes -> merge conserved regions
# -> count as conserved in one dataset at the exprds level
# -> min conserved region

SSHFS=F

# Rscript signif_purityFlagged_pkg.R CPE
# Rscript signif_purityFlagged_pkg.R ESTIMATE


purity_ds <- "aran_pkg"
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  pm <- args[1]  
} else{
  stop("error\n")
}
purity_plot_name <- paste0("aran pkg - ", pm)

script0_name <- "0_prepGeneData"



### HARD-CODED - MAIN SETTINGS

corMet <- "pearson"
transfExpr <- "log10"
signifThresh <- 0.01
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifcol <- paste0(signif_column, "_", signifThresh)

# args <- commandArgs(trailingOnly = TRUE)
# if(length(args) == 1) {
#   purity_ds <- args[1]  
#   purity_plot_name <- "EPIC"
# } else{
#   purity_ds <- ""
#   purity_plot_name <- "aran"
# }

script_name <- "signif_purityFlagged.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

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


outFolder <- file.path("SIGNIF_PURITYFLAGGED_PKG", purity_plot_name, transfExpr)
dir.create(outFolder, recursive = TRUE)

purity_file <- file.path("ALLTADS_AND_PURITY_PKG", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata")
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

agg_purity$regID <- file.path(agg_purity$dataset, agg_purity$region)

result_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))
merge_dt$purityFlagged <- merge_dt$purityCorr <= purityCorrThresh
merge_dt$signifFlagged <- merge_dt$signif & merge_dt$purityFlagged

aggSignif_merge_dt <- aggregate(signif~dataset, FUN=sum, data=merge_dt)
colnames(aggSignif_merge_dt)[2] <- "nSignif"
stopifnot(sum(aggSignif_merge_dt$nSignif) ==sum(merge_dt$signif))
aggFlagged_merge_dt <- aggregate(purityFlagged~dataset, FUN=sum, data=merge_dt)
colnames(aggFlagged_merge_dt)[2] <- "nPurityFlagged"
stopifnot(sum(aggFlagged_merge_dt$nPurityFlagged) ==sum(merge_dt$purityFlagged))
aggSignifFlagged_merge_dt <- aggregate(signifFlagged~dataset, FUN=sum, data=merge_dt)
colnames(aggSignifFlagged_merge_dt)[2] <- "nSignifAndFlagged"
stopifnot(sum(aggSignifFlagged_merge_dt$nSignifAndFlagged) ==sum(merge_dt$signif & merge_dt$purityFlagged ))

all_dt <- merge(merge(aggSignif_merge_dt, aggFlagged_merge_dt, by="dataset", all=TRUE ),aggSignifFlagged_merge_dt,by="dataset", all=TRUE)
all_dt$ratioSignifFlagged <- all_dt$nSignifAndFlagged/all_dt$nSignif
stopifnot(all_dt$ratioSignifFlagged >= 0 & all_dt$ratioSignifFlagged <= 1)
all_dt <- all_dt[order(all_dt$ratioSignifFlagged, decreasing = TRUE),]          

all_dt$ratioSignifFlagged <- round(all_dt$ratioSignifFlagged,4)

outFile <- file.path(outFolder, "all_dt_signif_flagged.txt")
write.table(all_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
cat(paste0("... written: ", outFile, "\n"))

resultData$regID <- file.path(resultData$hicds, resultData$exprds, resultData$region)
flagged_signif_dt <- resultData[resultData$regID %in% merge_dt$regID[merge_dt$signifFlagged],c("regID", "region_genes")]
flagged_signif_dt <- flagged_signif_dt[order(flagged_signif_dt$regID),]

flagged_signif_genes_dt <- do.call(rbind, apply(flagged_signif_dt, 1, function(x) data.frame(conserved_region=unique(x["regID"]), 
                                                                                                                           symbol=unlist(strsplit(x["region_genes"], ",")),
                                                                                                                           stringsAsFactors = FALSE)))
rownames(flagged_signif_genes_dt) <- NULL

nFlagged_genes_dt <- data.frame(
  symbol = names(table(flagged_signif_genes_dt$symbol)),
  nSignifFlagged=as.numeric(table(flagged_signif_genes_dt$symbol)),
  stringsAsFactors = FALSE
)
nFlagged_genes_dt <- nFlagged_genes_dt[order(nFlagged_genes_dt$nSignifFlagged, decreasing = TRUE),]

outFile <- file.path(outFolder, "nFlagged_genes.txt")
write.table(nFlagged_genes_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
cat(paste0("... written: ", outFile, "\n"))












