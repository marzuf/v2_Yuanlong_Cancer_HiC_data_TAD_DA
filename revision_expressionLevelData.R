require(foreach)
require(doMC)
registerDoMC(40)
require(ggsci)
require(ggpubr)
require(ggplot2)

# Rscript revision_expressionLevelData.R

setDir <- "/media/electron"
setDir <- ""

plotType <- "svg"
myHeightGG <- 5
myWidthGG <- 7

outFolder <- "REVISION_EXPRESSION_LEVELDATA"
dir.create(outFolder)

buildTable <- TRUE

aggFun <- "median"
log10offset <- 0.01

nHistBreaks1 <- 5
nHistBreaks2 <- 10

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$symbol))
symb2entrez <- setNames( gff_dt$entrezID, gff_dt$symbol)


pipFolder <- file.path("PIPELINE/OUTPUT_FOLDER/")
settingFolder <- file.path("PIPELINE/INPUT_FILES/")

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_dt <- get(load(final_table_file))
final_table_DT <- final_dt
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)

final_table_DT$dataset <- file.path(final_table_DT$hicds, final_table_DT$exprds)

all_ds <- unique(final_table_DT$dataset)

# all_ds = all_ds[7]

mainFolder <- "."

all_exprLevel_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  cat(paste0("> Start ", ds, "\n"))
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
    
  
  g2t_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
  g2t_dt <- read.delim(g2t_file, col.names=c("entrezID", "chromo", "start", "end", "region"), header=FALSE, stringsAsFactors = FALSE)
  g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
  stopifnot(!duplicated(g2t_dt$entrezID))
  stopifnot(g2t_dt$entrezID %in% names(entrez2symb))
  g2t_dt$symbol <- entrez2symb[paste0(g2t_dt$entrezID)]
  stopifnot(!is.na(g2t_dt$symbol))
  stopifnot(!duplicated(g2t_dt$symbol))
  symbol2tad <- setNames(g2t_dt$region, g2t_dt$symbol)
  # load the rnaseq data prep exprlevel
  # for each gene quantile of expr (-> something with hist ?)
  
  settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  
  samp1 <- get(load(file.path(setDir, sample1_file)))
  samp2 <- get(load(file.path(setDir, sample2_file)))
  
  geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
  stopifnot(geneList %in% gff_dt$entrezID)
  
  count_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")  # corrected here 11.03.2020 discussion to marco which data to plot
  stopifnot(file.exists(count_file))
  fpkm_dt <- get(load(count_file))
  # changed here 11.03.2020 -> should be normalized sample-wise
  fpkm_dt2 <- apply(fpkm_dt, 2, function(x)x/sum(x))
  # stopifnot(colSums(fpkm_dt2) == 1)
  stopifnot(abs(colSums(fpkm_dt2) - 1) <= 10^-4)
  # and then multiply by 10^6 to have FPKM
  fpkm_dt2 <- fpkm_dt2*10^6
  fpkm_dt2 <- data.frame(fpkm_dt2, check.names = FALSE)
  stopifnot(dim(fpkm_dt) == dim(fpkm_dt2))
  stopifnot(rownames(fpkm_dt) == rownames(fpkm_dt2))
  stopifnot(colnames(fpkm_dt) == colnames(fpkm_dt2))
  
  fpkm_dt <- fpkm_dt2
  stopifnot(names(geneList) %in% rownames(fpkm_dt))
  
  fpkm_plot_dt <- fpkm_dt[rownames(fpkm_dt) %in% names(geneList),]
  fpkm_plot_dt$entrezID <- geneList[paste0(rownames(fpkm_plot_dt))]
  stopifnot(!is.na(fpkm_plot_dt$entrez))
  stopifnot(fpkm_plot_dt$entrez %in% gff_dt$entrezID)
  fpkm_plot_dt$symbol <- entrez2symb[fpkm_plot_dt$entrezID]
  stopifnot(!is.na(fpkm_plot_dt$symbol))
  stopifnot(fpkm_plot_dt$symbol %in% gff_dt$symbol)
  
  stopifnot(samp1 %in% colnames(fpkm_plot_dt))
  stopifnot(samp2 %in% colnames(fpkm_plot_dt))
  
  fpkm_plot_dt_sub <- fpkm_plot_dt[,c(samp1, samp2)]
  fpkm_plot_dt_sub <- log10(fpkm_plot_dt_sub + log10offset)
  
  aggExpr <- apply(fpkm_plot_dt_sub, 1, FUN=aggFun)
  stopifnot( do.call(aggFun, list(as.numeric(fpkm_plot_dt_sub[1,]))) == aggExpr[1] )
  
  zscoreAggExpr <- as.numeric(scale(aggExpr))
  qqnormAggExpr <- qqnorm(aggExpr, plot.it=F)$x
  
  # first retrieve the breaks range - 1 for 5 breaks
  # all_aggExpr_hist <- hist(aggExpr, breaks=nHistBreaks1, plot = FALSE) ## !!! will not work, because if single value, passed to pretty
  all_aggExpr_hist <- hist(aggExpr, breaks = seq(min(aggExpr), max(aggExpr), length.out=nHistBreaks1+1), plot=FALSE)# $breaks
  
  # then for each value find in which break it falls
  histqt1AggExpr <- sapply(aggExpr, function(x) { 
    xbreak <- which(hist(x, breaks = all_aggExpr_hist$breaks, plot=FALSE)$counts == 1);
    stopifnot(length(xbreak) == 1); 
    xbreak})
  check_histqt <- factor(histqt1AggExpr, levels = seq_along(all_aggExpr_hist$counts))
  stopifnot(table(check_histqt) == all_aggExpr_hist$counts)
  stopifnot(max(histqt1AggExpr) <= nHistBreaks1)
  stopifnot(min(histqt1AggExpr) >= 1)
  
  # first retrieve the breaks range - 2 for 10 breaks
  # all_aggExpr_hist <- hist(aggExpr, breaks=nHistBreaks2, plot = FALSE)
  all_aggExpr_hist <- hist(aggExpr, breaks = seq(min(aggExpr), max(aggExpr), length.out=nHistBreaks2+1), plot=FALSE)# $breaks
  
  # then for each value find in which break it falls
  histqt2AggExpr <- sapply(aggExpr, function(x) { 
    xbreak <- which(hist(x, breaks = all_aggExpr_hist$breaks, plot=FALSE)$counts == 1);
    stopifnot(length(xbreak) == 1); 
    xbreak})
  check_histqt <- factor(histqt2AggExpr, levels = seq_along(all_aggExpr_hist$counts))
  stopifnot(table(check_histqt) == all_aggExpr_hist$counts)
  stopifnot(max(histqt2AggExpr) <= nHistBreaks2)
  stopifnot(min(histqt2AggExpr) >= 1)
  
  stopifnot(names(aggExpr) == rownames(fpkm_plot_dt))
  
  stopifnot(length(aggExpr) == nrow(qqnormAggExpr))
  stopifnot(length(aggExpr) == nrow(zscoreAggExpr))
  stopifnot(length(aggExpr) == nrow(histqt1AggExpr))
  stopifnot(length(aggExpr) == nrow(histqt2AggExpr))
  
  all_entrez <- symb2entrez[paste0(fpkm_plot_dt$symbol)]
  stopifnot(!is.na(all_entrez))
  
  out_dt <- data.frame(
    hicds=hicds, 
    exprds=exprds,
    symbol = fpkm_plot_dt$symbol,
    entrezID = all_entrez,
    region = as.character(symbol2tad[paste0(fpkm_plot_dt$symbol)]),
    aggLog10Expr = aggExpr,
    zscoreAggExpr=zscoreAggExpr,
    qqnormAggExpr=qqnormAggExpr,
    histqt1AggExpr=histqt1AggExpr,
    histqt2AggExpr = histqt2AggExpr,
    stringsAsFactors = FALSE
  )
  rownames(out_dt) <- NULL
  out_dt
}

outFile <- file.path(outFolder, "all_exprLevel_dt.Rdata")
save(all_exprLevel_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile,"\n"))

# since I have filtered by pipeline gene, should be true
stopifnot(!is.na(all_exprLevel_dt$region))

# x1 <- c(5,28,75,90)
# x1
# # [1]  5 28 75 90
# x2 <- c(50, 280, 750, 900)
# 
# # first retrieve the breaks range
# all_x1_hist <- hist(x1, breaks=10, plot = FALSE)
# x1_histBreaks <- sapply(x1, function(x) { 
#   xbreak <- which(hist(x, breaks = all_x1_hist$breaks, plot=FALSE)$counts == 1);
#     stopifnot(length(xbreak) == 1); 
#   xbreak})
# check_x1 <- factor(x1_histBreaks, levels = seq_along(all_x1_hist$counts))
# stopifnot(table(check_x1) == all_x1_hist$counts)
# 
# all_x2_hist <- hist(x2, breaks=10, plot = FALSE)
# x2_histBreaks <- sapply(x2, function(x) { xbreak <- which(hist(x, breaks = all_x2_hist$breaks, plot=FALSE)$counts == 1); stopifnot(length(xbreak) == 1); xbreak})
# 
# stopifnot(x1_histBreaks==x2_histBreaks)
# # 
# # > x1o=
# # > str(x10)
# Error in str(x10) : object 'x10' not found
# > str(x1o)
# List of 6
# $ breaks  : num [1:10] 0 10 20 30 40 50 60 70 80 90
# 
