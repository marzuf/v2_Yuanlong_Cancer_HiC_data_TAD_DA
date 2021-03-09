# Rscript revision_lookGOgenes_withPurity.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- "REVISION_LOOKGOGENES_WITH_PURITY"
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.2
plotType <- "png"
myHeight <- myWidth <-400

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

goSignif_thresh <- 0.05
#### #### #### #### #### #### #### #### #### #### #### #### 
#### retrieve the purity tagged tads

runFolder <- "."

purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- paste0("Aran - ", pm)

corMet <- "pearson"
transfExpr <- "log10"
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)


purity_file <- file.path(runFolder,"ALLTADS_AND_PURITY_FINAL", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)
agg_purity$region_ID <- file.path(agg_purity$dataset, agg_purity$region)
stopifnot(!duplicated(agg_purity$region_ID))

result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt$region_ID <- file.path(merge_dt$dataset, merge_dt$region)
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

# purity_flagged_tads <- merge_dt$region_ID[merge_dt$purityCorr <= purityCorrThresh]
tokeep_tads <- merge_dt$region_ID[merge_dt$purityCorr > purityCorrThresh]
cat(paste0( "purityCorrThresh = ", round(purityCorrThresh, 4), "\n"))


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 


all_infiles <-  "REVISION_GO_BYCPTMT/B22_keepnotPF_entrez_signif_enrich.Rdata"
infile <- "REVISION_GO_BYCPTMT/B22_keepnotPF_entrez_signif_enrich.Rdata"

# all_infiles <- c("REVISION_GO_SIGNIF/keepnotPF_entrez_signif_enrich.Rdata",
#                  "REVISION_GO_SIGNIF/discardPF_entrez_signif_enrich.Rdata",
#                  "REVISION_GO_BYCPTMT/B22_discardPF_entrez_signif_enrich.Rdata",
#                  "REVISION_GO_BYCPTMT/B22_keepnotPF_entrez_signif_enrich.Rdata",
#                  "REVISION_GO_BYCPTMT_WITHOUTPURITYFILTER/B22_withoutPF_entrez_signif_enrich.Rdata")

all_infiles <- infile <- "REVISION_GO_BYCPTMT_WITHOUTPURITYFILTER/B22_withoutPF_entrez_signif_enrich.Rdata"

         
for(infile in all_infiles) {
  
  
  if(infile == "REVISION_GO_BYCPTMT/B22_keepnotPF_entrez_signif_enrich.Rdata") {
    datatype <- "withPF_notPF"
    
    purityData$region_ID <- file.path(purityData$dataset, purityData$region)
    stopifnot(tokeep_tads %in% purityData$region_ID)
    filtered_pd_dt <- purityData[purityData$region_ID %in% tokeep_tads,]
   
    
  } else  if(infile == "REVISION_GO_BYCPTMT_WITHOUTPURITYFILTER/B22_withoutPF_entrez_signif_enrich.Rdata") {
     datatype <- "withoutPF"
     filtered_pd_dt <- purityData
  } else {
    stop("-ok\n")
  }
  stopifnot(nrow(filtered_pd_dt)>0)
  
  stopifnot(filtered_pd_dt$entrezID %in% names(entrez2symb))
  filtered_pd_dt$symbol <- entrez2symb[paste0(filtered_pd_dt$entrezID)]
  stopifnot(!is.na(filtered_pd_dt$symbol))
  gene_agg_dt <- aggregate(purityCorr~symbol, FUN=mean, data=filtered_pd_dt)
  geneAggPurity_values <- setNames(gene_agg_dt$purityCorr, gene_agg_dt$symbol)
  
  
  # outfile <- gsub(".Rdata", "_topTable.txt", infile)
  
  dt=get(load(infile))
  result_dt <- dt@result
  
  nTop <- 10
  sort_col <- "p.adjust"
  
  result_dt <- result_dt[order(result_dt[,paste0(sort_col)]),]
  result_dt$go_rank <- 1:nrow(result_dt)
  
  # x=result_dt$geneID[1]
  result_dt$geneSymbol <- sapply(result_dt$geneID, function(x) {
    all_entrez <- unlist(strsplit(x, split="/"))
    stopifnot(all_entrez %in% names(entrez2symb))
    all_symbols <- entrez2symb[paste0(all_entrez)]
    paste0(all_symbols, collapse="/")
  })
  
  # x=result_dt$geneID[1]
  result_dt$meanPurityScore <- sapply(result_dt$geneID, function(x) {
    all_entrez <- unlist(strsplit(x, split="/"))
    stopifnot(all_entrez %in% names(entrez2symb))
    all_symbols <- entrez2symb[paste0(all_entrez)]
    
    if(datatype == "withPF_notPF") {
      stopifnot(all_symbols %in% names(geneAggPurity_values))
      mp <- mean(as.numeric(geneAggPurity_values[paste0(all_symbols)]))
      stopifnot(!is.na(mp))
    } else if(datatype == "withoutPF") {
      mp <- mean(as.numeric(geneAggPurity_values[paste0(all_symbols)]), na.rm=TRUE)
    } else {
      stop("error\n")
    }
    mp
  })
  
  plotTit <- paste0("GO signif. and averaged gene purity scores")
  
  for(toplot in c("all", "signif")) {
    
    if(toplot=="all"){
      plot_dt <- result_dt
      mySub <- paste0("all GOs - n=", nrow(plot_dt), " - ", datatype)
    } else if(toplot=="signif"){
        plot_dt <- result_dt[result_dt$p.adjust <= goSignif_thresh,]
        mySub <- paste0("signif GOs - n=", nrow(plot_dt), " (p.adjust <= ",goSignif_thresh, ")", " - ", datatype )
    } else {
        stop("--error\n")
      }
    
    myy <- -log10(plot_dt$p.adjust)
    myx <- plot_dt$meanPurityScore
    
    outFile <- file.path(outFolder,paste0("meanPurity_byGO_", toplot, "GOs_", datatype, "_densplot.", plotType))
    do.call(plotType, list(outFile, height=myWidth, width=myWidth))
    
    densplot(x=myx,y=myy,
             xlab="meanPurityScore",
             ylab="p.adjust [-log10]",
             cex.axis=plotCex,
             cex.main = plotCex,
             cex.lab = plotCex,
             main=plotTit)
    mtext(side=3, text=mySub)
    addCorr(x=myx,y=myy, legPos="topright", bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    myy <- plot_dt$go_rank
    myx <- plot_dt$meanPurityScore
    
    outFile <- file.path(outFolder,paste0("meanPurity_byGO_", toplot, "GOs_", datatype, "_ranks_densplot.", plotType))
    do.call(plotType, list(outFile, height=myWidth, width=myWidth))
    
    densplot(x=myx,y=myy,
             xlab="meanPurityScore",
             ylab="GO rank (sorted by p.adjust)",
             cex.axis=plotCex,
             cex.main = plotCex,
             cex.lab = plotCex,
             main=plotTit)
    mtext(side=3, text=mySub)
    addCorr(x=myx,y=myy, legPos="topright", bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
  
  # crop_dt <- result_dt[1:nTop,c("Description", sort_col, "geneID")]
  # out_dt <- crop_dt
  # out_dt$geneID <- NULL
  # out_dt[,paste0(sort_col)] <- round(out_dt[,paste0(sort_col)] ,4) 
  # write.table(out_dt, file = outfile, col.names=F, row.names=F, sep="\t", quote=F)
  # cat(paste0("... written ", outfile, "\n"))
  
}
