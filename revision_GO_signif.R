# GO ANALYSIS ON THE GENES FROM SIGNIF TADs


# Rscript revision_GO_signif.R

printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}


script_name <- "purityFilter_density.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(clusterProfiler)
require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")


buildTable <- TRUE

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 5)
myWidth <- 8
axisCex <- 1.4

### HARD-CODED - MAIN SETTINGS
# _final:  discussion 04.08.2020 Giovanni - keep aran CPE data, first vial only

purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- paste0("Aran - ", pm)

corMet <- "pearson"
transfExpr <- "log10"
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)


source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 5)
myWidth <- 8

script0_name <- "0_prepGeneData"

####################################################################### FOR THE GO
source("enricher_settings.R")



padjVarGO <- "p.adjust" # p.adjust or qvalue ???


barplot_vars <- c("foldEnrichment", "geneRatio", "log10_pval")
barplot_vars_tit <- setNames(c("Fold enrichment", "Gene ratio", paste0("-log10(", padjVarGO,  ")")), barplot_vars)
barplot_vars <- c( "log10_pval")
barplot_vars_tit <- setNames(c(paste0("-log10(", padjVarGO,  ")")), barplot_vars)

plotMaxBars <- 10

par_bot <- 8

padjVarGO_plotThresh <- 0.05

# GO for BP nad MF [do not take c5_CC]
if(enricher_ontologyType == "BP" | enricher_ontologyType == "MF" | enricher_ontologyType == "BP_MF"){
  gmtFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2", paste0("c5.", tolower(enricher_ontologyType), ".v6.1.entrez.gmt"))
} else {
  stop(paste0(enricher_ontologyType, " is not a valid ontologyType\n"))
}
stopifnot(file.exists(gmtFile))
c5_msigdb <- read.gmt(gmtFile)


txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... gmtFile:\t", gmtFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher background genes:\t", "universe", "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pvalueCutoff:\t", enricher_pvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pAdjustMethod:\t", enricher_pAdjustMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher minGSSize:\t", enricher_minGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher maxGSSize:\t", enricher_maxGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher qvalueCutoff:\t", enricher_qvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher ontologyType:\t", enricher_ontologyType, "\n")
printAndLog(txt, logFile)
# txt <- paste0("... pval thresh select GOs:\t", pvalSelectGO, "\n")
# printAndLog(txt, logFile)
txt <- paste0("... enricher_results_sortGOby:\t", enricher_results_sortGOby, "\n")
printAndLog(txt, logFile)

txt <- paste0("... padjVarGO_plotThresh:\t", padjVarGO_plotThresh, "\n")
printAndLog(txt, logFile)
txt <- paste0("... plotMaxBars:\t", plotMaxBars, "\n")
printAndLog(txt, logFile)




outFolder <- file.path("REVISION_GO_SIGNIF")
dir.create(outFolder, recursive = TRUE)


setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


purity_file <- file.path(runFolder,"ALLTADS_AND_PURITY_FINAL", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt$region_id <- file.path(merge_dt$dataset, merge_dt$region)
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

purity_flagged_tads <- merge_dt$region_id[merge_dt$purityCorr <= purityCorrThresh]

# take those that have puritycorr annot + not flagged
tads_signif_notPurityFlagged <- merge_dt$region_id[ (! merge_dt$region_id %in% purity_flagged_tads) &
                                                             merge_dt$adjPvalComb <= signifThresh]
stopifnot(!purity_flagged_tads %in% tads_signif_notPurityFlagged)

# discard the flagged -> means that I keep those for which no puritycorr annot
merge_dt_all <- merge(agg_purity, resultData, by=c("dataset", "region"), all=TRUE)  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt_all$region_id <- file.path(merge_dt_all$dataset, merge_dt_all$region)
stopifnot(purity_flagged_tads %in%  merge_dt_all$region_id)
tads_signif_notPurityFlaggedAll <- merge_dt_all$region_id[ (! merge_dt_all$region_id %in% purity_flagged_tads ) &
                                                                merge_dt_all$adjPvalComb <= signifThresh]

stopifnot(tads_signif_notPurityFlagged %in% tads_signif_notPurityFlaggedAll)
stopifnot(! purity_flagged_tads  %in% tads_signif_notPurityFlaggedAll)
stopifnot(! purity_flagged_tads %in% tads_signif_notPurityFlagged)


tad_list  = tads_signif_notPurityFlagged
data_dt = merge_dt
stopifnot(tad_list %in% merge_dt$region_id)
signif_notPF_genes <- unique(unlist(strsplit(as.character(data_dt$region_genes[data_dt$region_id %in% tad_list]), split=",")))
stopifnot(signif_notPF_genes %in% gff_dt$symbol)
signif_notPF_entrez <- gff_dt$entrezID[gff_dt$symbol %in% signif_notPF_genes]

univers_genes <- unique(unlist(strsplit(as.character(data_dt$region_genes), split=",")))
stopifnot(univers_genes %in% gff_dt$symbol)
univers_entrez <- gff_dt$entrezID[gff_dt$symbol %in% univers_genes]




go_univers_entrez <- as.character(univers_entrez)[as.character(univers_entrez) %in% as.character(c5_msigdb$gene)]
stopifnot(length(go_univers_entrez)  > 0)
go_signif_notPF_entrez <- as.character(signif_notPF_entrez)[as.character(signif_notPF_entrez) %in% as.character(c5_msigdb$gene)]
stopifnot(length(go_signif_notPF_entrez)  > 0)


cat("... available annot. for univers_entrez:\t", length(go_univers_entrez) , "/", length(univers_entrez), "\n")
cat("... available annot. for univers_entrez:\t", length(go_signif_notPF_entrez) , "/", length(signif_notPF_entrez), "\n")


stopifnot(go_signif_notPF_entrez %in% go_univers_entrez)


  
  
  #***** 1) conserved_signif
  cat(paste0(">  start enricher for conserved_signif \n"))
  
  if(length(go_all_conserved_signif_tads_genes) > 0) {
    
    conserved_signif_enrich <- enricher(gene = go_all_conserved_signif_tads_genes, 
                                        TERM2GENE=c5_msigdb,
                                        universe = go_all_universe_genes,
                                        pvalueCutoff = enricher_pvalueCutoff, 
                                        pAdjustMethod = enricher_pAdjustMethod, 
                                        minGSSize = enricher_minGSSize, 
                                        maxGSSize = enricher_maxGSSize, 
                                        qvalueCutoff =enricher_qvalueCutoff)
    
    save(conserved_signif_enrich, version=2, file=file.path(outFolder, "conserved_signif_enrich.Rdata"))
    
    conserved_signif_enrich_resultDT <- conserved_signif_enrich@result
    conserved_signif_enrich_resultDT <- conserved_signif_enrich_resultDT[order(conserved_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
    conserved_signif_enrich_resultDT$log10_pval <- -log10(conserved_signif_enrich_resultDT[,paste0(padjVarGO)])
    conserved_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(conserved_signif_enrich_resultDT$GeneRatio, function(x) {
      gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    conserved_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(conserved_signif_enrich_resultDT$BgRatio, function(x) {
      gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    conserved_signif_enrich_resultDT$foldEnrichment <- conserved_signif_enrich_resultDT$geneRatio/conserved_signif_enrich_resultDT$bgRatio
    
    
    save(conserved_signif_enrich_resultDT, version=2, file=file.path(outFolder, "conserved_signif_enrich_resultDT_0.Rdata"))
    
    
    genes_signif_plotMax <- min(c(plotMaxBars, nrow(conserved_signif_enrich_resultDT)))
    if(genes_signif_plotMax > 0) {
      for(var_plot in barplot_vars) {
        myTit <- paste0(barplot_vars_tit[var_plot], " ", data_cmpType, " (conserved_signif)")
        conserved_signif_enrich_resultDT <- conserved_signif_enrich_resultDT[order(conserved_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
        
        my_x <- gsub("GO_", "", conserved_signif_enrich_resultDT$Description[1:genes_signif_plotMax])
        
        
        conserved_signif_dt <- conserved_signif_enrich_resultDT[1:genes_signif_plotMax,c(var_plot, "Description"), drop=FALSE]
        conserved_signif_dt$labs <-  gsub("GO_", "", conserved_signif_dt$Description)
        conserved_signif_dt$plot_labs <-  unlist(lapply(strwrap(gsub("_", " ", conserved_signif_dt$labs), 
                                                                width = strwdth, simplify=FALSE), 
                                                        function(x) paste0(x, collapse="\n")))
        
        outFile <- file.path(outFolder,paste0(file_prefix, "conserved_signif", var_plot, "_GO_",
                                              enricher_ontologyType, "_", var_plot, "_barplot_dt.Rdata"))
        save(conserved_signif_dt, file=outFile, version=2)
        
        
        outFile <- file.path(outFolder,paste0(file_prefix, "conserved_signif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
        do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.5))    
        par(oma=c(par_bot,1,1,1))
        
        
        
        barplot(conserved_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                main = myTit, 
                # ylab = paste0("-log10 (", padjVarGO, ")"),
                ylab = paste0(barplot_vars_tit[var_plot]),
                names.arg =unlist(lapply(strwrap(gsub("_", " ", my_x), width = strwdth, 
                                                 simplify=FALSE), function(x) paste0(x, collapse="\n"))), 
                las=2,
                cex.lab = plotCex,
                cex.main = plotCex,
                cex.axis=plotCex,
                cex.names = 0.9
        )
        mtext(side=3, text=paste0("(# datasets = ", nDS, ")"), font=3)
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
      }
    }
    conserved_signif_enrich_resultDT <- data.frame(conserved_signif_enrich_resultDT)
    
  } else {
    conserved_signif_enrich_resultDT <- NULL
  }
  
  
  
  

