options(scipen=100)

SSHFS=F


# Rscript go_signif_geneLevel_tadLevel.R

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "go_signif_geneLevel_tadLevel.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(clusterProfiler)
require(foreach)
require(doMC)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("plot_lolliTAD_funct.R")
source("my_heatmap.2.R")


buildTable <- TRUE

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightHeat <- myHeight * 1.8
myWidthHeat <- myWidth * 1.8

myWidthGG <- 12
myHeightGG <- 12

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("enricher_settings.R")
# from enricher_settings, load:
# enricher_ontologyType <- "BP"
# enricher_pvalueCutoff <- 1
# enricher_pAdjustMethod <- "BH"
# enricher_minGSSize <- 1
# enricher_maxGSSize <- 500
# enricher_qvalueCutoff <- 1
# enricher_results_sortGOby <- "p.adjust"


script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- file.path("GO_SIGNIF_GENELEVEL_TADLEVEL")
dir.create(outFolder, recursive = TRUE)

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))


final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))
signif_column <- "adjPvalComb"
tads_signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", tads_signifThresh)
final_dt[, paste0(signifcol)] <- final_dt[, paste0(signif_column)] <= tads_signifThresh
cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> tads_signifThresh\t=\t", tads_signifThresh, "\n"))


genes_signifTresh <- 0.05

padjVarGO <- "p.adjust" # p.adjust or qvalue ???

logFile <- file.path(outFolder, "go_signif_geneLevel_tadLevel_logFile.txt")
if(buildTable) file.remove(logFile)

barplot_vars <- c("foldEnrichment", "geneRatio", "log10_pval")
barplot_vars_tit <- setNames(c("Fold enrichment", "Gene ratio", paste0("-log10(", padjVarGO,  ")")), barplot_vars)
plotMaxBars <- 10
padjVarGO_plotThresh <- 0.05

# GO for BP nad MF [do not take c5_CC]
if(enricher_ontologyType == "BP" | enricher_ontologyType == "MF" | enricher_ontologyType == "BP_MF"){
  gmtFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2", paste0("c5.", tolower(enricher_ontologyType), ".v6.1.entrez.gmt"))
} else {
  stop(paste0(enricher_ontologyType, " is not a valid ontologyType\n"))
}
stopifnot(file.exists(gmtFile))
c5_msigdb <- read.gmt(gmtFile)


printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}


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


####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start preparing data before matching \n")
  
  hicds = all_hicds[1]
  all_go_enrich_list <- foreach(hicds = all_hicds) %dopar% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      
      
      
      geneList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      pipeline_geneList <- get(load(geneList_file))
      
      stopifnot(pipeline_geneList %in% g2t_dt$entrezID)
      
      exprds_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% pipeline_geneList,]
      
      ### retrieve genes from signif. tads
      signif_tads <- final_dt$region[final_dt$hicds == hicds &
                                       final_dt$exprds == exprds &
                                       final_dt[,paste0(signifcol)]]
      
      
      tads_signif_genes <- exprds_g2t_dt$entrezID[exprds_g2t_dt$region %in% signif_tads]
      stopifnot(tads_signif_genes %in% pipeline_geneList)
      
      ### retrieve limma signif genes
      topTable_DT_file <- file.path(pipFolder,  hicds, exprds, script1_name, "DE_topTable.Rdata")
      stopifnot(file.exists(topTable_DT_file))
      topTable_DT <- get(load(topTable_DT_file))
      topTable_DT$genes <- as.character(topTable_DT$genes)
      stopifnot(names(pipeline_geneList) %in% topTable_DT$genes)
      topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipeline_geneList),]
      topTable_DT$entrezID <- pipeline_geneList[paste0(topTable_DT$genes)]
      stopifnot(!is.na(topTable_DT$entrezID))
      stopifnot(!duplicated(topTable_DT$entrezID))
      stopifnot(topTable_DT$entrezID %in% pipeline_geneList)
      topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipeline_geneList),]
      topTable_DT$entrezID <- pipeline_geneList[paste0(topTable_DT$genes)]
      stopifnot(!is.na(topTable_DT$entrezID))
      stopifnot(!duplicated(topTable_DT$entrezID))
      
      limma_signif_genes <- topTable_DT$entrezID[topTable_DT$adj.P.Val <= genes_signifTresh]
      stopifnot(limma_signif_genes %in% pipeline_geneList)
      
    
      universe_genes <- pipeline_geneList
      
      go_universe_genes <- as.character(universe_genes)[as.character(universe_genes) %in% as.character(c5_msigdb$gene)]
      go_tads_signif_genes <- as.character(tads_signif_genes)[as.character(tads_signif_genes) %in% as.character(c5_msigdb$gene)]
      go_limma_signif_genes <- as.character(limma_signif_genes)[as.character(limma_signif_genes) %in% as.character(c5_msigdb$gene)]
      
      txt <- paste0("... ", hicds, " - ", exprs, " available annot. for universe_genes:\t", length(go_universe_genes) , "/", length(universe_genes), "\n")
      txt <- paste0("... ", hicds, " - ", exprs, " available annot. for tads_signif_genes:\t", length(go_tads_signif_genes) , "/", length(tads_signif_genes), "\n")
      txt <- paste0("... ", hicds, " - ", exprs, " available annot. for limma_signif_genes:\t", length(go_limma_signif_genes) , "/", length(limma_signif_genes), "\n")
      
      stopifnot(go_tads_signif_genes %in% go_universe_genes)
      stopifnot(go_limma_signif_genes %in% go_universe_genes)
      
      #***** 1) TADs signif genes
      cat(paste0(">  start enricher for TADs signif genes \n"))
      
      if(length(go_tads_signif_genes) > 0) {
        
        tad_signif_enrich <- enricher(gene = go_tads_signif_genes, 
                                            TERM2GENE=c5_msigdb,
                                            universe = go_universe_genes,
                                            pvalueCutoff = enricher_pvalueCutoff, 
                                            pAdjustMethod = enricher_pAdjustMethod, 
                                            minGSSize = enricher_minGSSize, 
                                            maxGSSize = enricher_maxGSSize, 
                                            qvalueCutoff =enricher_qvalueCutoff)
        
        tad_signif_enrich_resultDT <- tad_signif_enrich@result
        tad_signif_enrich_resultDT <- tad_signif_enrich_resultDT[order(tad_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        tad_signif_enrich_resultDT$log10_pval <- -log10(tad_signif_enrich_resultDT[,paste0(padjVarGO)])
        tad_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(tad_signif_enrich_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        tad_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(tad_signif_enrich_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        tad_signif_enrich_resultDT$foldEnrichment <- tad_signif_enrich_resultDT$geneRatio/tad_signif_enrich_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(tad_signif_enrich_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds,  " (TADs signif. genes)")
            tad_signif_enrich_resultDT <- tad_signif_enrich_resultDT[order(tad_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_tads_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(tad_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", tad_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        tad_signif_enrich_resultDT <- data.frame(tad_signif_enrich_resultDT)
        
      } else {
        tad_signif_enrich_resultDT <- NULL
      }
      
      
      #***** 2) limma signif. genes
      
      cat(paste0("> start enricher for not_conserved_signif \n"))
      
      
      if(length(go_limma_signif_genes) > 0) {
        
        limma_signif_enrich <- enricher(gene = go_limma_signif_genes, 
                                                TERM2GENE=c5_msigdb,
                                                universe = go_universe_genes,
                                                pvalueCutoff = enricher_pvalueCutoff, 
                                                pAdjustMethod = enricher_pAdjustMethod, 
                                                minGSSize = enricher_minGSSize, 
                                                maxGSSize = enricher_maxGSSize, 
                                                qvalueCutoff =enricher_qvalueCutoff)
        
        limma_signif_enrich_resultDT <- limma_signif_enrich@result
        limma_signif_enrich_resultDT <- limma_signif_enrich_resultDT[order(limma_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        limma_signif_enrich_resultDT$log10_pval <- -log10(limma_signif_enrich_resultDT[,paste0(padjVarGO)])
        limma_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(limma_signif_enrich_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        limma_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(limma_signif_enrich_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        limma_signif_enrich_resultDT$foldEnrichment <- limma_signif_enrich_resultDT$geneRatio/limma_signif_enrich_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(limma_signif_enrich_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds, " (limma signif. genes)")
            limma_signif_enrich_resultDT <- limma_signif_enrich_resultDT[order(limma_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_limma_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(limma_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", limma_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        limma_signif_enrich_resultDT <- data.frame(limma_signif_enrich_resultDT)
        # limma_signif_enrich_resultDT$hicds <- hicds
        # limma_signif_enrich_resultDT$exprds <- exprds
        
        
      } else {
        limma_signif_enrich_resultDT <- NULL
      }
      list(
        tad_signif_enrich_resultDT=tad_signif_enrich_resultDT,
        limma_signif_enrich_resultDT=limma_signif_enrich_resultDT
      )
    } # end-foreach iterating over exprds
    names(exprds_list) <- file.path(hicds, all_exprds[[paste0(hicds)]])
    exprds_list
  } # end-foreach iterating over hicds
  all_go_enrich_list <- unlist(all_go_enrich_list, recursive=FALSE)
  outFile <- file.path(outFolder, paste0("all_go_enrich_list.Rdata"))
  save(all_go_enrich_list, file = outFile, version=2)
  stopifnot(length(all_go_enrich_list) == length(all_datasets))
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0( "all_go_enrich_list.Rdata"))
  cat("... load data\n")
  all_go_enrich_list <- get(load(outFile))
}



######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


  
  