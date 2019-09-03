startTime <- Sys.time()
cat(paste0("> Rscript cmp_DE_TADs_genes_GO.R\n"))

# Rscript cmp_DE_TADs_genes_GO.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS
# Rscript cmp_DE_TADs_genes_GO.R


options(scipen=100)

# suppressPackageStartupMessages(library(reshape, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(ontologyIndex, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
# suppressPackageStartupMessages(library(AnnotationDbi, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
# suppressPackageStartupMessages(library(GO.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # data(go) and minimal_set
# suppressPackageStartupMessages(library(topGO, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) #inducedGraph
# suppressPackageStartupMessages(library(igraph, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(VennDiagram, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(gridExtra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(grid, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

buildData <- FALSE

script0_name <- "0_prepGeneData"

outFolder <- file.path("CMP_DE_TADS_GENES_GO")
dir.create(outFolder, recursive = TRUE)

pipFolder <- file.path(".")

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

hicds="Panc1_rep12_40kb"
exprds="TCGApaad_wt_mutKRAS"

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10
myHeight <- ifelse(plotType == "png", 300, 7)
myWidth <- myHeight

  
TAD_pvalThresh <- 0.01
gene_pvalThresh <- 0.05

padjVarGO <- "p.adjust" # p.adjust or qvalue ???

logFile <- file.path(outFolder, "cmp_DE_TADs_genes_GO_logFile.txt")
if(buildData) file.remove(logFile)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}


final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))


setDir <- "/media/electron"
setDir <- ""

enricher_ontologyType <- "BP"
enricher_pvalueCutoff <- 1
enricher_pAdjustMethod <- "BH"
enricher_minGSSize <- 1
enricher_maxGSSize <- 500
enricher_qvalueCutoff <- 1
enricher_results_sortGOby <- "p.adjust"


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
txt <- paste0("... TAD_pvalThresh:\t", TAD_pvalThresh, "\n")
printAndLog(txt, logFile)
txt <- paste0("... gene_pvalThresh:\t", gene_pvalThresh, "\n")
printAndLog(txt, logFile)

txt <- paste0("... plotMaxBars:\t", plotMaxBars, "\n")
printAndLog(txt, logFile)



if(buildData) {

  
  all_enricher_result <- foreach(hicds = all_hicds) %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    
    exprds_data <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
  
      
      cat("... start for :", hicds, " - ", exprds, "\n")
      
      ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
      
      stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
      geneListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneListFile))
      pipeline_geneList <- eval(parse(text = load(geneListFile))) 
      
      
      
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
      g2t_DT <- g2t_DT[ g2t_DT$entrezID %in% pipeline_geneList,]
  
          
      
      # signif. TADs pval comb 0.01 & FDR=0.2
      signifTADs <- final_dt$region[final_dt$hicds == hicds & final_dt$exprds == exprds & final_dt$signifFDR_0.2 & final_dt$adjPvalComb <= TAD_pvalThresh]
      
      genes_signifTADs <- g2t_DT$entrezID[g2t_DT$region %in% signifTADs]
      
      
      
      # retrieve signif. genes
      script1_name <- "1_runGeneDE"
      topTable_DT_file <- file.path(pipOutFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")
      topTable_DT <- get(load(topTable_DT_file))
      topTable_DT$genes <- as.character(topTable_DT$genes)
      stopifnot(names(pipeline_geneList) %in% topTable_DT$genes)
      topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipeline_geneList),]
      genes_signif <- pipeline_geneList[topTable_DT$genes[topTable_DT$adj.P.Val <= gene_pvalThresh]]
      stopifnot(!is.na(genes_signif))
      
      
      
      go_genes_signifTADs <- as.character(genes_signifTADs)[as.character(genes_signifTADs) %in% as.character(c5_msigdb$gene)]
      
  
      
      go_genes_signif <- as.character(genes_signif)[as.character(genes_signif) %in% as.character(c5_msigdb$gene)]
      
      
      universe_genes <- as.character(pipeline_geneList)[as.character(pipeline_geneList) %in% as.character(pipeline_geneList)]
      
      cat("... available annot. for genes_signifTADs:\t", length(go_genes_signifTADs) , "/", length(genes_signifTADs), "\n")
      cat("... available annot. for genes_signif:\t", length(go_genes_signif) , "/", length(genes_signif), "\n")
      cat("... available annot. for universe_genes:\t", length(universe_genes) , "/", length(pipeline_geneList), "\n")
      
      
      stopifnot(go_genes_signif %in% universe_genes)
      stopifnot(go_genes_signifTADs %in% universe_genes)
      stopifnot( length(go_genes_signif) <= length(universe_genes) )
      stopifnot( length(go_genes_signifTADs) <= length(universe_genes) )
      
      
      #***** 0b) genes_onlySignif
      
      go_genes_onlySignif <- go_genes_signif[! go_genes_signif %in% go_genes_signifTADs]
      
      cat(paste0("... ", hicds, " - ", exprds, " - start enricher for genes_signifTADs \n"))
      
      if(length(go_genes_onlySignif) > 0) {
        
        genes_onlySignif_enrich <- enricher(gene = go_genes_onlySignif, 
                                            TERM2GENE=c5_msigdb,
                                            universe = universe_genes,
                                            pvalueCutoff = enricher_pvalueCutoff, 
                                            pAdjustMethod = enricher_pAdjustMethod, 
                                            minGSSize = enricher_minGSSize, 
                                            maxGSSize = enricher_maxGSSize, 
                                            qvalueCutoff =enricher_qvalueCutoff)
        
        
        
        
        genes_onlySignif_enrich_resultDT <- genes_onlySignif_enrich@result
        genes_onlySignif_enrich_resultDT <- genes_onlySignif_enrich_resultDT[order(genes_onlySignif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        genes_onlySignif_enrich_resultDT$log10_pval <- -log10(genes_onlySignif_enrich_resultDT[,paste0(padjVarGO)])
        genes_onlySignif_enrich_resultDT$geneRatio <- as.numeric(sapply(genes_onlySignif_enrich_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        genes_onlySignif_enrich_resultDT$bgRatio <- as.numeric(sapply(genes_onlySignif_enrich_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        genes_onlySignif_enrich_resultDT$foldEnrichment <- genes_onlySignif_enrich_resultDT$geneRatio/genes_onlySignif_enrich_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(genes_onlySignif_enrich_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds, " (genes_onlySignif)")
            genes_onlySignif_enrich_resultDT <- genes_onlySignif_enrich_resultDT[order(genes_onlySignif_enrich_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_genes_onlySignif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(genes_onlySignif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", genes_onlySignif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        genes_onlySignif_enrich_resultDT <- data.frame(genes_onlySignif_enrich_resultDT)
        genes_onlySignif_enrich_resultDT$hicds <- hicds
        genes_onlySignif_enrich_resultDT$exprds <- exprds
        
        
      } else {
        genes_onlySignif_enrich_resultDT <- NULL
      }
      
      
      
      #***** 0a) genes_onlySignifTADs
      
      go_genes_onlySignifTADs <- go_genes_signifTADs[! go_genes_signifTADs %in% go_genes_signif]
      
      cat(paste0("... ", hicds, " - ", exprds, " - start enricher for genes_signifTADs \n"))
      
      if(length(go_genes_onlySignifTADs) > 0) {
        
      
        
        genes_onlySignifTADs_enrich <- enricher(gene = go_genes_onlySignifTADs, 
                                                TERM2GENE=c5_msigdb,
                                                universe = universe_genes,
                                                pvalueCutoff = enricher_pvalueCutoff, 
                                                pAdjustMethod = enricher_pAdjustMethod, 
                                                minGSSize = enricher_minGSSize, 
                                                maxGSSize = enricher_maxGSSize, 
                                                qvalueCutoff =enricher_qvalueCutoff)
        
        
        
        
        genes_onlySignifTADs_enrich_resultDT <- genes_onlySignifTADs_enrich@result
        genes_onlySignifTADs_enrich_resultDT <- genes_onlySignifTADs_enrich_resultDT[order(genes_onlySignifTADs_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        genes_onlySignifTADs_enrich_resultDT$log10_pval <- -log10(genes_onlySignifTADs_enrich_resultDT[,paste0(padjVarGO)])
        genes_onlySignifTADs_enrich_resultDT$geneRatio <- as.numeric(sapply(genes_onlySignifTADs_enrich_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        genes_onlySignifTADs_enrich_resultDT$bgRatio <- as.numeric(sapply(genes_onlySignifTADs_enrich_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        genes_onlySignifTADs_enrich_resultDT$foldEnrichment <- genes_onlySignifTADs_enrich_resultDT$geneRatio/genes_onlySignifTADs_enrich_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(genes_onlySignifTADs_enrich_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds, " (genes_onlySignifTADs)")
            genes_onlySignifTADs_enrich_resultDT <- genes_onlySignifTADs_enrich_resultDT[order(genes_onlySignifTADs_enrich_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_genes_onlySignifTADs", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(genes_onlySignifTADs_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", genes_onlySignifTADs_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        genes_onlySignifTADs_enrich_resultDT <- data.frame(genes_onlySignifTADs_enrich_resultDT)
        genes_onlySignifTADs_enrich_resultDT$hicds <- hicds
        genes_onlySignifTADs_enrich_resultDT$exprds <- exprds
        
        
      } else {
        genes_onlySignifTADs_enrich_resultDT <- NULL
      }
      
      
      #***** 1) genes_signifTADs
  
      cat(paste0("... ", hicds, " - ", exprds, " - start enricher for genes_signifTADs \n"))
      
      if(length(go_genes_signifTADs) > 0) {
        
        
        
        genes_signifTADs_enrich <- enricher(gene = go_genes_signifTADs, 
                                         TERM2GENE=c5_msigdb,
                                         universe = universe_genes,
                                         pvalueCutoff = enricher_pvalueCutoff, 
                                         pAdjustMethod = enricher_pAdjustMethod, 
                                         minGSSize = enricher_minGSSize, 
                                         maxGSSize = enricher_maxGSSize, 
                                         qvalueCutoff =enricher_qvalueCutoff)
        
        
        
        
        genes_signifTADs_enrich_resultDT <- genes_signifTADs_enrich@result
        genes_signifTADs_enrich_resultDT <- genes_signifTADs_enrich_resultDT[order(genes_signifTADs_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        genes_signifTADs_enrich_resultDT$log10_pval <- -log10(genes_signifTADs_enrich_resultDT[,paste0(padjVarGO)])
        genes_signifTADs_enrich_resultDT$geneRatio <- as.numeric(sapply(genes_signifTADs_enrich_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        genes_signifTADs_enrich_resultDT$bgRatio <- as.numeric(sapply(genes_signifTADs_enrich_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        genes_signifTADs_enrich_resultDT$foldEnrichment <- genes_signifTADs_enrich_resultDT$geneRatio/genes_signifTADs_enrich_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(genes_signifTADs_enrich_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds, " (genes_signifTADs)")
            genes_signifTADs_enrich_resultDT <- genes_signifTADs_enrich_resultDT[order(genes_signifTADs_enrich_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_genes_signifTADs", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(genes_signifTADs_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", genes_signifTADs_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        genes_signifTADs_enrich_resultDT <- data.frame(genes_signifTADs_enrich_resultDT)
        genes_signifTADs_enrich_resultDT$hicds <- hicds
        genes_signifTADs_enrich_resultDT$exprds <- exprds
        
      } else {
        genes_signifTADs_enrich_resultDT <- NULL
      }
        
        
      #***** 2) genes_signif
      
      cat(paste0("... ", hicds, " - ", exprds, " - start enricher for genes_signif \n"))
      
      
      if(length(go_genes_signif) > 0){
        genes_signif_enrich <- enricher(gene = go_genes_signif, 
                                            TERM2GENE=c5_msigdb,
                                            universe = universe_genes,
                                            pvalueCutoff = enricher_pvalueCutoff, 
                                            pAdjustMethod = enricher_pAdjustMethod, 
                                            minGSSize = enricher_minGSSize, 
                                            maxGSSize = enricher_maxGSSize, 
                                            qvalueCutoff =enricher_qvalueCutoff)
        
        genes_signif_enrich_resultDT <- genes_signif_enrich@result
        genes_signif_enrich_resultDT <- genes_signif_enrich_resultDT[order(genes_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        genes_signif_enrich_resultDT$log10_pval <- -log10(genes_signif_enrich_resultDT[,paste0(padjVarGO)])
        genes_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(genes_signif_enrich_resultDT$GeneRatio, function(x) {
                                                                          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        genes_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(genes_signif_enrich_resultDT$BgRatio, function(x) {
                                                                          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        genes_signif_enrich_resultDT$foldEnrichment <- genes_signif_enrich_resultDT$geneRatio/genes_signif_enrich_resultDT$bgRatio
            
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(genes_signif_enrich_resultDT)))
        if(genes_signif_plotMax > 0) {
            for(var_plot in barplot_vars) {
              myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds, " (genes_signif)")
              genes_signif_enrich_resultDT <- genes_signif_enrich_resultDT[order(genes_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
              outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_genes_signif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
              do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
              par(oma=c(10,1,1,1))
              barplot(genes_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                      main = myTit, 
                      # ylab = paste0("-log10 (", padjVarGO, ")"),
                      ylab = paste0(barplot_vars_tit[var_plot]),
                      names.arg = gsub("GO_", "", genes_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
                      cex.names = 0.6
              )
              foo <- dev.off()
              cat(paste0("... written: ", outFile, "\n"))
            }
        }
        genes_signif_enrich_resultDT <- data.frame(genes_signif_enrich_resultDT)
        genes_signif_enrich_resultDT$hicds <- hicds
        genes_signif_enrich_resultDT$exprds <- exprds
        save(genes_signif_enrich_resultDT, file="genes_signif_enrich_resultDT.Rdata", version=2)
        
      } else {
        genes_signif_enrich_resultDT <- NULL 
      }
        
      
      
      list(
        genes_signifTADs_enrich_resultDT=genes_signifTADs_enrich_resultDT,
        genes_onlySignifTADs_enrich_resultDT=genes_onlySignifTADs_enrich_resultDT,
        genes_onlySignif_enrich_resultDT= genes_onlySignif_enrich_resultDT,
        genes_signif_enrich_resultDT= genes_signif_enrich_resultDT
      )
      
      
      
      
    } # end-for iterating over hicds
    names(exprds_data) <- all_exprds[[paste0(hicds)]]
    exprds_data
    
  } # end-foreach iterating over exprcds
  
  
  names(all_enricher_result) <- all_hicds
  
  outFile <- file.path(outFolder, "all_enricher_result.Rdata")
  save(all_enricher_result, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_enricher_result.Rdata")  
  all_enricher_result <- get(load(outFile))
}


signif_go_pvalThresh <- -log10(0.05)

all_dt=c(
"genes_signifTADs_enrich_resultDT",
"genes_onlySignifTADs_enrich_resultDT",
"genes_onlySignif_enrich_resultDT",
"genes_signif_enrich_resultDT")

hicds=all_hicds[1]
dt = all_dt[1]

topCommonBars <- 10

for(dt in all_dt) {
  all_go_categories <- foreach(hicds = all_hicds) %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_go_categories <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      curr_dt <- all_enricher_result[[paste0(hicds)]][[paste0(exprds)]][[paste0(dt)]]
      rownames(curr_dt)[curr_dt[,paste0(padjVarGO)] <= padjVarGO_plotThresh]
    } 
    names(exprds_go_categories) <- all_exprds[[paste0(hicds)]]
    exprds_go_categories
  }
  go_categories <- unlist(all_go_categories)
  go_categories_count <- setNames(as.numeric(table(go_categories)), names(table(go_categories)))
  go_categories_count <- sort(go_categories_count, decreasing = TRUE)
  outFile <- file.path(outFolder,paste0("all_ds_", dt, "_intersect_",padjVarGO, "_barplot", ".", plotType))
  do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
  par(oma=c(10,1,1,1))
  barplot(go_categories_count[1:topCommonBars], las=2, 
          ylab="# of datasets",
          names.arg=gsub("GO_", "", names(go_categories_count[1:topCommonBars])),
          main=paste0(gsub("_resultDT", "", dt), " - intersect across DS"),
          cex.names=0.6
          )
  
  mtext(side=3, text=paste0(padjVarGO, "<=",padjVarGO_plotThresh, " (top ", topCommonBars, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

all_cmp_names <- unique(all_cmps)
cmp=all_cmp_names[1]
for(cmp in all_cmp_names){
  
  for(dt in all_dt) {
    all_go_categories <- foreach(hicds = all_hicds) %dopar% {
      exprds = all_exprds[[paste0(hicds)]][1]
      exprds_go_categories <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
        curr_dt <- all_enricher_result[[paste0(hicds)]][[paste0(exprds)]][[paste0(dt)]]
        cmpType <- all_cmps[paste0(exprds)]
        stopifnot(!is.na(cmpType))
        if(cmpType != cmp) return(NULL)
        
        rownames(curr_dt)[curr_dt[,paste0(padjVarGO)] <= padjVarGO_plotThresh]
      } 
      names(exprds_go_categories) <- all_exprds[[paste0(hicds)]]
      exprds_go_categories
    }
    go_categories <- unlist(all_go_categories)
    go_categories_count <- setNames(as.numeric(table(go_categories)), names(table(go_categories)))
    go_categories_count <- sort(go_categories_count, decreasing = TRUE)
    outFile <- file.path(outFolder,paste0("all_ds_",cmp, "_", dt, "_intersect_",padjVarGO, "_barplot", ".", plotType))
    do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
    par(oma=c(10,1,1,1))
    barplot(go_categories_count[1:topCommonBars], las=2, 
            names.arg=gsub("GO_", "", names(go_categories_count[1:topCommonBars])),
            ylab="# of datasets",
            main=paste0(gsub("_resultDT", "", dt), " - ", cmp),
            cex.names=0.6
    )
    mtext(side=3, text=paste0(padjVarGO, "<=",padjVarGO_plotThresh, " (top ", topCommonBars, ")"))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

