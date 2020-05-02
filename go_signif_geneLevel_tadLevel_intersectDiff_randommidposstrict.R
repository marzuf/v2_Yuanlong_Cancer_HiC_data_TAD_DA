options(scipen=100)

SSHFS=F

# Rscript go_signif_geneLevel_tadLevel_intersectDiff_randommidposstrict.R <tadpval> <genepval>
# Rscript go_signif_geneLevel_tadLevel_intersectDiff_randommidposstrict.R 0.01 0.05 
# Rscript go_signif_geneLevel_tadLevel_intersectDiff_randommidposstrict.R 0.05 0.05 
# Rscript go_signif_geneLevel_tadLevel_intersectDiff_randommidposstrict.R 0.01 0.01 

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "go_signif_geneLevel_tadLevel_intersectDiff_randommidposstrict.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(clusterProfiler)
require(foreach)
require(doMC)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)
require(ggpubr)
require(ggplot2)
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

minOcc <- 5


col1 <- get_palette("Dark2", 5)[1]
col2 <- get_palette("Dark2", 5)[2]
col3 <- get_palette("Dark2", 5)[3]
col4 <- get_palette("Dark2", 5)[4]
col5 <- get_palette("Dark2", 5)[5]



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

all_hicds <- all_hicds[  grepl("RANDOMMIDPOSSTRICT", all_hicds)]

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds


all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

tads_signifThresh <- 0.05
genes_signifThresh <- 0.05
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
tads_signifThresh <- as.numeric(args[1])
genes_signifThresh <- as.numeric(args[2])
stopifnot(!is.na(genes_signifThresh))
stopifnot(!is.na(tads_signifThresh))

final_dt_file <- file.path("CREATE_FINAL_TABLE_RANDOM", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))
signif_column <- "adjPvalComb"

signifcol <- paste0(signif_column, "_", tads_signifThresh)
final_dt[, paste0(signifcol)] <- final_dt[, paste0(signif_column)] <= tads_signifThresh
cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> tads_signifThresh\t=\t", tads_signifThresh, "\n"))

padjVarGO <- "p.adjust" # p.adjust or qvalue ???

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

cat(paste0("tads_signifThresh = ", tads_signifThresh, "\n"))
cat(paste0("genes_signifThresh = ", genes_signifThresh, "\n"))
#stop("-ok\n")

outFolder <- file.path("GO_SIGNIF_GENELEVEL_TADLEVEL_INTERSECTDIFF_RANDOMMIDPOSSTRICT", paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifThresh))
dir.create(outFolder, recursive = TRUE)
logFile <- file.path(outFolder, "go_signif_geneLevel_tadLevel_intersectDiff_logFile.txt")
if(buildTable) file.remove(logFile)

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
      
      limma_signif_genes <- topTable_DT$entrezID[topTable_DT$adj.P.Val <= genes_signifThresh]
      stopifnot(limma_signif_genes %in% pipeline_geneList)
      
      
      universe_genes <- pipeline_geneList
      
      go_universe_genes <- as.character(universe_genes)[as.character(universe_genes) %in% as.character(c5_msigdb$gene)]
      go_tads_signif_genes <- as.character(tads_signif_genes)[as.character(tads_signif_genes) %in% as.character(c5_msigdb$gene)]
      go_limma_signif_genes <- as.character(limma_signif_genes)[as.character(limma_signif_genes) %in% as.character(c5_msigdb$gene)]
      
      intersect_signif_genes <- intersect(tads_signif_genes, limma_signif_genes)
      limmaOnly_signif_genes <- setdiff(limma_signif_genes, tads_signif_genes); stopifnot(limmaOnly_signif_genes %in% limma_signif_genes)
      tadsOnly_signif_genes <- setdiff(tads_signif_genes, limma_signif_genes); stopifnot(tadsOnly_signif_genes %in% tads_signif_genes)
      
      go_intersect_signif_genes <- intersect(go_tads_signif_genes, go_limma_signif_genes)
      go_limmaOnly_signif_genes <- setdiff(go_limma_signif_genes, go_tads_signif_genes); stopifnot(go_limmaOnly_signif_genes %in% go_limma_signif_genes)
      go_tadsOnly_signif_genes <- setdiff(go_tads_signif_genes, go_limma_signif_genes); stopifnot(go_tadsOnly_signif_genes %in% go_tads_signif_genes)
      stopifnot(setequal(go_limmaOnly_signif_genes, limmaOnly_signif_genes[limmaOnly_signif_genes %in% as.character(c5_msigdb$gene)]))
      stopifnot(setequal(go_tadsOnly_signif_genes, tadsOnly_signif_genes[tadsOnly_signif_genes %in% as.character(c5_msigdb$gene)]))
      
      txt <- paste0("... ", hicds, " - ", exprds, " available annot. for universe_genes:\t", length(go_universe_genes) , "/", length(universe_genes), "\n")
      txt <- paste0("... ", hicds, " - ", exprds, " available annot. for tads_signif_genes:\t", length(go_tads_signif_genes) , "/", length(tads_signif_genes), "\n")
      txt <- paste0("... ", hicds, " - ", exprds, " available annot. for limma_signif_genes:\t", length(go_limma_signif_genes) , "/", length(limma_signif_genes), "\n")
      txt <- paste0("... ", hicds, " - ", exprds, " available annot. for intersect_signif_genes:\t", length(intersect_signif_genes) , "/", length(go_intersect_signif_genes), "\n")
      txt <- paste0("... ", hicds, " - ", exprds, " available annot. for limmaOnly_signif_genes:\t", length(go_limmaOnly_signif_genes) , "/", length(limmaOnly_signif_genes), "\n")
      txt <- paste0("... ", hicds, " - ", exprds, " available annot. for tadsOnly_signif_genes:\t", length(go_tadsOnly_signif_genes) , "/", length(tadsOnly_signif_genes), "\n")
      
      stopifnot(go_tads_signif_genes %in% go_universe_genes)
      stopifnot(go_limma_signif_genes %in% go_universe_genes)
      stopifnot(go_intersect_signif_genes %in% go_universe_genes)
      stopifnot(go_limmaOnly_signif_genes %in% go_universe_genes)
      stopifnot(go_tadsOnly_signif_genes %in% go_universe_genes)
      
      
      #***** A) limma only signif genes
      cat(paste0(">  start enricher for limma only signif genes \n"))
      
      if(length(go_limmaOnly_signif_genes) > 0) {
        
        limmaOnly_signif_enrich <- enricher(gene = go_limmaOnly_signif_genes, 
                                            TERM2GENE=c5_msigdb,
                                            universe = go_universe_genes,
                                            pvalueCutoff = enricher_pvalueCutoff, 
                                            pAdjustMethod = enricher_pAdjustMethod, 
                                            minGSSize = enricher_minGSSize, 
                                            maxGSSize = enricher_maxGSSize, 
                                            qvalueCutoff =enricher_qvalueCutoff)
        
        limmaOnly_signif_enrich_resultDT <- limmaOnly_signif_enrich@result
        limmaOnly_signif_enrich_resultDT <- limmaOnly_signif_enrich_resultDT[order(limmaOnly_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        limmaOnly_signif_enrich_resultDT$log10_pval <- -log10(limmaOnly_signif_enrich_resultDT[,paste0(padjVarGO)])
        limmaOnly_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(limmaOnly_signif_enrich_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        limmaOnly_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(limmaOnly_signif_enrich_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        limmaOnly_signif_enrich_resultDT$foldEnrichment <- limmaOnly_signif_enrich_resultDT$geneRatio/limmaOnly_signif_enrich_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(limmaOnly_signif_enrich_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds,  " (limma only signif. genes)")
            limmaOnly_signif_enrich_resultDT <- limmaOnly_signif_enrich_resultDT[order(limmaOnly_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_limmaOnly_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(limmaOnly_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", limmaOnly_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        limmaOnly_signif_enrich_resultDT <- data.frame(limmaOnly_signif_enrich_resultDT)
        
      } else {
        limmaOnly_signif_enrich_resultDT <- NULL
      }
      
      
      #***** B) TADs only signif genes
      cat(paste0(">  start enricher for TADs only signif genes \n"))
      
      if(length(go_tadsOnly_signif_genes) > 0) {
        
        tadsOnly_signif_enrich <- enricher(gene = go_tadsOnly_signif_genes, 
                                           TERM2GENE=c5_msigdb,
                                           universe = go_universe_genes,
                                           pvalueCutoff = enricher_pvalueCutoff, 
                                           pAdjustMethod = enricher_pAdjustMethod, 
                                           minGSSize = enricher_minGSSize, 
                                           maxGSSize = enricher_maxGSSize, 
                                           qvalueCutoff =enricher_qvalueCutoff)
        
        tadsOnly_signif_enrich_resultDT <- tadsOnly_signif_enrich@result
        tadsOnly_signif_enrich_resultDT <- tadsOnly_signif_enrich_resultDT[order(tadsOnly_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        tadsOnly_signif_enrich_resultDT$log10_pval <- -log10(tadsOnly_signif_enrich_resultDT[,paste0(padjVarGO)])
        tadsOnly_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(tadsOnly_signif_enrich_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        tadsOnly_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(tadsOnly_signif_enrich_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        tadsOnly_signif_enrich_resultDT$foldEnrichment <- tadsOnly_signif_enrich_resultDT$geneRatio/tadsOnly_signif_enrich_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(tadsOnly_signif_enrich_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds,  " (TADs only signif. genes)")
            tadsOnly_signif_enrich_resultDT <- tadsOnly_signif_enrich_resultDT[order(tadsOnly_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_tadsOnly_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(tadsOnly_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", tadsOnly_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        tadsOnly_signif_enrich_resultDT <- data.frame(tadsOnly_signif_enrich_resultDT)
        
      } else {
        tadsOnly_signif_enrich_resultDT <- NULL
      }
      
      
      #***** C) intersect signif genes
      cat(paste0(">  start enricher for intersect signif genes \n"))
      
      if(length(go_intersect_signif_genes) > 0) {
        
        intersect_signif_enrich <- enricher(gene = go_intersect_signif_genes, 
                                            TERM2GENE=c5_msigdb,
                                            universe = go_universe_genes,
                                            pvalueCutoff = enricher_pvalueCutoff, 
                                            pAdjustMethod = enricher_pAdjustMethod, 
                                            minGSSize = enricher_minGSSize, 
                                            maxGSSize = enricher_maxGSSize, 
                                            qvalueCutoff =enricher_qvalueCutoff)
        
        intersect_signif_enrich_resultDT <- intersect_signif_enrich@result
        intersect_signif_enrich_resultDT <- intersect_signif_enrich_resultDT[order(intersect_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        intersect_signif_enrich_resultDT$log10_pval <- -log10(intersect_signif_enrich_resultDT[,paste0(padjVarGO)])
        intersect_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(intersect_signif_enrich_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        intersect_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(intersect_signif_enrich_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        intersect_signif_enrich_resultDT$foldEnrichment <- intersect_signif_enrich_resultDT$geneRatio/intersect_signif_enrich_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(intersect_signif_enrich_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds,  " (intersect signif. genes)")
            intersect_signif_enrich_resultDT <- intersect_signif_enrich_resultDT[order(intersect_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_intersect_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(intersect_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", intersect_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        intersect_signif_enrich_resultDT <- data.frame(intersect_signif_enrich_resultDT)
        
      } else {
        intersect_signif_enrich_resultDT <- NULL
      }
      
      
      
      
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
        tads_signif_genes=tads_signif_genes,
        limma_signif_genes=limma_signif_genes,
        tadsOnly_signif_genes=tadsOnly_signif_genes,
        limmaOnly_signif_genes=limmaOnly_signif_genes,
        universe_genes = universe_genes,
        tad_signif_enrich_resultDT=tad_signif_enrich_resultDT,
        tadsOnly_signif_enrich_resultDT=tadsOnly_signif_enrich_resultDT,
        limma_signif_enrich_resultDT=limma_signif_enrich_resultDT,
        limmaOnly_signif_enrich_resultDT=limmaOnly_signif_enrich_resultDT,
        intersect_signif_enrich_resultDT=intersect_signif_enrich_resultDT
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
# outFile="GO_SIGNIF_GENELEVEL_TADLEVEL_INTERSECTDIFF/tadPvalThresh0.05_genePvalThresh0.05/all_go_enrich_list.Rdata"
# outFile="GO_SIGNIF_GENELEVEL_TADLEVEL_INTERSECTDIFF/tadPvalThresh0.01_genePvalThresh0.05/all_go_enrich_list.Rdata"
# all_go_enrich_list <- get(load(outFile))



### GO ENRICHMENT ALL DS LEVEL
# all_tads_signif_genes <- lapply(all_go_enrich_list, function(x) x[["tads_signif_genes"]])
# all_limma_signif_genes <- lapply(all_go_enrich_list, function(x) x[["limma_signif_genes"]])
# all_tadsOnly_signif_genes <- lapply(all_go_enrich_list, function(x) x[["tadsOnly_signif_genes"]])
# all_limmaOnly_signif_genes <- lapply(all_go_enrich_list, function(x) x[["limmaOnly_signif_genes"]])

all_universe_genes <- lapply(all_go_enrich_list, function(x) x[["universe_genes"]])
all_universe_genes <- unique(unlist(all_universe_genes))

all_tads_signif_genes_dt <- do.call(rbind, lapply(1:length(all_go_enrich_list), function(x)  {
  matching_genes <- as.character(all_go_enrich_list[[x]][["tads_signif_genes"]])
  if(length(matching_genes) == 0) {
    dataset = character(0)
    hicds = character(0)
    exprds = character(0)
  } else {
    dataset = names(all_go_enrich_list)[x]
    hicds = dirname(names(all_go_enrich_list)[x])
    exprds = basename(names(all_go_enrich_list)[x])
  }
  data.frame(
  dataset = dataset,
  hicds = hicds,
  exprds = exprds,
  entrezID = matching_genes,
  stringsAsFactors = FALSE
  )}))

all_limma_signif_genes_dt <- do.call(rbind, lapply(1:length(all_go_enrich_list), function(x)  {
  matching_genes <- as.character(all_go_enrich_list[[x]][["limma_signif_genes"]])
  if(length(matching_genes) == 0) {
    dataset = character(0)
    hicds = character(0)
    exprds = character(0)
  } else {
    dataset = names(all_go_enrich_list)[x]
    hicds = dirname(names(all_go_enrich_list)[x])
    exprds = basename(names(all_go_enrich_list)[x])
  }
  data.frame(
    dataset = dataset,
    hicds = hicds,
    exprds = exprds,
    entrezID = matching_genes,
    stringsAsFactors = FALSE
  )}))


all_tadsOnly_signif_genes_dt <- do.call(rbind, lapply(1:length(all_go_enrich_list), function(x)  {
  matching_genes <- as.character(all_go_enrich_list[[x]][["tadsOnly_signif_genes"]])
  if(length(matching_genes) == 0) {
    dataset = character(0)
    hicds = character(0)
    exprds = character(0)
  } else {
    dataset = names(all_go_enrich_list)[x]
    hicds = dirname(names(all_go_enrich_list)[x])
    exprds = basename(names(all_go_enrich_list)[x])
  }
  data.frame(
    dataset = dataset,
    hicds = hicds,
    exprds = exprds,
    entrezID = matching_genes,
    stringsAsFactors = FALSE
  )}))


all_limmaOnly_signif_genes_dt <- do.call(rbind, lapply(1:length(all_go_enrich_list), function(x)  {
  matching_genes <- as.character(all_go_enrich_list[[x]][["limmaOnly_signif_genes"]])
  if(length(matching_genes) == 0) {
    dataset = character(0)
    hicds = character(0)
    exprds = character(0)
  } else {
    dataset = names(all_go_enrich_list)[x]
    hicds = dirname(names(all_go_enrich_list)[x])
    exprds = basename(names(all_go_enrich_list)[x])
  }
  data.frame(
    dataset = dataset,
    hicds = hicds,
    exprds = exprds,
    entrezID = matching_genes,
    stringsAsFactors = FALSE
  )}))



all_intersect_signif_genes_dt <- do.call(rbind, lapply(1:length(all_go_enrich_list), function(x)  {
  limma_genes <- as.character(all_go_enrich_list[[x]][["limma_signif_genes"]])
  tads_genes <- as.character(all_go_enrich_list[[x]][["tads_signif_genes"]])
  matching_genes <- intersect(tads_genes, limma_genes)
  if(length(matching_genes) == 0) {
    dataset = character(0)
    hicds = character(0)
    exprds = character(0)
  } else {
    dataset = names(all_go_enrich_list)[x]
    hicds = dirname(names(all_go_enrich_list)[x])
    exprds = basename(names(all_go_enrich_list)[x])
  }
  data.frame(
    dataset = dataset,
    hicds = hicds,
    exprds = exprds,
    entrezID = matching_genes,
    stringsAsFactors = FALSE
  )}))

all_tads_signif_genes_dt_count <- aggregate(exprds~entrezID, FUN=function(x)length(unique(x)), data = all_tads_signif_genes_dt)
stopifnot(setequal(all_tads_signif_genes_dt$entrezID, all_tads_signif_genes_dt_count$entrezID))
range(all_tads_signif_genes_dt_count$exprds)
all_tads_signif_genes_dt_count <- all_tads_signif_genes_dt_count[all_tads_signif_genes_dt_count$exprds >= minOcc,]

all_tadsOnly_signif_genes_dt_count <- aggregate(exprds~entrezID, FUN=function(x)length(unique(x)), data = all_tadsOnly_signif_genes_dt)
stopifnot(setequal(all_tadsOnly_signif_genes_dt$entrezID, all_tadsOnly_signif_genes_dt_count$entrezID))
range(all_tadsOnly_signif_genes_dt_count$exprds)
all_tadsOnly_signif_genes_dt_count <- all_tadsOnly_signif_genes_dt_count[all_tadsOnly_signif_genes_dt_count$exprds >= minOcc,]

all_limma_signif_genes_dt_count <- aggregate(exprds~entrezID, FUN=function(x)length(unique(x)), data = all_limma_signif_genes_dt)
stopifnot(setequal(all_limma_signif_genes_dt$entrezID, all_limma_signif_genes_dt_count$entrezID))
range(all_limma_signif_genes_dt_count$exprds)
all_limma_signif_genes_dt_count <- all_limma_signif_genes_dt_count[all_limma_signif_genes_dt_count$exprds >= minOcc,]

all_limmaOnly_signif_genes_dt_count <- aggregate(exprds~entrezID, FUN=function(x)length(unique(x)), data = all_limmaOnly_signif_genes_dt)
stopifnot(setequal(all_limmaOnly_signif_genes_dt$entrezID, all_limmaOnly_signif_genes_dt_count$entrezID))
range(all_limmaOnly_signif_genes_dt_count$exprds)
all_limmaOnly_signif_genes_dt_count <- all_limmaOnly_signif_genes_dt_count[all_limmaOnly_signif_genes_dt_count$exprds >= minOcc,]

all_intersect_signif_genes_dt_count <- aggregate(exprds~entrezID, FUN=function(x)length(unique(x)), data = all_intersect_signif_genes_dt)
stopifnot(setequal(all_intersect_signif_genes_dt$entrezID, all_intersect_signif_genes_dt_count$entrezID))
range(all_intersect_signif_genes_dt_count$exprds)
all_intersect_signif_genes_dt_count <- all_intersect_signif_genes_dt_count[all_intersect_signif_genes_dt_count$exprds >= minOcc,]


all_tads_signif_genes <- as.character(all_tads_signif_genes_dt_count$entrezID)
go_all_tads_signif_genes <- all_tads_signif_genes[all_tads_signif_genes %in% as.character(c5_msigdb$gene)]

all_tadsOnly_signif_genes <- as.character(all_tadsOnly_signif_genes_dt_count$entrezID)
go_all_tadsOnly_signif_genes <- all_tadsOnly_signif_genes[all_tadsOnly_signif_genes %in% as.character(c5_msigdb$gene)]

all_limma_signif_genes <- as.character(all_limma_signif_genes_dt_count$entrezID)
go_all_limma_signif_genes <- all_limma_signif_genes[all_limma_signif_genes %in% as.character(c5_msigdb$gene)]

all_limmaOnly_signif_genes <- as.character(all_limmaOnly_signif_genes_dt_count$entrezID)
go_all_limmaOnly_signif_genes <- all_limmaOnly_signif_genes[all_limmaOnly_signif_genes%in% as.character(c5_msigdb$gene)]

all_intersect_signif_genes <- as.character(all_intersect_signif_genes_dt_count$entrezID)
go_all_intersect_signif_genes <- all_intersect_signif_genes[all_intersect_signif_genes %in% as.character(c5_msigdb$gene)]

go_all_universe_genes <- as.character(all_universe_genes)[as.character(all_universe_genes) %in% as.character(c5_msigdb$gene)]

txt <- paste0("... all DS",  " available annot. for universe_genes:\t", length(go_all_universe_genes) , "/", length(all_universe_genes), "\n")
txt <- paste0("... all DS",  " available annot. for all_tads_signif_genes:\t", length(go_all_tads_signif_genes) , "/", length(all_tads_signif_genes), "\n")
txt <- paste0("... all DS", " available annot. for tadsOnly_signif_genes:\t", length(go_all_tadsOnly_signif_genes) , "/", length(all_tadsOnly_signif_genes), "\n")
txt <- paste0("... all DS",  " available annot. for limma_signif_genes:\t", length(go_all_limma_signif_genes) , "/", length(all_limma_signif_genes), "\n")
txt <- paste0("... all DS",  " available annot. for limmaOnly_signif_genes:\t", length(go_all_limmaOnly_signif_genes) , "/", length(all_limmaOnly_signif_genes), "\n")
txt <- paste0("... all DS",  " available annot. for intersect_signif_genes:\t", length(go_all_intersect_signif_genes) , "/", length(all_intersect_signif_genes), "\n")

stopifnot(go_all_tads_signif_genes %in% go_all_universe_genes)
stopifnot(go_all_limma_signif_genes %in% go_all_universe_genes)
stopifnot(go_all_intersect_signif_genes %in% go_all_universe_genes)
stopifnot(go_all_limmaOnly_signif_genes %in% go_all_universe_genes)
stopifnot(go_all_tadsOnly_signif_genes %in% go_all_universe_genes)



#***** A) limma only signif genes
cat(paste0(">  start enricher for limma only signif genes \n"))

if(length(go_all_limmaOnly_signif_genes) > 0) {
  
  all_limmaOnly_signif_enrich <- enricher(gene = go_all_limmaOnly_signif_genes, 
                                          TERM2GENE=c5_msigdb,
                                          universe = go_all_universe_genes,
                                          pvalueCutoff = enricher_pvalueCutoff, 
                                          pAdjustMethod = enricher_pAdjustMethod, 
                                          minGSSize = enricher_minGSSize, 
                                          maxGSSize = enricher_maxGSSize, 
                                          qvalueCutoff =enricher_qvalueCutoff)
  
  all_limmaOnly_signif_enrich_resultDT <- all_limmaOnly_signif_enrich@result
  all_limmaOnly_signif_enrich_resultDT <- all_limmaOnly_signif_enrich_resultDT[order(all_limmaOnly_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
  all_limmaOnly_signif_enrich_resultDT$log10_pval <- -log10(all_limmaOnly_signif_enrich_resultDT[,paste0(padjVarGO)])
  all_limmaOnly_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(all_limmaOnly_signif_enrich_resultDT$GeneRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  all_limmaOnly_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(all_limmaOnly_signif_enrich_resultDT$BgRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  all_limmaOnly_signif_enrich_resultDT$foldEnrichment <- all_limmaOnly_signif_enrich_resultDT$geneRatio/all_limmaOnly_signif_enrich_resultDT$bgRatio
  
  genes_signif_plotMax <- min(c(plotMaxBars, nrow(all_limmaOnly_signif_enrich_resultDT)))
  if(genes_signif_plotMax > 0) {
    for(var_plot in barplot_vars) {
      myTit <- paste0(barplot_vars_tit[var_plot], " - all DS",  " (limma only signif. genes)")
      all_limmaOnly_signif_enrich_resultDT <- all_limmaOnly_signif_enrich_resultDT[order(all_limmaOnly_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
      outFile <- file.path(outFolder,paste0("allDS_all_limmaOnly_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
      do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
      par(oma=c(10,1,1,1))
      barplot(all_limmaOnly_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
              main = myTit, 
              # ylab = paste0("-log10 (", padjVarGO, ")"),
              ylab = paste0(barplot_vars_tit[var_plot]),
              names.arg = gsub("GO_", "", all_limmaOnly_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
              cex.names = 0.6
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
  }
  all_limmaOnly_signif_enrich_resultDT <- data.frame(all_limmaOnly_signif_enrich_resultDT)
  
} else {
  all_limmaOnly_signif_enrich_resultDT <- NULL
}


#***** B) TADs only signif genes
cat(paste0(">  start enricher for TADs only signif genes \n"))

if(length(go_all_tadsOnly_signif_genes) > 0) {
  
  all_tadsOnly_signif_enrich <- enricher(gene = go_all_tadsOnly_signif_genes, 
                                         TERM2GENE=c5_msigdb,
                                         universe = go_all_universe_genes,
                                         pvalueCutoff = enricher_pvalueCutoff, 
                                         pAdjustMethod = enricher_pAdjustMethod, 
                                         minGSSize = enricher_minGSSize, 
                                         maxGSSize = enricher_maxGSSize, 
                                         qvalueCutoff =enricher_qvalueCutoff)
  
  all_tadsOnly_signif_enrich_resultDT <- all_tadsOnly_signif_enrich@result
  all_tadsOnly_signif_enrich_resultDT <- all_tadsOnly_signif_enrich_resultDT[order(all_tadsOnly_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
  all_tadsOnly_signif_enrich_resultDT$log10_pval <- -log10(all_tadsOnly_signif_enrich_resultDT[,paste0(padjVarGO)])
  all_tadsOnly_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(all_tadsOnly_signif_enrich_resultDT$GeneRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  all_tadsOnly_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(all_tadsOnly_signif_enrich_resultDT$BgRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  all_tadsOnly_signif_enrich_resultDT$foldEnrichment <- all_tadsOnly_signif_enrich_resultDT$geneRatio/all_tadsOnly_signif_enrich_resultDT$bgRatio
  
  genes_signif_plotMax <- min(c(plotMaxBars, nrow(all_tadsOnly_signif_enrich_resultDT)))
  if(genes_signif_plotMax > 0) {
    for(var_plot in barplot_vars) {
      myTit <- paste0(barplot_vars_tit[var_plot], " - all DS",  " (TADs only signif. genes)")
      all_tadsOnly_signif_enrich_resultDT <- all_tadsOnly_signif_enrich_resultDT[order(all_tadsOnly_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
      outFile <- file.path(outFolder,paste0("allDS_all_tadsOnly_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
      do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
      par(oma=c(10,1,1,1))
      barplot(all_tadsOnly_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
              main = myTit, 
              # ylab = paste0("-log10 (", padjVarGO, ")"),
              ylab = paste0(barplot_vars_tit[var_plot]),
              names.arg = gsub("GO_", "", all_tadsOnly_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
              cex.names = 0.6
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
  }
  all_tadsOnly_signif_enrich_resultDT <- data.frame(all_tadsOnly_signif_enrich_resultDT)
  
} else {
  all_tadsOnly_signif_enrich_resultDT <- NULL
}


#***** C) all_intersect signif genes
cat(paste0(">  start enricher for all_intersect signif genes \n"))

if(length(go_all_intersect_signif_genes) > 0) {
  
  all_intersect_signif_enrich <- enricher(gene = go_all_intersect_signif_genes, 
                                          TERM2GENE=c5_msigdb,
                                          universe = go_all_universe_genes,
                                          pvalueCutoff = enricher_pvalueCutoff, 
                                          pAdjustMethod = enricher_pAdjustMethod, 
                                          minGSSize = enricher_minGSSize, 
                                          maxGSSize = enricher_maxGSSize, 
                                          qvalueCutoff =enricher_qvalueCutoff)
  
  all_intersect_signif_enrich_resultDT <- all_intersect_signif_enrich@result
  all_intersect_signif_enrich_resultDT <- all_intersect_signif_enrich_resultDT[order(all_intersect_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
  all_intersect_signif_enrich_resultDT$log10_pval <- -log10(all_intersect_signif_enrich_resultDT[,paste0(padjVarGO)])
  all_intersect_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(all_intersect_signif_enrich_resultDT$GeneRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  all_intersect_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(all_intersect_signif_enrich_resultDT$BgRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  all_intersect_signif_enrich_resultDT$foldEnrichment <- all_intersect_signif_enrich_resultDT$geneRatio/all_intersect_signif_enrich_resultDT$bgRatio
  
  genes_signif_plotMax <- min(c(plotMaxBars, nrow(all_intersect_signif_enrich_resultDT)))
  if(genes_signif_plotMax > 0) {
    for(var_plot in barplot_vars) {
      myTit <- paste0(barplot_vars_tit[var_plot], " - all DS",  " (all_intersect signif. genes)")
      all_intersect_signif_enrich_resultDT <- all_intersect_signif_enrich_resultDT[order(all_intersect_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
      outFile <- file.path(outFolder,paste0("allDS_all_intersect_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
      do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
      par(oma=c(10,1,1,1))
      barplot(all_intersect_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
              main = myTit, 
              # ylab = paste0("-log10 (", padjVarGO, ")"),
              ylab = paste0(barplot_vars_tit[var_plot]),
              names.arg = gsub("GO_", "", all_intersect_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
              cex.names = 0.6
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
  }
  all_intersect_signif_enrich_resultDT <- data.frame(all_intersect_signif_enrich_resultDT)
  
} else {
  all_intersect_signif_enrich_resultDT <- NULL
}




#***** 1) all_tads signif genes
cat(paste0(">  start enricher for all_tads signif genes \n"))

if(length(go_all_tads_signif_genes) > 0) {
  
  tad_signif_enrich <- enricher(gene = go_all_tads_signif_genes, 
                                TERM2GENE=c5_msigdb,
                                universe = go_all_universe_genes,
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
      myTit <- paste0(barplot_vars_tit[var_plot], " - all DS",  " (all_tads signif. genes)")
      tad_signif_enrich_resultDT <- tad_signif_enrich_resultDT[order(tad_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
      outFile <- file.path(outFolder,paste0("allDS_all_tads_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
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






#***** 2) all_limma signif. genes

cat(paste0("> start enricher for not_conserved_signif \n"))


if(length(go_all_limma_signif_genes) > 0) {
  
  all_limma_signif_enrich <- enricher(gene = go_all_limma_signif_genes, 
                                      TERM2GENE=c5_msigdb,
                                      universe = go_all_universe_genes,
                                      pvalueCutoff = enricher_pvalueCutoff, 
                                      pAdjustMethod = enricher_pAdjustMethod, 
                                      minGSSize = enricher_minGSSize, 
                                      maxGSSize = enricher_maxGSSize, 
                                      qvalueCutoff =enricher_qvalueCutoff)
  
  all_limma_signif_enrich_resultDT <- all_limma_signif_enrich@result
  all_limma_signif_enrich_resultDT <- all_limma_signif_enrich_resultDT[order(all_limma_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
  all_limma_signif_enrich_resultDT$log10_pval <- -log10(all_limma_signif_enrich_resultDT[,paste0(padjVarGO)])
  all_limma_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(all_limma_signif_enrich_resultDT$GeneRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  all_limma_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(all_limma_signif_enrich_resultDT$BgRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  all_limma_signif_enrich_resultDT$foldEnrichment <- all_limma_signif_enrich_resultDT$geneRatio/all_limma_signif_enrich_resultDT$bgRatio
  
  genes_signif_plotMax <- min(c(plotMaxBars, nrow(all_limma_signif_enrich_resultDT)))
  if(genes_signif_plotMax > 0) {
    for(var_plot in barplot_vars) {
      myTit <- paste0(barplot_vars_tit[var_plot], " - all DS", " (all_limma signif. genes)")
      all_limma_signif_enrich_resultDT <- all_limma_signif_enrich_resultDT[order(all_limma_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
      outFile <- file.path(outFolder,paste0("allDS_all_limma_signif_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
      do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
      par(oma=c(10,1,1,1))
      barplot(all_limma_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
              main = myTit, 
              # ylab = paste0("-log10 (", padjVarGO, ")"),
              ylab = paste0(barplot_vars_tit[var_plot]),
              names.arg = gsub("GO_", "", all_limma_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
              cex.names = 0.6
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
  }
  all_limma_signif_enrich_resultDT <- data.frame(all_limma_signif_enrich_resultDT)
  
  
} else {
  all_limma_signif_enrich_resultDT <- NULL
}


############################################################################################################################################################################ specificity
library(ontologySimilarity)
data(GO_IC)
library(org.Hs.eg.db)
### PREPARE GO DATA
ontologyType <- "BP"

## Bimap interface:
x <- org.Hs.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
all_genes_GO_list <- as.list(x[mapped_genes])

all_genes_GO_list_filter <- lapply( all_genes_GO_list, function(x) {
  Filter(function(k) k[["Ontology"]] == ontologyType, x)
})

all_tads_signif_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% all_tads_signif_genes]
# stopifnot(length(all_tads_genes_GO_list) > 0)
all_tads_signif_genes_GO_list_filter <- lapply(all_tads_signif_genes_GO_list, function(x) {
  Filter(function(k) k[["Ontology"]] == ontologyType, x)
})
all_tads_signif_genes_GO_list_filter_mostSpec <- lapply(all_tads_signif_genes_GO_list_filter, function(x) {
  if(length(x) == 0) return(NA)
  stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
  x[[1]][["GOID"]] # in the list, the 1st is the most specific
})

all_tadsOnly_signif_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% all_tadsOnly_signif_genes]
# stopifnot(length(all_tadsOnly_genes_GO_list) > 0)
all_tadsOnly_signif_genes_GO_list_filter <- lapply(all_tadsOnly_signif_genes_GO_list, function(x) {
  Filter(function(k) k[["Ontology"]] == ontologyType, x)
})
all_tadsOnly_signif_genes_GO_list_filter_mostSpec <- lapply(all_tadsOnly_signif_genes_GO_list_filter, function(x) {
  if(length(x) == 0) return(NA)
  stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
  x[[1]][["GOID"]] # in the list, the 1st is the most specific
})


all_limma_signif_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% all_limma_signif_genes]
# stopifnot(length(all_limma_genes_GO_list) > 0)
all_limma_signif_genes_GO_list_filter <- lapply(all_limma_signif_genes_GO_list, function(x) {
  Filter(function(k) k[["Ontology"]] == ontologyType, x)
})
all_limma_signif_genes_GO_list_filter_mostSpec <- lapply(all_limma_signif_genes_GO_list_filter, function(x) {
  if(length(x) == 0) return(NA)
  stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
  x[[1]][["GOID"]] # in the list, the 1st is the most specific
})


all_limmaOnly_signif_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% all_limmaOnly_signif_genes]
# stopifnot(length(all_limmaOnly_genes_GO_list) > 0)
all_limmaOnly_signif_genes_GO_list_filter <- lapply(all_limmaOnly_signif_genes_GO_list, function(x) {
  Filter(function(k) k[["Ontology"]] == ontologyType, x)
})
all_limmaOnly_signif_genes_GO_list_filter_mostSpec <- lapply(all_limmaOnly_signif_genes_GO_list_filter, function(x) {
  if(length(x) == 0) return(NA)
  stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
  x[[1]][["GOID"]] # in the list, the 1st is the most specific
})


all_intersect_signif_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% all_intersect_signif_genes]
# stopifnot(length(all_intersect_genes_GO_list) > 0)
all_intersect_signif_genes_GO_list_filter <- lapply(all_intersect_signif_genes_GO_list, function(x) {
  Filter(function(k) k[["Ontology"]] == ontologyType, x)
})
all_intersect_signif_genes_GO_list_filter_mostSpec <- lapply(all_intersect_signif_genes_GO_list_filter, function(x) {
  if(length(x) == 0) return(NA)
  stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
  x[[1]][["GOID"]] # in the list, the 1st is the most specific
})



all_intersect_signif_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% all_intersect_signif_genes_GO_list_filter_mostSpec]
all_tads_signif_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% all_tads_signif_genes_GO_list_filter_mostSpec]
all_tadsOnly_signif_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% all_tadsOnly_signif_genes_GO_list_filter_mostSpec]
all_limma_signif_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% all_limma_signif_genes_GO_list_filter_mostSpec]
all_limmaOnly_signif_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% all_limmaOnly_signif_genes_GO_list_filter_mostSpec]

plot_dt <- rbind(
  data.frame(
    signif_type = "tad_signif",
    ic_values = all_tads_signif_genes_GO_list_filter_mostSpec_ic),
  data.frame(
    signif_type = "tadOnly_signif",
    ic_values = all_tadsOnly_signif_genes_GO_list_filter_mostSpec_ic),
  data.frame(
    signif_type = "limma_signif",
    ic_values = all_limma_signif_genes_GO_list_filter_mostSpec_ic),
  data.frame(
    signif_type = "limmaOnly_signif",
    ic_values = all_limmaOnly_signif_genes_GO_list_filter_mostSpec_ic),
  data.frame(
    signif_type = "intersect",
    ic_values = all_intersect_signif_genes_GO_list_filter_mostSpec_ic)
  )


myLevels <- c("tad_signif", "tadOnly_signif", "limma_signif", "limmaOnly_signif", "intersect")
plot_dt$signif_type <- factor(plot_dt$signif_type, levels=myLevels)
stopifnot(!is.na(plot_dt$signif_type))


nLimma <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="limma"]))
nTAD <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="tad"]))
nLimmaOnly <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="limmaOnly"]))
nTADonly <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="tadOnly"]))
nIntersect <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="intersect"]))


subTit <- paste0("# values: TAD=", nTAD,"; limma=", nLimma, "; tadOnly=", nTADonly, "; limmaOnly=", nLimmaOnly, "; intersect=", nIntersect)


p <- ggdensity(plot_dt, 
               title = paste0("signif. gene GO IC"),
               subtitle=paste0(subTit, " (exprds minOcc >=", minOcc, ")"),
               x = "ic_values", 
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               xlab = "signif. genes IC",
               palette = c(col1,col2,col3,col4,col5))

outFile <- file.path(outFolder, paste0("all_gene_ic_values_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))


p <- ggviolin(plot_dt, 
              y = "ic_values", 
              x="signif_type",
              title = paste0("signif. gene GO IC"),
              subtitle=paste0(subTit, " (exprds minOcc >=", minOcc, ")"),
              color = "signif_type", 
              # fill = "signif_type",
              add = "mean", rug = TRUE,
              xlab = "signif. genes IC",
              palette = c(col1,col2,col3,col4,col5))

outFile <- file.path(outFolder, paste0("all_gene_ic_values_boxplot.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myHeightGG)
cat(paste0("... written: ", outFile,"\n"))




######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



