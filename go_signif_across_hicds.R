options(scipen=100)

SSHFS=F

buildData <- TRUE

# Rscript tad_matching_across_hicds_allMatch.R

script_name <- "tad_matching_signif_across_hicds_allMatch.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
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

script0_name <- "0_prepGeneData"

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

outFolder <- file.path("GO_SIGNIF_ACROSS_HICDS")
dir.create(outFolder, recursive = TRUE)

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

# => best TAD matching
# in # of genes
# in bp

final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))

signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)

final_dt[, paste0(signifcol)] <- final_dt[, paste0(signif_column)] <= signifThresh

minOverlapBpRatio <- 0.8
minIntersectGenes <- 3

nRegionLolli <- 10

cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> signifThresh\t=\t", signifThresh, "\n"))
cat(paste0("> minOverlapBpRatio\t=\t", minOverlapBpRatio, "\n"))
cat(paste0("> minIntersectGenes\t=\t", minIntersectGenes, "\n"))
cat(paste0("> nRegionLolli\t=\t", nRegionLolli, "\n"))


inFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH")
stopifnot(dir.exists(inFolder))

inFile <- file.path(inFolder, paste0("conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(inFile))
conserved_signif_tads <- get(load(inFile))


inFile <- file.path(inFolder, paste0("signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(inFile))
signif_tads <- get(load(inFile))

all_conserved_signif_tads <- unique(unlist(conserved_signif_tads))
stopifnot(all_conserved_signif_tads %in% signif_tads)

all_not_conserved_signif_tads <- signif_tads[! signif_tads %in% all_conserved_signif_tads]

####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start preparing data before matching \n")
  
  hicds = all_hicds[1]
  all_gene_list <- foreach(hicds = all_hicds) %dopar% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(hicds_file))
    tadpos_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
    tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]
    stopifnot(nrow(tadpos_dt) > 0 )
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    
    stopifnot(g2t_dt$entrezID %in% gff_dt$entrezID)
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      
      signif_tads <- final_dt$region[final_dt$hicds == hicds & final_dt$exprds == exprds & final_dt[, paste0(signifcol)] ]
      
      gene_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(gene_file))  
      
      region_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
      stopifnot(file.exists(region_file))  
      
      geneList <- get(load(gene_file))
      regionList <- get(load(region_file))
      
      stopifnot(geneList %in% g2t_dt$entrezID)
      ref_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      stopifnot(setequal(ref_g2t_dt$region, regionList))
      
      stopifnot(regionList %in% tadpos_dt$region)
      ref_tadpos_dt <- tadpos_dt[tadpos_dt$region %in% regionList,]
      stopifnot(setequal(ref_g2t_dt$region, ref_tadpos_dt$region))
      
      # take only the data from signif TADs
      stopifnot(setequal(ref_g2t_dt$entrezID, geneList))
      
      
      all_gene_tads <- setNames(file.path(hicds, exprds, ref_g2t_dt$region), ref_g2t_dt$entrezID)
      
      
      conserved_signif_tads_genes <- names(all_gene_tads[all_gene_tads %in% all_conserved_signif_tads]) 
      
      not_conserved_signif_tads_genes <- names(all_gene_tads[all_gene_tads %in% all_not_conserved_signif_tads]) 
      
      stopifnot(length(intersect(conserved_signif_tads_genes, not_conserved_signif_tads_genes)) == 0)
      
      not_signif_tads_genes <- geneList[!geneList %in% c(conserved_signif_tads_genes, not_conserved_signif_tads_genes) ]
      
      
      length(conserved_signif_tads_genes) + length(not_conserved_signif_tads_genes) + length(not_signif_tads_genes) == length(geneList)
      
      
      
      list(
        universe_genes = as.character(geneList),
        conserved_signif_tads_genes = as.character(conserved_signif_tads_genes),
        not_conserved_signif_tads_genes = as.character(not_conserved_signif_tads_genes) #,
        # not_signif_tads_genes = as.character(not_signif_tads_genes)
        
        
      )
    }
    names(exprds_list) <- file.path(hicds, all_exprds[[paste0(hicds)]])
    exprds_list
  }
  all_gene_list <- unlist(all_gene_list, recursive=FALSE)
  stopifnot(length(all_gene_list) == length(all_datasets))
  outFile <- file.path(outFolder, paste0("all_gene_list.Rdata"))
  
  save(all_gene_list, file = file.path(outFolder, "all_gene_list.Rdata"), version=2)
  
} else {
  outFile <- file.path(outFolder, paste0("all_gene_list.Rdata"))
  cat("... load data\n")
  all_gene_list <- get(load(outFile))
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


padjVarGO <- "p.adjust" # p.adjust or qvalue ???

logFile <- file.path(outFolder, "go_signif_across_hicds_logFile.txt")
if(buildData) file.remove(logFile)

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

txt <- paste0("... padjVarGO_plotThresh:\t", padjVarGO_plotThresh, "\n")
printAndLog(txt, logFile)
txt <- paste0("... plotMaxBars:\t", plotMaxBars, "\n")
printAndLog(txt, logFile)



all_universe_genes <- unique(unlist(lapply(all_gene_list, function(x) x[["universe_genes"]])))
all_conserved_signif_tads_genes <- unique(unlist(lapply(all_gene_list, function(x) x[["conserved_signif_tads_genes"]]))) 
all_not_conserved_signif_tads_genes <- unique(unlist(lapply(all_gene_list, function(x) x[["not_conserved_signif_tads_genes"]])))
# all_not_signif_tads_genes 



go_all_universe_genes <- as.character(all_universe_genes)[as.character(all_universe_genes) %in% as.character(c5_msigdb$gene)]

      
go_all_conserved_signif_tads_genes <- as.character(all_conserved_signif_tads_genes)[as.character(all_conserved_signif_tads_genes) %in% as.character(c5_msigdb$gene)]
      
      
      
go_all_not_conserved_signif_tads_genes <- as.character(all_not_conserved_signif_tads_genes)[as.character(all_not_conserved_signif_tads_genes) %in% as.character(c5_msigdb$gene)]
      


# go_all_not_signif_tads_genes <- as.character(all_not_signif_tads_genes)[as.character(all_not_signif_tads_genes) %in% as.character(c5_msigdb$gene)]



cat("... available annot. for all_universe_genes:\t", length(go_all_universe_genes) , "/", length(all_universe_genes), "\n")
cat("... available annot. for all_conserved_signif_tads_genes:\t", length(go_all_conserved_signif_tads_genes) , "/", length(all_conserved_signif_tads_genes), "\n")
cat("... available annot. for all_not_conserved_signif_tads_genes:\t", length(go_all_not_conserved_signif_tads_genes) , "/", length(all_not_conserved_signif_tads_genes), "\n")

# cat("... available annot. for all_not_signif_tads_genes:\t", length(go_all_not_signif_tads_genes) , "/", length(all_not_signif_tads_genes), "\n")

stopifnot(go_all_conserved_signif_tads_genes %in% go_all_universe_genes)
stopifnot(go_all_not_conserved_signif_tads_genes %in% go_all_universe_genes)
# stopifnot(go_all_not_signif_tads_genes %in% go_all_universe_genes)
      
      
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
  
  
  
  
  conserved_signif_enrich_resultDT <- conserved_signif_enrich@result
  conserved_signif_enrich_resultDT <- conserved_signif_enrich_resultDT[order(conserved_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
  conserved_signif_enrich_resultDT$log10_pval <- -log10(conserved_signif_enrich_resultDT[,paste0(padjVarGO)])
  conserved_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(conserved_signif_enrich_resultDT$GeneRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  conserved_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(conserved_signif_enrich_resultDT$BgRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  conserved_signif_enrich_resultDT$foldEnrichment <- conserved_signif_enrich_resultDT$geneRatio/conserved_signif_enrich_resultDT$bgRatio
  
  genes_signif_plotMax <- min(c(plotMaxBars, nrow(conserved_signif_enrich_resultDT)))
  if(genes_signif_plotMax > 0) {
    for(var_plot in barplot_vars) {
      myTit <- paste0(barplot_vars_tit[var_plot], " (conserved_signif)")
      conserved_signif_enrich_resultDT <- conserved_signif_enrich_resultDT[order(conserved_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
      outFile <- file.path(outFolder,paste0("conserved_signif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
      do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
      par(oma=c(10,1,1,1))
      barplot(conserved_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
              main = myTit, 
              # ylab = paste0("-log10 (", padjVarGO, ")"),
              ylab = paste0(barplot_vars_tit[var_plot]),
              names.arg = gsub("GO_", "", conserved_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
              cex.names = 0.6
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
  }
  conserved_signif_enrich_resultDT <- data.frame(conserved_signif_enrich_resultDT)
  # conserved_signif_enrich_resultDT$hicds <- hicds
  # conserved_signif_enrich_resultDT$exprds <- exprds
  
  
} else {
  conserved_signif_enrich_resultDT <- NULL
}



#***** 2) not_conserved_signif

cat(paste0("> start enricher for not_conserved_signif \n"))

if(length(go_all_not_conserved_signif_tads_genes) > 0) {
  
  
  
  not_conserved_signif_enrich <- enricher(gene = go_all_not_conserved_signif_tads_genes, 
                                          TERM2GENE=c5_msigdb,
                                          universe = go_all_universe_genes,
                                          pvalueCutoff = enricher_pvalueCutoff, 
                                          pAdjustMethod = enricher_pAdjustMethod, 
                                          minGSSize = enricher_minGSSize, 
                                          maxGSSize = enricher_maxGSSize, 
                                          qvalueCutoff =enricher_qvalueCutoff)
  
  
  
  
  not_conserved_signif_enrich_resultDT <- not_conserved_signif_enrich@result
  not_conserved_signif_enrich_resultDT <- not_conserved_signif_enrich_resultDT[order(not_conserved_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
  not_conserved_signif_enrich_resultDT$log10_pval <- -log10(not_conserved_signif_enrich_resultDT[,paste0(padjVarGO)])
  not_conserved_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(not_conserved_signif_enrich_resultDT$GeneRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  not_conserved_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(not_conserved_signif_enrich_resultDT$BgRatio, function(x) {
    gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
  not_conserved_signif_enrich_resultDT$foldEnrichment <- not_conserved_signif_enrich_resultDT$geneRatio/not_conserved_signif_enrich_resultDT$bgRatio
  
  genes_signif_plotMax <- min(c(plotMaxBars, nrow(not_conserved_signif_enrich_resultDT)))
  if(genes_signif_plotMax > 0) {
    for(var_plot in barplot_vars) {
      myTit <- paste0(barplot_vars_tit[var_plot], " (not_conserved_signif)")
      not_conserved_signif_enrich_resultDT <- not_conserved_signif_enrich_resultDT[order(not_conserved_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
      outFile <- file.path(outFolder,paste0("not_conserved_signif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
      do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
      par(oma=c(10,1,1,1))
      barplot(not_conserved_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
              main = myTit, 
              # ylab = paste0("-log10 (", padjVarGO, ")"),
              ylab = paste0(barplot_vars_tit[var_plot]),
              names.arg = gsub("GO_", "", not_conserved_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
              cex.names = 0.6
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
  }
  not_conserved_signif_enrich_resultDT <- data.frame(not_conserved_signif_enrich_resultDT)
  # not_conserved_signif_enrich_resultDT$hicds <- hicds
  # not_conserved_signif_enrich_resultDT$exprds <- exprds
  
  
} else {
  not_conserved_signif_enrich_resultDT <- NULL
}

#       
#       #***** 1) not_signif
#       
#       cat(paste0("... ", hicds, " - ", exprds, " - start enricher for genes_signifTADs \n"))
#       
#       if(length(go_genes_signifTADs) > 0) {
#         
#         
#         
#         genes_signifTADs_enrich <- enricher(gene = go_genes_signifTADs, 
#                                             TERM2GENE=c5_msigdb,
#                                             universe = universe_genes,
#                                             pvalueCutoff = enricher_pvalueCutoff, 
#                                             pAdjustMethod = enricher_pAdjustMethod, 
#                                             minGSSize = enricher_minGSSize, 
#                                             maxGSSize = enricher_maxGSSize, 
#                                             qvalueCutoff =enricher_qvalueCutoff)
#         
#         
#         
#         
#         genes_signifTADs_enrich_resultDT <- genes_signifTADs_enrich@result
#         genes_signifTADs_enrich_resultDT <- genes_signifTADs_enrich_resultDT[order(genes_signifTADs_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
#         genes_signifTADs_enrich_resultDT$log10_pval <- -log10(genes_signifTADs_enrich_resultDT[,paste0(padjVarGO)])
#         genes_signifTADs_enrich_resultDT$geneRatio <- as.numeric(sapply(genes_signifTADs_enrich_resultDT$GeneRatio, function(x) {
#           gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
#         genes_signifTADs_enrich_resultDT$bgRatio <- as.numeric(sapply(genes_signifTADs_enrich_resultDT$BgRatio, function(x) {
#           gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
#         genes_signifTADs_enrich_resultDT$foldEnrichment <- genes_signifTADs_enrich_resultDT$geneRatio/genes_signifTADs_enrich_resultDT$bgRatio
#         
#         genes_signif_plotMax <- min(c(plotMaxBars, nrow(genes_signifTADs_enrich_resultDT)))
#         if(genes_signif_plotMax > 0) {
#           for(var_plot in barplot_vars) {
#             myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds, " (genes_signifTADs)")
#             genes_signifTADs_enrich_resultDT <- genes_signifTADs_enrich_resultDT[order(genes_signifTADs_enrich_resultDT[,var_plot], decreasing=TRUE),]
#             outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_genes_signifTADs", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
#             do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
#             par(oma=c(10,1,1,1))
#             barplot(genes_signifTADs_enrich_resultDT[1:genes_signif_plotMax,var_plot],
#                     main = myTit, 
#                     # ylab = paste0("-log10 (", padjVarGO, ")"),
#                     ylab = paste0(barplot_vars_tit[var_plot]),
#                     names.arg = gsub("GO_", "", genes_signifTADs_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
#                     cex.names = 0.6
#             )
#             foo <- dev.off()
#             cat(paste0("... written: ", outFile, "\n"))
#           }
#         }
#         genes_signifTADs_enrich_resultDT <- data.frame(genes_signifTADs_enrich_resultDT)
#         genes_signifTADs_enrich_resultDT$hicds <- hicds
#         genes_signifTADs_enrich_resultDT$exprds <- exprds
#         
#       } else {
#         genes_signifTADs_enrich_resultDT <- NULL
#       }
#       
#       
#       #***** 2) genes_signif
#       
#       cat(paste0("... ", hicds, " - ", exprds, " - start enricher for genes_signif \n"))
#       
#       
#       if(length(go_genes_signif) > 0){
#         genes_signif_enrich <- enricher(gene = go_genes_signif, 
#                                         TERM2GENE=c5_msigdb,
#                                         universe = universe_genes,
#                                         pvalueCutoff = enricher_pvalueCutoff, 
#                                         pAdjustMethod = enricher_pAdjustMethod, 
#                                         minGSSize = enricher_minGSSize, 
#                                         maxGSSize = enricher_maxGSSize, 
#                                         qvalueCutoff =enricher_qvalueCutoff)
#         
#         genes_signif_enrich_resultDT <- genes_signif_enrich@result
#         genes_signif_enrich_resultDT <- genes_signif_enrich_resultDT[order(genes_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
#         genes_signif_enrich_resultDT$log10_pval <- -log10(genes_signif_enrich_resultDT[,paste0(padjVarGO)])
#         genes_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(genes_signif_enrich_resultDT$GeneRatio, function(x) {
#           gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
#         genes_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(genes_signif_enrich_resultDT$BgRatio, function(x) {
#           gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
#         genes_signif_enrich_resultDT$foldEnrichment <- genes_signif_enrich_resultDT$geneRatio/genes_signif_enrich_resultDT$bgRatio
#         
#         genes_signif_plotMax <- min(c(plotMaxBars, nrow(genes_signif_enrich_resultDT)))
#         if(genes_signif_plotMax > 0) {
#           for(var_plot in barplot_vars) {
#             myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds, " (genes_signif)")
#             genes_signif_enrich_resultDT <- genes_signif_enrich_resultDT[order(genes_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
#             outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_genes_signif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
#             do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
#             par(oma=c(10,1,1,1))
#             barplot(genes_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
#                     main = myTit, 
#                     # ylab = paste0("-log10 (", padjVarGO, ")"),
#                     ylab = paste0(barplot_vars_tit[var_plot]),
#                     names.arg = gsub("GO_", "", genes_signif_enrich_resultDT$Description[1:genes_signif_plotMax]), las=2,
#                     cex.names = 0.6
#             )
#             foo <- dev.off()
#             cat(paste0("... written: ", outFile, "\n"))
#           }
#         }
#         genes_signif_enrich_resultDT <- data.frame(genes_signif_enrich_resultDT)
#         genes_signif_enrich_resultDT$hicds <- hicds
#         genes_signif_enrich_resultDT$exprds <- exprds
#         save(genes_signif_enrich_resultDT, file="genes_signif_enrich_resultDT.Rdata", version=2)
#         
#       } else {
#         genes_signif_enrich_resultDT <- NULL 
#       }
#       
#       
#       
#       list(
#         genes_signifTADs_enrich_resultDT=genes_signifTADs_enrich_resultDT,
#         not_conserved_signif_enrich_resultDT=not_conserved_signif_enrich_resultDT,
#         conserved_signif_enrich_resultDT= conserved_signif_enrich_resultDT,
#         genes_signif_enrich_resultDT= genes_signif_enrich_resultDT
#       )
#       
#       
#       
#       
#     } # end-for iterating over hicds
#     names(exprds_data) <- all_exprds[[paste0(hicds)]]
#     exprds_data
#     
#   } # end-foreach iterating over exprcds
#   
#   
#   names(all_enricher_result) <- all_hicds
#   
#   outFile <- file.path(outFolder, "all_enricher_result.Rdata")
#   save(all_enricher_result, file = outFile, version=2)
#   cat(paste0("... written: ", outFile, "\n"))
# } else {
#   outFile <- file.path(outFolder, "all_enricher_result.Rdata")  
#   all_enricher_result <- get(load(outFile))
# }
# 
# 
# signif_go_pvalThresh <- -log10(0.05)
# 
# all_dt=c(
#   "genes_signifTADs_enrich_resultDT",
#   "not_conserved_signif_enrich_resultDT",
#   "conserved_signif_enrich_resultDT",
#   "genes_signif_enrich_resultDT")
# 
# hicds=all_hicds[1]
# dt = all_dt[1]
# 
# topCommonBars <- 10
# 
# for(dt in all_dt) {
#   all_go_categories <- foreach(hicds = all_hicds) %dopar% {
#     exprds = all_exprds[[paste0(hicds)]][1]
#     exprds_go_categories <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
#       curr_dt <- all_enricher_result[[paste0(hicds)]][[paste0(exprds)]][[paste0(dt)]]
#       rownames(curr_dt)[curr_dt[,paste0(padjVarGO)] <= padjVarGO_plotThresh]
#     } 
#     names(exprds_go_categories) <- all_exprds[[paste0(hicds)]]
#     exprds_go_categories
#   }
#   go_categories <- unlist(all_go_categories)
#   go_categories_count <- setNames(as.numeric(table(go_categories)), names(table(go_categories)))
#   go_categories_count <- sort(go_categories_count, decreasing = TRUE)
#   outFile <- file.path(outFolder,paste0("all_ds_", dt, "_intersect_",padjVarGO, "_barplot", ".", plotType))
#   do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
#   par(oma=c(10,1,1,1))
#   barplot(go_categories_count[1:topCommonBars], las=2, 
#           ylab="# of datasets",
#           names.arg=gsub("GO_", "", names(go_categories_count[1:topCommonBars])),
#           main=paste0(gsub("_resultDT", "", dt), " - intersect across DS"),
#           cex.names=0.6
#   )
#   
#   mtext(side=3, text=paste0(padjVarGO, "<=",padjVarGO_plotThresh, " (top ", topCommonBars, ")"))
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))
#   
# }
# 
# all_cmp_names <- unique(all_cmps)
# cmp=all_cmp_names[1]
# for(cmp in all_cmp_names){
#   
#   for(dt in all_dt) {
#     all_go_categories <- foreach(hicds = all_hicds) %dopar% {
#       exprds = all_exprds[[paste0(hicds)]][1]
#       exprds_go_categories <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
#         curr_dt <- all_enricher_result[[paste0(hicds)]][[paste0(exprds)]][[paste0(dt)]]
#         cmpType <- all_cmps[paste0(exprds)]
#         stopifnot(!is.na(cmpType))
#         if(cmpType != cmp) return(NULL)
#         
#         rownames(curr_dt)[curr_dt[,paste0(padjVarGO)] <= padjVarGO_plotThresh]
#       } 
#       names(exprds_go_categories) <- all_exprds[[paste0(hicds)]]
#       exprds_go_categories
#     }
#     go_categories <- unlist(all_go_categories)
#     go_categories_count <- setNames(as.numeric(table(go_categories)), names(table(go_categories)))
#     go_categories_count <- sort(go_categories_count, decreasing = TRUE)
#     outFile <- file.path(outFolder,paste0("all_ds_",cmp, "_", dt, "_intersect_",padjVarGO, "_barplot", ".", plotType))
#     do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
#     par(oma=c(10,1,1,1))
#     barplot(go_categories_count[1:topCommonBars], las=2, 
#             names.arg=gsub("GO_", "", names(go_categories_count[1:topCommonBars])),
#             ylab="# of datasets",
#             main=paste0(gsub("_resultDT", "", dt), " - ", cmp),
#             cex.names=0.6
#     )
#     mtext(side=3, text=paste0(padjVarGO, "<=",padjVarGO_plotThresh, " (top ", topCommonBars, ")"))
#     foo <- dev.off()
#     cat(paste0("... written: ", outFile, "\n"))
#     
#   }
# }
# 

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


  
  