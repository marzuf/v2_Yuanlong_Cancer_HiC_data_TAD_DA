options(scipen=100)

SSHFS=F

buildData <- TRUE

# Rscript go_signif_conserved_notConserved_MF.R
# Rscript go_signif_conserved_notConserved_MF.R norm_vs_tumor
# Rscript go_signif_conserved_notConserved_MF.R subtypes
# Rscript go_signif_conserved_notConserved_MF.R wt_vs_mut

script_name <- "go_signif_conserved_notConserved_MF.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(clusterProfiler)
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


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  data_cmpType <- ""
  file_prefix <- ""
} else if(length(args) == 1) {
  data_cmpType <- args[1]  
  stopifnot(data_cmpType %in% c("norm_vs_tumor", "subtypes", "wt_vs_mut"))
  file_prefix <- paste0(data_cmpType, "_")
} else {
  stop("error\n") 
}


tads_signifThresh <- 0.01
genes_signifTresh <- 0.05


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

outFolder <- file.path("GO_SIGNIF_CONSERVED_NOTCONSERVED_MF", data_cmpType, paste0(signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes))
dir.create(outFolder, recursive = TRUE)


inFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", data_cmpType)
stopifnot(dir.exists(inFolder))

inFile <- file.path(inFolder, paste0(file_prefix, "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(inFile))
conserved_signif_tads <- get(load(inFile))

inFile <- file.path(inFolder, paste0(file_prefix, "signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(inFile))
signif_tads <- get(load(inFile))


all_conserved_signif_tads <- unique(unlist(conserved_signif_tads))
stopifnot(all_conserved_signif_tads %in% signif_tads)

all_not_conserved_signif_tads <- signif_tads[! signif_tads %in% all_conserved_signif_tads]


padjVarGO <- "p.adjust" # p.adjust or qvalue ???

logFile <- file.path(outFolder, "go_signif_conserved_notConserved_logFile.txt")
if(buildData) file.remove(logFile)

barplot_vars <- c("foldEnrichment", "geneRatio", "log10_pval")
barplot_vars_tit <- setNames(c("Fold enrichment", "Gene ratio", paste0("-log10(", padjVarGO,  ")")), barplot_vars)

plotMaxBars <- 10

padjVarGO_plotThresh <- 0.05



# GO for BP nad MF [do not take c5_CC]
enricher_ontologyType="MF"
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








####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start preparing data before matching \n")
  
  hicds = all_hicds[1]
  all_go_enrich_list <- foreach(hicds = all_hicds) %dopar% {
    
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
      
      # list(
      #   universe_genes = as.character(geneList),
      #   conserved_signif_tads_genes = as.character(conserved_signif_tads_genes),
      #   not_conserved_signif_tads_genes = as.character(not_conserved_signif_tads_genes) #,
      # )
      universe_genes <- as.character(geneList)
      go_universe_genes <- as.character(universe_genes)[as.character(universe_genes) %in% as.character(c5_msigdb$gene)]
      go_not_conserved_signif_tads_genes <- as.character(not_conserved_signif_tads_genes)[as.character(not_conserved_signif_tads_genes) %in% as.character(c5_msigdb$gene)]
      go_conserved_signif_tads_genes <- as.character(conserved_signif_tads_genes)[as.character(conserved_signif_tads_genes) %in% as.character(c5_msigdb$gene)]
      
      txt <- paste0("... ", hicds, " - ", exprds, " available annot. for universe_genes:\t", length(go_universe_genes) , "/", length(universe_genes), "\n")
      txt <- paste0("... ", hicds, " - ", exprds, " available annot. for not_conserved_signif_tads_genes:\t", length(go_not_conserved_signif_tads_genes) , "/", length(not_conserved_signif_tads_genes), "\n")
      txt <- paste0("... ", hicds, " - ", exprds, " available annot. for conserved_signif_tads_genes:\t", length(go_conserved_signif_tads_genes) , "/", length(conserved_signif_tads_genes), "\n")
      
      stopifnot(go_not_conserved_signif_tads_genes %in% go_universe_genes)
      stopifnot(go_conserved_signif_tads_genes %in% go_universe_genes)
      
      #***** 1) not_conserved_signif_tads_genes
      
      cat(paste0("> start enricher for not_conserved_signif \n"))
      
      
      if(length(go_not_conserved_signif_tads_genes) > 0) {
        
        not_conserved_signif_tads_genes <- enricher(gene = go_not_conserved_signif_tads_genes, 
                                                    TERM2GENE=c5_msigdb,
                                                    universe = go_universe_genes,
                                                    pvalueCutoff = enricher_pvalueCutoff, 
                                                    pAdjustMethod = enricher_pAdjustMethod, 
                                                    minGSSize = enricher_minGSSize, 
                                                    maxGSSize = enricher_maxGSSize, 
                                                    qvalueCutoff =enricher_qvalueCutoff)
        
        not_conserved_signif_tads_genes_resultDT <- not_conserved_signif_tads_genes@result
        not_conserved_signif_tads_genes_resultDT <- not_conserved_signif_tads_genes_resultDT[order(not_conserved_signif_tads_genes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        not_conserved_signif_tads_genes_resultDT$log10_pval <- -log10(not_conserved_signif_tads_genes_resultDT[,paste0(padjVarGO)])
        not_conserved_signif_tads_genes_resultDT$geneRatio <- as.numeric(sapply(not_conserved_signif_tads_genes_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        not_conserved_signif_tads_genes_resultDT$bgRatio <- as.numeric(sapply(not_conserved_signif_tads_genes_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        not_conserved_signif_tads_genes_resultDT$foldEnrichment <- not_conserved_signif_tads_genes_resultDT$geneRatio/not_conserved_signif_tads_genes_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(not_conserved_signif_tads_genes_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds, " (not conserved signif. TADs genes)")
            not_conserved_signif_tads_genes_resultDT <- not_conserved_signif_tads_genes_resultDT[order(not_conserved_signif_tads_genes_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_not_conserved_signif_TADs_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(not_conserved_signif_tads_genes_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", not_conserved_signif_tads_genes_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        not_conserved_signif_tads_genes_resultDT <- data.frame(not_conserved_signif_tads_genes_resultDT)
        # not_conserved_signif_tads_genes_resultDT$hicds <- hicds
        # not_conserved_signif_tads_genes_resultDT$exprds <- exprds
        
        
      } else {
        not_conserved_signif_tads_genes_resultDT <- NULL
      }
      
      
      #***** 2) conserved_signif_tads_genes
      
      cat(paste0("> start enricher for not_conserved_signif \n"))
      
      
      if(length(go_conserved_signif_tads_genes) > 0) {
        
        conserved_signif_tads_genes <- enricher(gene = go_conserved_signif_tads_genes, 
                                        TERM2GENE=c5_msigdb,
                                        universe = go_universe_genes,
                                        pvalueCutoff = enricher_pvalueCutoff, 
                                        pAdjustMethod = enricher_pAdjustMethod, 
                                        minGSSize = enricher_minGSSize, 
                                        maxGSSize = enricher_maxGSSize, 
                                        qvalueCutoff =enricher_qvalueCutoff)
        
        conserved_signif_tads_genes_resultDT <- conserved_signif_tads_genes@result
        conserved_signif_tads_genes_resultDT <- conserved_signif_tads_genes_resultDT[order(conserved_signif_tads_genes_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
        conserved_signif_tads_genes_resultDT$log10_pval <- -log10(conserved_signif_tads_genes_resultDT[,paste0(padjVarGO)])
        conserved_signif_tads_genes_resultDT$geneRatio <- as.numeric(sapply(conserved_signif_tads_genes_resultDT$GeneRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        conserved_signif_tads_genes_resultDT$bgRatio <- as.numeric(sapply(conserved_signif_tads_genes_resultDT$BgRatio, function(x) {
          gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
        conserved_signif_tads_genes_resultDT$foldEnrichment <- conserved_signif_tads_genes_resultDT$geneRatio/conserved_signif_tads_genes_resultDT$bgRatio
        
        genes_signif_plotMax <- min(c(plotMaxBars, nrow(conserved_signif_tads_genes_resultDT)))
        if(genes_signif_plotMax > 0) {
          for(var_plot in barplot_vars) {
            myTit <- paste0(barplot_vars_tit[var_plot], " - ", hicds, " - ", exprds, " (conserved signif. TADs genes)")
            conserved_signif_tads_genes_resultDT <- conserved_signif_tads_genes_resultDT[order(conserved_signif_tads_genes_resultDT[,var_plot], decreasing=TRUE),]
            outFile <- file.path(outFolder,paste0(hicds, "_", exprds, "_conserved_signif_TADs_genes_", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
            do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.2))    
            par(oma=c(10,1,1,1))
            barplot(conserved_signif_tads_genes_resultDT[1:genes_signif_plotMax,var_plot],
                    main = myTit, 
                    # ylab = paste0("-log10 (", padjVarGO, ")"),
                    ylab = paste0(barplot_vars_tit[var_plot]),
                    names.arg = gsub("GO_", "", conserved_signif_tads_genes_resultDT$Description[1:genes_signif_plotMax]), las=2,
                    cex.names = 0.6
            )
            foo <- dev.off()
            cat(paste0("... written: ", outFile, "\n"))
          }
        }
        conserved_signif_tads_genes_resultDT <- data.frame(conserved_signif_tads_genes_resultDT)
        # conserved_signif_tads_genes_resultDT$hicds <- hicds
        # conserved_signif_tads_genes_resultDT$exprds <- exprds
        
        
      } else {
        conserved_signif_tads_genes_resultDT <- NULL
      }
      
      
      
      
      
      
      
      
      
      
      list(
        conserved_signif_tads_genes_resultDT=conserved_signif_tads_genes_resultDT,
        not_conserved_signif_tads_genes_resultDT=not_conserved_signif_tads_genes_resultDT
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
# outFile="GO_SIGNIF_GENELEVEL_TADLEVEL/tadPvalThresh0.05_genePvalThresh0.05/all_go_enrich_list.Rdata"
# all_go_enrich_list <- get(load(outFile))

######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

