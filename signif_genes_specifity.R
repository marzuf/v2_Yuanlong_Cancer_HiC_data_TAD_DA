startTime <- Sys.time()
cat(paste0("> Rscript signif_genes_specifity.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
library(ontologySimilarity)
data(GO_IC)
require(reshape2)

buildTable <- TRUE

printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}

options(scipen=100)

# Rscript signif_genes_specifity.R 

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

### PREPARE GO DATA
ontologyType <- "BP"

## Bimap interface:
x <- org.Hs.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
all_genes_GO_list <- as.list(x[mapped_genes])

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

outFolder <- file.path("SIGNIF_GENES_SPECIFICITY")
dir.create(outFolder, recursive = TRUE)

logFile <- file.path(outFolder, "signif_genes_specifity_logFile.txt")
file.remove(logFile)

mainFolder <- file.path(".")

pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"

all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

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

####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start preparing data before matching \n")
  
  hicds = all_hicds[1]
  all_signif_genes_meanIC_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat(paste0("... start ", hicds, " - ", exprds, "\n"))
      
      dataset_pipDir <- file.path(pipFolder, hicds, exprds)
      stopifnot(dir.exists(dataset_pipDir))
      
      
      
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
      
      
      signif_tad_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% tads_signif_genes]
      stopifnot(length(signif_tad_genes_GO_list) > 0)
      signif_tad_genes_GO_list_filter <- lapply(signif_tad_genes_GO_list, function(x) {
        Filter(function(k) k[["Ontology"]] == ontologyType, x)
      })
      signif_tad_genes_GO_list_filter_mostSpec <- lapply(signif_tad_genes_GO_list_filter, function(x) {
        if(length(x) == 0) return(NA)
        stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
        x[[1]][["GOID"]] # in the list, the 1st is the most specific
      })
      signif_tad_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% signif_tad_genes_GO_list_filter_mostSpec]
      txt <- paste0(hicds, " - ", exprds, " - av. IC for signif. TADs genes:\t", length(signif_tad_genes_GO_list_filter_mostSpec_ic), "/", length(signif_tad_genes_GO_list_filter_mostSpec), "\n")
      printAndLog(txt, logFile)
      mean_ic_signif_tad_genes <- mean(signif_tad_genes_GO_list_filter_mostSpec_ic)
      
      tad_genes_GO <- names(GO_IC)[names(GO_IC) %in% signif_tad_genes_GO_list_filter_mostSpec]
      sigif_tad_nGO <- length(tad_genes_GO)
      sigif_tad_nUniqueGO <- length(unique(tad_genes_GO))
      
      
      signif_limma_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% limma_signif_genes]
      stopifnot(length(signif_limma_genes_GO_list) > 0)
      signif_limma_genes_GO_list_filter <- lapply(signif_limma_genes_GO_list, function(x) {
        Filter(function(k) k[["Ontology"]] == ontologyType, x)
      })
      signif_limma_genes_GO_list_filter_mostSpec <- lapply(signif_limma_genes_GO_list_filter, function(x) {
        if(length(x) == 0) return(NA)
        stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
        x[[1]][["GOID"]] # in the list, the 1st is the most specific
      })
      signif_limma_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% signif_limma_genes_GO_list_filter_mostSpec]
      txt <- paste0(hicds, " - ", exprds, " - av. IC for signif. TADs genes:\t", length(signif_limma_genes_GO_list_filter_mostSpec_ic), "/", length(signif_limma_genes_GO_list_filter_mostSpec), "\n")
      printAndLog(txt, logFile)
      mean_ic_signif_limma_genes <- mean(signif_limma_genes_GO_list_filter_mostSpec_ic)
      
      
      limma_genes_GO <- names(GO_IC)[names(GO_IC) %in% signif_limma_genes_GO_list_filter_mostSpec]
      sigif_limma_nGO <- length(limma_genes_GO)
      sigif_limma_nUniqueGO <- length(unique(limma_genes_GO))
      
      
      data.frame(hicds=hicds,
                 exprds=exprds,
                 sigif_tad_nGO=sigif_tad_nGO,
                 sigif_limma_nGO=sigif_limma_nGO,
                 sigif_tad_nUniqueGO=sigif_tad_nUniqueGO,
                 sigif_limma_nUniqueGO=sigif_limma_nUniqueGO,
                 mean_ic_signif_tad_genes=mean_ic_signif_tad_genes,
                 mean_ic_signif_limma_genes=mean_ic_signif_limma_genes,
                 stringsAsFactors = FALSE)
      
    }# end-foreach iterating over exprds
    
    exprds_dt
    
  }# end-foreach iterating over hicds

  outFile <- file.path(outFolder, paste0("all_signif_genes_meanIC_dt.Rdata"))
  save(all_signif_genes_meanIC_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0( "all_signif_genes_meanIC_dt.Rdata"))
  cat("... load data\n")
  all_signif_genes_meanIC_dt <- get(load(outFile))
}


all_vars_to_plot <- c("mean_ic", "nGO", "nUniqueGO")

plot_dt <- melt(all_signif_genes_meanIC_dt, id=c("hicds", "exprds"))

nDS <- length(unique(paste0(all_signif_genes_meanIC_dt$hicds, all_signif_genes_meanIC_dt$exprds)))

var_to_plot="mean_ic"
for(var_to_plot in all_vars_to_plot) {
  
  var_plot_dt <- plot_dt[grepl(var_to_plot, plot_dt$variable),]
  var_plot_dt$variable <- as.character(var_plot_dt$variable)
  
  var_plot_dt$variable <- gsub(var_to_plot, "", var_plot_dt$variable)
  
  outFile <- file.path(outFolder,paste0(var_to_plot, "_signifTAD_signifLimma_genes_boxplots.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))    
  boxplot(value~variable, data=var_plot_dt, las=1, main = var_to_plot)  
  mtext(side=3, text=paste0("(nDS = ", nDS, ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

all(all_signif_genes_meanIC_dt$sigif_limma_nUniqueGO == all_signif_genes_meanIC_dt$sigif_limma_nGO)
# [1] TRUE
all(all_signif_genes_meanIC_dt$sigif_tad_nUniqueGO == all_signif_genes_meanIC_dt$sigif_tad_nGO)
# [1] TRUE


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




    
      
      