startTime <- Sys.time()
cat(paste0("> Rscript signif_genes_specifity_intersectDiff.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
library(ontologySimilarity)
data(GO_IC)
require(reshape2)

require(ggpubr)

buildTable <- TRUE

printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}

options(scipen=100)

# Rscript signif_genes_specifity_intersectDiff.R 0.01 0.05 # <tadpval> <genepval>

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4


myHeightGG <- 7
myWidthGG <- myHeightGG*1.2

limmaCol <- "#00AFBB"
tadCol <- "#FC4E07"

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

mainFolder <- file.path(".")

pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"

all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

args <- commandArgs(trailingOnly = TRUE)
tads_signifThresh <- 0.01
tads_signifThresh <- args[1]

genes_signifTresh <- 0.05
genes_signifTresh <- args[2]

file_suffix <- paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifTresh)

outFolder <- file.path("SIGNIF_GENES_SPECIFICITY_INTERSECTDIFF", file_suffix)
dir.create(outFolder, recursive = TRUE)

logFile <- file.path(outFolder, "signif_genes_specifity_intersectDiff_logFile.txt")
file.remove(logFile)

final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))
signif_column <- "adjPvalComb"
signifcol <- paste0(signif_column, "_", tads_signifThresh)
final_dt[, paste0(signifcol)] <- final_dt[, paste0(signif_column)] <= tads_signifThresh
cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> tads_signifThresh\t=\t", tads_signifThresh, "\n"))


all_genes_GO_list_filter <- lapply( all_genes_GO_list, function(x) {
  Filter(function(k) k[["Ontology"]] == ontologyType, x)
})
all_genes_GO_list_filter_ul <- as.character(unlist(lapply(unlist(all_genes_GO_list_filter, recursive=FALSE), function(x)x[["GOID"]])))

GO_freq <- setNames(
  as.numeric(table(all_genes_GO_list_filter_ul)),
  names(table(all_genes_GO_list_filter_ul))
)
GO_freq <- GO_freq/length(all_genes_GO_list_filter_ul)


####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start preparing data before matching \n")
  
  hicds = all_hicds[1]
  all_signif_genes_ic <- foreach(hicds = all_hicds) %dopar% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
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
      
      tadsOnly_signif_genes <- setdiff(tads_signif_genes, limma_signif_genes)
      stopifnot(tadsOnly_signif_genes %in% tads_signif_genes)
      
      limmaOnly_signif_genes <- setdiff(limma_signif_genes, tads_signif_genes)
      stopifnot(limmaOnly_signif_genes %in% limma_signif_genes)

      intersect_signif_genes <- intersect(limma_signif_genes, tads_signif_genes)
      stopifnot(intersect_signif_genes %in% limma_signif_genes)
      stopifnot(intersect_signif_genes %in% tads_signif_genes)
      
      
      ####### intersect_genes GO
      
      
      signif_intersect_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% intersect_signif_genes]
      stopifnot(length(signif_intersect_genes_GO_list) > 0)
      signif_intersect_genes_GO_list_filter <- lapply(signif_intersect_genes_GO_list, function(x) {
        Filter(function(k) k[["Ontology"]] == ontologyType, x)
      })
      signif_intersect_genes_GO_list_filter_mostSpec <- lapply(signif_intersect_genes_GO_list_filter, function(x) {
        if(length(x) == 0) return(NA)
        stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
        x[[1]][["GOID"]] # in the list, the 1st is the most specific
      })
      signif_intersect_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% signif_intersect_genes_GO_list_filter_mostSpec]
      txt <- paste0(hicds, " - ", exprds, " - av. IC for signif. intersect genes:\t", length(signif_intersect_genes_GO_list_filter_mostSpec_ic), "/", length(signif_intersect_genes_GO_list_filter_mostSpec), "\n")
      printAndLog(txt, logFile)
      
      
      signif_intersect_genes_GO_list_filter_mostSpec_freq <- GO_freq[names(GO_freq) %in% signif_intersect_genes_GO_list_filter_mostSpec]
      
      
      intersect_genes_GO <- names(GO_IC)[names(GO_IC) %in% signif_intersect_genes_GO_list_filter_mostSpec]
      sigif_intersect_nGO <- length(intersect_genes_GO)
      sigif_intersect_nUniqueGO <- length(unique(intersect_genes_GO))
      
    
      
      
      ####### limmaOnly_genes GO
      
      
      signif_limmaOnly_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% limmaOnly_signif_genes]
      stopifnot(length(signif_limmaOnly_genes_GO_list) > 0)
      signif_limmaOnly_genes_GO_list_filter <- lapply(signif_limmaOnly_genes_GO_list, function(x) {
        Filter(function(k) k[["Ontology"]] == ontologyType, x)
      })
      signif_limmaOnly_genes_GO_list_filter_mostSpec <- lapply(signif_limmaOnly_genes_GO_list_filter, function(x) {
        if(length(x) == 0) return(NA)
        stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
        x[[1]][["GOID"]] # in the list, the 1st is the most specific
      })
      signif_limmaOnly_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% signif_limmaOnly_genes_GO_list_filter_mostSpec]
      txt <- paste0(hicds, " - ", exprds, " - av. IC for signif. limma only genes:\t", length(signif_limmaOnly_genes_GO_list_filter_mostSpec_ic), "/", length(signif_limmaOnly_genes_GO_list_filter_mostSpec), "\n")
      printAndLog(txt, logFile)
      
      
      signif_limmaOnly_genes_GO_list_filter_mostSpec_freq <- GO_freq[names(GO_freq) %in% signif_limmaOnly_genes_GO_list_filter_mostSpec]
      
      
      limmaOnly_genes_GO <- names(GO_IC)[names(GO_IC) %in% signif_limmaOnly_genes_GO_list_filter_mostSpec]
      sigif_limmaOnly_nGO <- length(limmaOnly_genes_GO)
      sigif_limmaOnly_nUniqueGO <- length(unique(limmaOnly_genes_GO))
      
      
      
      
      
      ####### tadOnly_genes GO
      
      
      signif_tadOnly_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% tadsOnly_signif_genes]
      stopifnot(length(signif_tadOnly_genes_GO_list) > 0)
      signif_tadOnly_genes_GO_list_filter <- lapply(signif_tadOnly_genes_GO_list, function(x) {
        Filter(function(k) k[["Ontology"]] == ontologyType, x)
      })
      signif_tadOnly_genes_GO_list_filter_mostSpec <- lapply(signif_tadOnly_genes_GO_list_filter, function(x) {
        if(length(x) == 0) return(NA)
        stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
        x[[1]][["GOID"]] # in the list, the 1st is the most specific
      })
      signif_tadOnly_genes_GO_list_filter_mostSpec_ic <- GO_IC[names(GO_IC) %in% signif_tadOnly_genes_GO_list_filter_mostSpec]
      txt <- paste0(hicds, " - ", exprds, " - av. IC for signif. TADs only genes:\t", length(signif_tadOnly_genes_GO_list_filter_mostSpec_ic), "/", length(signif_tadOnly_genes_GO_list_filter_mostSpec), "\n")
      printAndLog(txt, logFile)
      
      
      signif_tadOnly_genes_GO_list_filter_mostSpec_freq <- GO_freq[names(GO_freq) %in% signif_tadOnly_genes_GO_list_filter_mostSpec]
      
      
      tadOnly_genes_GO <- names(GO_IC)[names(GO_IC) %in% signif_tadOnly_genes_GO_list_filter_mostSpec]
      sigif_tadOnly_nGO <- length(tadOnly_genes_GO)
      sigif_tadOnly_nUniqueGO <- length(unique(tadOnly_genes_GO))
      
      
      
      ####### tad_genes GO
      
      
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

            
      signif_tad_genes_GO_list_filter_mostSpec_freq <- GO_freq[names(GO_freq) %in% signif_tad_genes_GO_list_filter_mostSpec]
      
      
      tad_genes_GO <- names(GO_IC)[names(GO_IC) %in% signif_tad_genes_GO_list_filter_mostSpec]
      sigif_tad_nGO <- length(tad_genes_GO)
      sigif_tad_nUniqueGO <- length(unique(tad_genes_GO))
      
      ####### limma_genes GO
      
      
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
      
      signif_limma_genes_GO_list_filter_mostSpec_freq <- GO_freq[names(GO_freq) %in% signif_limma_genes_GO_list_filter_mostSpec]
      
      
      limma_genes_GO <- names(GO_IC)[names(GO_IC) %in% signif_limma_genes_GO_list_filter_mostSpec]
      sigif_limma_nGO <- length(limma_genes_GO)
      sigif_limma_nUniqueGO <- length(unique(limma_genes_GO))
      
      
      list( ic_signif_tad_genes=signif_tad_genes_GO_list_filter_mostSpec_ic,
            ic_signif_tadOnly_genes=signif_tadOnly_genes_GO_list_filter_mostSpec_ic,
            ic_signif_limma_genes=signif_limma_genes_GO_list_filter_mostSpec_ic,
            ic_signif_limmaOnly_genes=signif_limmaOnly_genes_GO_list_filter_mostSpec_ic,
            ic_signif_intersect_genes=signif_intersect_genes_GO_list_filter_mostSpec_ic,
            
            freq_signif_tad_genes=signif_tad_genes_GO_list_filter_mostSpec_freq,
            freq_signif_tadOnly_genes=signif_tadOnly_genes_GO_list_filter_mostSpec_freq,
            freq_signif_limma_genes=signif_limma_genes_GO_list_filter_mostSpec_freq,
            freq_signif_limmaOnly_genes=signif_limmaOnly_genes_GO_list_filter_mostSpec_freq,
            freq_signif_intersect_genes=signif_intersect_genes_GO_list_filter_mostSpec_freq
            
            )
      
    }# end-foreach iterating over exprds
    
    exprds_dt
    
  }# end-foreach iterating over hicds

  outFile <- file.path(outFolder, paste0("all_signif_genes_ic.Rdata"))
  save(all_signif_genes_ic, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
    outFile <- file.path(outFolder, paste0( "all_signif_genes_ic.Rdata"))
    cat("... load data\n")
    all_signif_genes_ic <- get(load(outFile))
}

all_signif_genes_ic_ul <- unlist(all_signif_genes_ic, recursive = FALSE)

#####################################################################################################################################

plot_dt <- rbind(
  data.frame(
    signif_type="limma",
  ic_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["ic_signif_limma_genes"]])),
  stringsAsFactors=FALSE),
  data.frame(
    signif_type="tad",
    ic_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["ic_signif_tad_genes"]])),
    stringsAsFactors=FALSE),
  data.frame(
    signif_type="tadOnly",
    ic_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["ic_signif_tadOnly_genes"]])),
    stringsAsFactors=FALSE),
  data.frame(
    signif_type="limmaOnly",
    ic_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["ic_signif_limmaOnly_genes"]])),
    stringsAsFactors=FALSE),
  data.frame(
    signif_type="intersect",
    ic_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["ic_signif_intersect_genes"]])),
    stringsAsFactors=FALSE)
  )


nLimma <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="limma"]))
nTAD <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="tad"]))
nLimmaOnly <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="limmaOnly"]))
nTADonly <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="tadOnly"]))
nIntersect <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="intersect"]))


subTit <- paste0("# values: TAD=", nTAD,"; limma=", nLimma, "; tadOnly=", nTADonly, "; limmaOnly=", nLimmaOnly, "; intersect=", nIntersect)


p <- ggdensity(plot_dt, 
               title = paste0("signif. gene GO IC"),
               subtitle=paste0(subTit),
               x = "ic_values", 
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               xlab = "signif. genes IC",
               palette = c(limmaCol, tadCol))

outFile <- file.path(outFolder, paste0("all_gene_ic_values_", file_suffix, "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))


p <- ggviolin(plot_dt, 
               y = "ic_values", 
               x="signif_type",
              title = paste0("signif. gene GO IC"),
              subtitle=paste0(subTit),
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               xlab = "signif. genes IC",
               palette = c(limmaCol, tadCol))

outFile <- file.path(outFolder, paste0("all_gene_ic_values_", file_suffix, "_boxplot.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myHeightGG)
cat(paste0("... written: ", outFile,"\n"))


#####################################################################################################################################

plot_dt <- rbind(
  data.frame(
    signif_type="limma",
    freq_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["freq_signif_limma_genes"]])),
    stringsAsFactors=FALSE),
  data.frame(
    signif_type="tad",
    freq_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["freq_signif_tad_genes"]])),
    stringsAsFactors=FALSE)    ,
  data.frame(
    signif_type="limmaOnly",
    freq_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["freq_signif_limmaOnly_genes"]])),
    stringsAsFactors=FALSE)    ,
  data.frame(
    signif_type="tadOnly",
    freq_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["freq_signif_tadOnly_genes"]])),
    stringsAsFactors=FALSE)    ,
  data.frame(
    signif_type="intersect",
    freq_values=unlist(lapply(all_signif_genes_ic_ul, function(x) x[["freq_signif_intersect_genes"]])),
    stringsAsFactors=FALSE)    
)


nLimma <- sum(!is.na(plot_dt$freq_values[plot_dt$signif_type=="limma"]))
nTAD <- sum(!is.na(plot_dt$freq_values[plot_dt$signif_type=="tad"]))
nLimmaOnly <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="limmaOnly"]))
nTADonly <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="tadOnly"]))
nIntersect <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="intersect"]))

subTit <- paste0("# values: TAD=", nTAD,"; limma=", nLimma, "; tadOnly=", nTADonly, "; limmaOnly=", nLimmaOnly, "; intersect=", nIntersect)


p <- ggdensity(plot_dt, 
               title = paste0("signif. gene GO freq"),
               subtitle=paste0(subTit)
               x = "freq_values", 
               color = "signif_type", fill = "signif_type",
               add = "mean", rug = TRUE,
               xlab = "signif. genes freq",
               palette = c(limmaCol, tadCol))

outFile <- file.path(outFolder, paste0("all_gene_freq_values_", file_suffix, "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))


p <- ggviolin(plot_dt, 
              y = "freq_values", 
              x="signif_type",
              title = paste0("signif. gene GO freq"),
              subtitle=paste0(subTit),
              color = "signif_type", fill = "signif_type",
              add = "mean", rug = TRUE,
              xlab = "signif. genes freq",
              palette = c(limmaCol, tadCol))

outFile <- file.path(outFolder, paste0("all_gene_freq_values_", file_suffix, "_boxplot.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myHeightGG)
cat(paste0("... written: ", outFile,"\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




    
      
      