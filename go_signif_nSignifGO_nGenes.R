script_name <- "go_signif_nSignifGO_nGenes.R"
options(scipen=100)

SSHFS=F

cat("> START ", script_name, "\n")

# Rscript go_signif_nSignifGO_nGenes.R 0.01 0.05


library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)
library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)
library(ggpubr)

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))
startTime <- Sys.time()


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightGG <- 7
myWidthGG <- myHeightGG*1.2


col1 <- get_palette("Dark2", 5)[1]
col2 <- get_palette("Dark2", 5)[2]
col3 <- get_palette("Dark2", 5)[3]
col4 <- get_palette("Dark2", 5)[4]
col5 <- get_palette("Dark2", 5)[5]

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


buildTable <- TRUE

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

args <- commandArgs(trailingOnly = TRUE)

tads_signifThresh <- args[1]
genes_signifTresh <- args[2]

outFolder <- file.path("GO_SIGNIF_NSIGNIFGO_NGENES", paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifTresh))
dir.create(outFolder, recursive = TRUE)

logFile <- file.path(outFolder, "go_signif_nSignifGO_nGenes_logFile.txt")
if(buildTable) file.remove(logFile)

printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}


go_signif_col <- "p.adjust"
go_signifThresh <- 0.05


txt <- paste0("... go_signif_col\t=\t", go_signif_col, "\n")
printAndLog(txt, logFile)
txt <- paste0("... go_signifThresh\t=\t", go_signifThresh, "\n")
printAndLog(txt, logFile)

file_suffix <- paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifTresh)

inFile <- file.path("GO_SIGNIF_GENELEVEL_TADLEVEL_INTERSECTDIFF", paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifTresh), "all_go_enrich_list.Rdata")
stopifnot(file.exists(inFile))
all_go_enrich_list <- get(load(inFile))

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

if(buildTable) {
  
  hicds = all_hicds[1]
  all_nGO_nGenes_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat(paste0("... start: ", hicds, " - ", exprds, "\n"))
      
      ############# intersect
      go_signif_intersect_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["intersect_signif_enrich_resultDT"]]
      if(!is.null(go_signif_intersect_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - intersect signif.: # annot. GO:\t", nrow(go_signif_intersect_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_intersect_dt <- go_signif_intersect_dt[go_signif_intersect_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - intersect signif.: # signif. annot. GO:\t", nrow(go_signif_intersect_dt), "\n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_intersect <- nrow(go_signif_intersect_dt)
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... intersect signif.: NULL \n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_intersect <- NA
      }
      ############# limma only
      go_signif_limmaOnly_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["limmaOnly_signif_enrich_resultDT"]]
      if(!is.null(go_signif_limmaOnly_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - limma only signif.: # annot. GO:\t", nrow(go_signif_limmaOnly_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_limmaOnly_dt <- go_signif_limmaOnly_dt[go_signif_limmaOnly_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - limma only signif.: # signif. annot. GO:\t", nrow(go_signif_limmaOnly_dt), "\n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_limmaOnly <- nrow(go_signif_limmaOnly_dt)
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... limmaOnly signif.: NULL \n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_limmaOnly <- NA
      }
      ############# TAD only
      go_signif_tadsOnly_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["tadsOnly_signif_enrich_resultDT"]]
      if(!is.null(go_signif_tadsOnly_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - TAD only signif: # annot. GO:\t", nrow(go_signif_tadsOnly_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_tadsOnly_dt <- go_signif_tadsOnly_dt[go_signif_tadsOnly_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - TAD only signif: # signif. annot. GO:\t", nrow(go_signif_tadsOnly_dt), "\n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_tadsOnly <- nrow(go_signif_tadsOnly_dt)
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... tadOnly signif.: NULL \n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_tadsOnly <- NA
      }
      ############# limma
      go_signif_limma_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["limma_signif_enrich_resultDT"]]
      if(!is.null(go_signif_limma_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - limma signif.: # annot. GO:\t", nrow(go_signif_limma_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_limma_dt <- go_signif_limma_dt[go_signif_limma_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - limma signif.: # signif. annot. GO:\t", nrow(go_signif_limma_dt), "\n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_limma <- nrow(go_signif_limma_dt)
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: NULL \n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_limma <- NA
      }
      ############# TAD
      go_signif_tads_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["tad_signif_enrich_resultDT"]]
      if(!is.null(go_signif_tads_dt)) {
        txt <- paste0(hicds, " - ", exprds, " - TAD signif: # annot. GO:\t", nrow(go_signif_tads_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_tads_dt <- go_signif_tads_dt[go_signif_tads_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - TAD signif: # signif. annot. GO:\t", nrow(go_signif_tads_dt), "\n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_tads <- nrow(go_signif_tads_dt)
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... TAD signif.: NULL \n")
        printAndLog(txt, logFile)
        nEnrichedGO_signif_tads <- NA
      }
      nGenes_tad_signif <- length(all_go_enrich_list[[file.path(hicds, exprds)]][["tads_signif_genes"]])
      nGenes_tadOnly_signif <- length(all_go_enrich_list[[file.path(hicds, exprds)]][["tadsOnly_signif_genes"]])
      nGenes_limma_signif <- length(all_go_enrich_list[[file.path(hicds, exprds)]][["limma_signif_genes"]])
      nGenes_limmaOnly_signif <- length(all_go_enrich_list[[file.path(hicds, exprds)]][["limmaOnly_signif_genes"]])
      nGenes_intersect_signif <- length(intersect(all_go_enrich_list[[file.path(hicds, exprds)]][["tads_signif_genes"]], all_go_enrich_list[[file.path(hicds, exprds)]][["limma_signif_genes"]]))
      
      data.frame(
        hicds=hicds,
        exprds=exprds,
        nEnrichedGO_signif_intersect=nEnrichedGO_signif_intersect,
        nEnrichedGO_signif_limmaOnly=nEnrichedGO_signif_limmaOnly,
        nEnrichedGO_signif_tadOnly=nEnrichedGO_signif_tadsOnly,
        nEnrichedGO_signif_limma=nEnrichedGO_signif_limma,
        nEnrichedGO_signif_tad=nEnrichedGO_signif_tads,
        nGenes_tad_signif=nGenes_tad_signif,
        nGenes_tadOnly_signif=nGenes_tadOnly_signif,
        nGenes_limma_signif=nGenes_limma_signif,
        nGenes_limmaOnly_signif=nGenes_limmaOnly_signif,
        nGenes_intersect_signif=nGenes_intersect_signif,
        stringsAsFactors = FALSE
      )
    } # end-foreach iterating over exprds
    exprds_dt
  } # end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, paste0("all_nGO_nGenes_dt.Rdata"))
  save(all_nGO_nGenes_dt, file = outFile, version=2)
  
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0( "all_nGO_nGenes_dt.Rdata"))
  cat("... load data\n")
  all_nGO_nGenes_dt <- get(load(outFile))
}


stopifnot(all_nGO_nGenes_dt$exprds %in% names(all_cmps))
all_nGO_nGenes_dt$cmptype <- all_cmps[all_nGO_nGenes_dt$exprds]
stopifnot(!is.na(all_nGO_nGenes_dt$cmptype))
all_nGO_nGenes_dt$cmptype_col <- all_cols[all_nGO_nGenes_dt$cmptype]
stopifnot(!is.na(all_nGO_nGenes_dt$cmptype_col))

dotCols <- all_nGO_nGenes_dt$cmptype_col

all_types <- c("intersect", "tad", "tadOnly", "limma", "limmaOnly")
signif_type = all_types[1]
for(signif_type in all_types) {
  
  y_var <- paste0("nEnrichedGO_signif_", signif_type)
  x_var <- paste0("nGenes_", signif_type, "_signif")
  stopifnot(x_var %in% colnames(all_nGO_nGenes_dt))
  stopifnot(y_var %in% colnames(all_nGO_nGenes_dt))
  
  myx <- all_nGO_nGenes_dt[,paste0(x_var)]
  myy <- all_nGO_nGenes_dt[,paste0(y_var)]
  
  
  outFile <- file.path(outFolder, paste0("nGO_vs_nSignifGenes_", signif_type, "_densplot.", plotType ))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = myx,
    y = myy,
    main=paste0("# GO vs. # genes - ", signif_type),
    sub=paste0(""),
    xlab=paste0(x_var),
    ylab=paste0(y_var),
    col = dotCols,
    cex=0.9,
    pch=16,
    cex.axis=axisCex,
    cex.lab=axisCex
  )
  addSubtypeLeg(mypos="topright",bty="n")
  addCorr(x = myx, y=myy, bty="n", legPos = "topleft")
  mtext(side=3, paste0("n=", length(!is.na(myx))))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}



######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


