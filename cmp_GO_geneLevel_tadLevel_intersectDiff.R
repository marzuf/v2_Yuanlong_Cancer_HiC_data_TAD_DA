startTime <- Sys.time()
cat(paste0("> Rscript cmp_GO_geneLevel_tadLevel_intersectDiff.R\n"))

# Rscript cmp_GO_geneLevel_tadLevel_intersectDiff.R 0.01 0.05

options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


outFolder <- file.path("CMP_GO_GENELEVEL_TADLEVEL_INTERSECTDIFF")
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


padjVarGO <- "p.adjust" # p.adjust or qvalue ???



TAD_pvalThresh <- 0.01
gene_pvalThresh <- 0.05
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
TAD_pvalThresh <- args[1]
gene_pvalThresh <- args[2]


setDir <- "/media/electron"
setDir <- ""


inFile <- file.path("GO_SIGNIF_GENELEVEL_TADLEVEL_INTERSECTDIFF", paste0("tadPvalThresh", TAD_pvalThresh, "_genePvalThresh", gene_pvalThresh), "all_go_enrich_list.Rdata")
all_go_enrich_list <- get(load(inFile))


signif_go_pvalThresh <- -log10(0.05)

padjVarGO_plotThresh <- 0.05

all_dt=c(
"tadsOnly_signif_enrich_resultDT",
"limmaOnly_signif_enrich_resultDT",
"tad_signif_enrich_resultDT",
"limma_signif_enrich_resultDT",
"intersect_signif_enrich_resultDT"
)

dt = all_dt[1]

topCommonBars <- 10

curr_dataset = names(all_go_enrich_list)[1]

for(dt in all_dt) {
  all_go_categories <- foreach(curr_dataset = names(all_go_enrich_list)) %dopar% {
      curr_dt <- all_go_enrich_list[[paste0(curr_dataset)]][[paste0(dt)]]
    rownames(curr_dt)[curr_dt[,paste0(padjVarGO)] <= padjVarGO_plotThresh]
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
  
  
  outFile <- file.path(outFolder,paste0("all_ds_", dt, "_intersect_",padjVarGO, "_textTable", ".", "txt"))
  write.table(data.frame(go=names(go_categories_count),count=go_categories_count,stringsAsFactors = FALSE), 
              file = outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=F, append=F)
  cat(paste0("... written: ", outFile, "\n"))
  
}

all_cmp_names <- unique(all_cmps)
cmp=all_cmp_names[2]
for(cmp in all_cmp_names){
  
  toKeepDS <- names(all_go_enrich_list)[basename(names(all_go_enrich_list)) %in% names(all_cmps)[all_cmps==cmp]]
  
  for(dt in all_dt) {
      all_go_categories <- foreach(curr_dataset = toKeepDS) %dopar% {
        curr_dt <- all_go_enrich_list[[paste0(curr_dataset)]][[paste0(dt)]]
        rownames(curr_dt)[curr_dt[,paste0(padjVarGO)] <= padjVarGO_plotThresh]
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
    
    
    outFile <- file.path(outFolder,paste0("all_ds_",cmp, "_", dt, "_intersect_",padjVarGO, "_textTable", ".", "txt"))
    write.table(data.frame(go=names(go_categories_count),count=go_categories_count,stringsAsFactors = FALSE), 
                file = outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=F, append=F)
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
}


######################################################################################
######################################################################################
######################################################################################

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

