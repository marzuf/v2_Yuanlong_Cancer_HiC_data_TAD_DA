startTime <- Sys.time()
cat(paste0("> Rscript cmp_GO_conserved_notConserved.R\n"))

# Rscript cmp_GO_conserved_notConserved.R 
# Rscript cmp_GO_conserved_notConserved.R norm_vs_tumor
# Rscript cmp_GO_conserved_notConserved.R subtypes
# Rscript cmp_GO_conserved_notConserved.R wt_vs_mut


options(scipen=100)

suppressPackageStartupMessages(library(clusterProfiler, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


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


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  cmpType <- args[1]
} else {
  cmpType <- ""
}

setDir <- "/media/electron"
setDir <- ""


signif_column <- "adjPvalComb"
signifThresh <- 0.01
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3
file_suffix <- paste0(signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes)

outFolder <- file.path("CMP_GO_CONSERVED_NOTCONSERVED", cmpType, file_suffix)
dir.create(outFolder, recursive = TRUE)



inFile <- file.path("GO_SIGNIF_CONSERVED_NOTCONSERVED", cmpType, file_suffix, "all_go_enrich_list.Rdata")
stopifnot(file.exists(inFile))
all_go_enrich_list <- get(load(inFile))

signif_go_pvalThresh <- -log10(0.05)

padjVarGO_plotThresh <- 0.05

all_dt=c(
"conserved_signif_tads_genes_resultDT",
"not_conserved_signif_tads_genes_resultDT"
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
  
}

all_cmp_names <- unique(all_cmps)
cmp=all_cmp_names[1]
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
    
  }
}


######################################################################################
######################################################################################
######################################################################################

cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

