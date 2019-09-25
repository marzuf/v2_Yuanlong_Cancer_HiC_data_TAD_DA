startTime <- Sys.time()
cat(paste0("> Rscript assemble_homer_results.R\n"))

# Rscript assemble_homer_results.R adjPvalComb_0.01 plusminus
# Rscript assemble_homer_results.R signifFDR_0.2 plusminus

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)

buildTable <- TRUE


outFolder <- file.path("ASSEMBLE_HOMER_RESULTS")
dir.create(outFolder, recursive = TRUE)

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))

pipFolder <- file.path("PIPELINE","OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))

setDir="/media/electron"
setDir=""

homerFolder <- file.path("HOMER")  
homerSignifCol <- "q.value..Benjamini."
homerSignifThresh <- 0.05
homerKeepCol <- "Motif.Name"

# bin/findMotifsGenome.pl Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/Barutcu_MCF-10A_40kb_TCGAbrca_lum_bas_adjPvalComb_0.01_plus.txt hg19r Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/adjPvalComb_0.01/MotifOutput_plus -size 200 -p 80

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
signifcol <- "adjPvalComb_0.01"
signifcol <- args[1]
strandselect <- "plusminus"
strandselect <- args[2]


all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

# all_hicds=all_hicds[1]

if(buildTable) {
  
  hicds="GSE105381_HepG2_40kb"
  all_results_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds="TCGAlihc_norm_lihc"  
    ds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat("... start ", hicds, " - ", exprds, "\n")
      
      result_file <- file.path(homerFolder, hicds, exprds, signifcol, paste0("MotifOutput_", strandselect), "knownResults.txt")
      #stopifnot(file.exists(result_file))
      if(!file.exists(result_file)) return(NULL)
      
      result_dt <- read.delim(result_file, header=TRUE, stringsAsFactors = FALSE)
      
      stopifnot(homerKeepCol %in% colnames(result_dt))
      stopifnot(homerSignifCol %in% colnames(result_dt))
      
      out_dt <- result_dt[result_dt[,paste0(homerSignifCol)] <= homerSignifThresh,c(homerKeepCol, homerSignifCol)]
      
      if(nrow(out_dt) == 0) return(NULL)
      
      out_dt$hicds <- hicds
      out_dt$exprds <- exprds
      out_dt[,c("hicds", "exprds", homerKeepCol, homerSignifCol)]
      
      
    } #end-foreach iterating  exprd  
    ds_dt
  } #end-foreach iterating  hicds

  outFile <- file.path(outFolder, paste0("homer_signifMotifs_", signifcol, "_", strandselect, ".Rdata"))
  save(all_results_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else { # end-if buildtable
  outFile <- file.path(outFolder, paste0("homer_signifMotifs_", signifcol, "_", strandselect, ".Rdata"))
  all_results_dt <- get(load(outFile))
}
all_results_dt$dataset <- file.path(all_results_dt$hicds, all_results_dt$exprds)
all_results_dt$motif <- gsub("(^.+?)/.+", "\\1", all_results_dt[,paste0(homerKeepCol)])

all_results_dt$cmpType <-  all_cmps[paste0(all_results_dt$exprds)]
stopifnot(!is.na(all_results_dt$cmpType))

nSignif_all_tmp1 <- aggregate( as.formula(paste0(homerSignifCol, " ~ motif")), FUN=length, data = all_results_dt)
colnames(nSignif_all_tmp1)[colnames(nSignif_all_tmp1) == paste0(homerSignifCol)] <- "nSignif"
nSignif_all_tmp2 <- aggregate( as.formula(paste0("dataset", " ~ motif")), FUN=function(x) paste0(x, collapse=","), data = all_results_dt)
colnames(nSignif_all_tmp2)[colnames(nSignif_all_tmp2) == paste0(homerSignifCol)] <- "datasets"
nSignif_all <- merge(nSignif_all_tmp1, nSignif_all_tmp2, by="motif")
nSignif_all <- nSignif_all[order(nSignif_all$nSignif, decreasing=TRUE),]
outFile <- file.path(outFolder, paste0("nSignif_all_motifCount_", signifcol, "_", strandselect, ".txt"))
write.table(nSignif_all, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

nSignif_byCmpType_tmp1 <- aggregate( as.formula(paste0(homerSignifCol, " ~ motif + cmpType")), FUN=length, data = all_results_dt)
colnames(nSignif_byCmpType_tmp1)[colnames(nSignif_byCmpType_tmp1) == paste0(homerSignifCol)] <- "nSignif"
nSignif_byCmpType_tmp2 <- aggregate( as.formula(paste0("dataset", " ~ motif + cmpType")), FUN=function(x) paste0(x, collapse=","), data = all_results_dt)
colnames(nSignif_byCmpType_tmp2)[colnames(nSignif_byCmpType_tmp2) == paste0(homerSignifCol)] <- "datasets"
nSignif_byCmpType <- merge(nSignif_byCmpType_tmp1, nSignif_byCmpType_tmp2, by=c("motif", "cmpType"))
nSignif_byCmpType <- nSignif_byCmpType[order(nSignif_byCmpType$nSignif, nSignif_byCmpType$cmpType, decreasing = TRUE),]
outFile <- file.path(outFolder, paste0("nSignif_byCmpType_motifCount_", signifcol, "_", strandselect, ".txt"))
write.table(nSignif_byCmpType, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))





