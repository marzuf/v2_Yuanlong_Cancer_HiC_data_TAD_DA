startTime <- Sys.time()
cat(paste0("> Rscript assemble_homer_results.R\n"))

# Rscript assemble_homer_results.R

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)

buildTable <- FALSE


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
      stopifnot(file.exists(result_file))
      
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

all_results_dt$motif <- gsub("(^.+?)/.+", "\\1", all_results_dt[,paste0(homerKeepCol)])

all_results_dt$cmpType <-  all_cmps[paste0(all_results_dt$exprds)]
stopifnot(!is.na(all_results_dt$cmpType))



nSignif_all <- aggregate( as.formula(paste0(homerSignifCol, " ~ motif")), FUN=length, data = all_results_dt)
nSignif_byCmpType <- aggregate( as.formula(paste0(homerSignifCol, " ~ motif + cmpType")), FUN=length, data = all_results_dt)






