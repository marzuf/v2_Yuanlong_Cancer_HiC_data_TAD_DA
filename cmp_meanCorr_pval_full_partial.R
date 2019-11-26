startTime <- Sys.time()
cat(paste0("> Rscript cmp_meanCorr_pval_full_partial.R\n"))

script_name <- "cmp_meanCorr_pval_full_partial.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript cmp_meanCorr_pval_full_partial.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- TRUE

rankTiesMet <- "min"

myHeight <- 400 
myWidth <- 400
plotType <- "png"
plotCex <- 1.4

outFolder <- file.path(paste0("CMP_MEANCORR_PVAL_FULL_PARTIAL"))
dir.create(outFolder, recursive = TRUE)

meanCorrPartial_pval_pattern <- file.path("10sameNbrPartial_runEmpPvalMeanTADCorr", "emp_pval_meanCorr.Rdata")
combPartial_pval_pattern  <- file.path("11sameNbrPartial_runEmpPvalCombined", "emp_pval_combined.Rdata")

meanCorr_pval_pattern <- file.path("10sameNbr_runEmpPvalMeanTADCorr", "emp_pval_meanCorr.Rdata")
comb_pval_pattern <- file.path("11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(pipFolder)
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

if(buildTable) {
  all_ds_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
  
  meanCorr_pval <- get(load(file.path(pipFolder, hicds, exprds, meanCorr_pval_pattern)))
  meanCorr_pval <- p.adjust(meanCorr_pval, method="BH")
  
  meanCorrPartial_pval <- get(load(file.path(pipFolder, hicds, exprds, meanCorrPartial_pval_pattern)))
  meanCorrPartial_pval <- p.adjust(meanCorrPartial_pval, method="BH")
  
  stopifnot(setequal(names(meanCorr_pval), names(meanCorrPartial_pval)))
  
  comb_pval <- get(load(file.path(pipFolder, hicds, exprds, comb_pval_pattern)))
  comb_pval <- p.adjust(comb_pval, method="BH")
  combPartial_pval <- get(load(file.path(pipFolder, hicds, exprds, combPartial_pval_pattern)))
  combPartial_pval <- p.adjust(combPartial_pval, method="BH")
  
  stopifnot(setequal(names(comb_pval), names(combPartial_pval)))
  
  stopifnot(setequal(names(meanCorr_pval), names(comb_pval)))
  
  all_regs <- names(comb_pval)
  
  
  data.frame(
    hicds=hicds,
    exprds=exprds,
    reg=all_regs,
    meanCorr_pval=meanCorr_pval[all_regs],
    meanCorr_pval_rank=rank(meanCorr_pval[all_regs], ties.method=rankTiesMet),
    meanCorrPartial_pval=meanCorrPartial_pval[all_regs],
    meanCorrPartial_pval_rank=rank(meanCorrPartial_pval[all_regs], ties.method=rankTiesMet),
    comb_pval=comb_pval[all_regs],
    comb_pval_rank=rank(comb_pval[all_regs], ties.method=rankTiesMet),
    combPartial_pval=combPartial_pval[all_regs],
    combPartial_pval_rank=rank(combPartial_pval[all_regs], ties.method=rankTiesMet),
    stringsAsFactors=FALSE
  )
  
  }
  outFile <- file.path(outFolder, "all_ds_dt.Rdata")
  save(all_ds_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {
  outFile <- file.path(outFolder, "all_ds_dt.Rdata")
  # outFile <- file.path(outFolder, "sub_coexpr_and_purity_dt.Rdata")
  cat(paste0("... load coexpr data - ", Sys.time()))
  all_ds_dt <- get(load(outFile))
  cat(paste0(" - ", Sys.time(), " - done\n"))
}



source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_x <- c("meanCorr_pval","meanCorr_pval_rank", "comb_pval", "comb_pval_rank")
all_y <- c("meanCorrPartial_pval", "meanCorrPartial_pval_rank", "combPartial_pval", "combPartial_pval_rank")

i=1
for(i in 1:length(all_x)){
  
  myx <- all_x[i]
  myy <- all_y[i]
  
  all_ds_dt[,paste0(myx, "_log10")] <- -log10(all_ds_dt[,paste0(myx)])
  all_ds_dt[,paste0(myy, "_log10")] <- -log10(all_ds_dt[,paste0(myy)])
  
  outFile <- file.path(outFolder, paste0("all_datasets_", myy, "_", myx, "_log10.", plotType))
  do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
  densplot(
    x = all_ds_dt[,paste0(myx, "_log10")] ,
    y = all_ds_dt[,paste0(myy, "_log10")] ,
    xlab = paste0(myx, " [-log10]"),
    ylab = paste0(myy, " [-log10]"),
    main=paste0("partial vs. full pvals"),
    cex.axis=plotCex,
    cex.lab=plotCex,
    pch=16,
    cex=0.7
  )
  cor_est <- cor(all_ds_dt[,paste0(myx, "_log10")], all_ds_dt[,paste0(myy, "_log10")], method="pearson")
  mtext(side=3, text = paste0("(adj. p-val.)"))
  # legend("topleft", legend = paste0("n=", nrow(all_ds_dt)), bty="n")
  legend("topleft", legend = c(paste0("n=", nrow(all_ds_dt)), paste0("PCC=", round(cor_est, 4))), bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))    
  
  outFile <- file.path(outFolder, paste0("all_datasets_", myy, "_", myx, ".", plotType))
  do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
  densplot(
    x = all_ds_dt[,paste0(myx)] ,
    y = all_ds_dt[,paste0(myy)] ,
    xlab = paste0(myx, ""),
    ylab = paste0(myy, ""),
    main=paste0("partial vs. full pvals"),
    cex.axis=plotCex,
    cex.lab=plotCex,
    pch=16,
    cex=0.7
  )
  cor_est <- cor(all_ds_dt[,paste0(myx)], all_ds_dt[,paste0(myy)], method="pearson")
  mtext(side=3, text = paste0("(adj. p-val.)"))
  legend("topleft", legend = c(paste0("n=", nrow(all_ds_dt)), paste0("PCC=", round(cor_est, 4))), bty="n")
  

    
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))    
  
  
}



##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
