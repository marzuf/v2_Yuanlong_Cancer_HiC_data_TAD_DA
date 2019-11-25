########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript immune_genes_tads_features.R\n"))

script_name <- "immune_genes_tads_features.R"


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript immune_genes_tads_features.R

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4
myHeightGG <- 7
myWidthGG <- 7

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))



outFolder <- file.path("IMMUNE_GENES_TADS_FEATURES")
dir.create(outFolder, recursive = TRUE)

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

immune_genes <- c("CD", "C1QA", "GIMAP", "CCR", "CCRL")
immune_genes <- c("C1QA", "GIMAP", "CCR", "CCRL")

matching_symbols <- sapply(final_dt$region_genes, function(x) {
  tad_genes <- unlist(strsplit(x=x, split=","))
  # if any of the immune gene match any of the tad genes
  any(unlist(sapply(immune_genes, function(patt) any(grepl(paste0(patt, "\\d+$"), tad_genes)))))
})

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")


all_hicds <- list.files(pipOutFolder)
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))



script10_name <- "10sameNbr_runEmpPvalMeanTADCorr"
script9_name <- "9_runEmpPvalMeanTADLogFC"

#all_ds=all_ds[1]
save(all_ds, file="all_ds.Rdata", version=2)
all_pval_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
  # adjust !
  
  pvalfc_file <- file.path(pipOutFolder, hicds, exprds, script9_name, "emp_pval_meanLogFC.Rdata")
  stopifnot(file.exists(pvalfc_file))
  pvalfc <- get(load(pvalfc_file))
  adj_pvalfc <- p.adjust(pvalfc, method="BH")
# cat("A\n")
  
  pvalcorr_file <- file.path(pipOutFolder, hicds, exprds, script10_name, "emp_pval_meanCorr.Rdata")
  stopifnot(file.exists(pvalcorr_file))
  pvalcorr <- get(load(pvalcorr_file))
  adj_pvalcorr <- p.adjust(pvalcorr, method="BH")
  # cat("B\n")
  
  stopifnot(setequal(names(adj_pvalcorr), names(adj_pvalfc)))
  
  all_regs <- names(adj_pvalcorr)
  
  outdt <- data.frame(
    hicds =  hicds,
    exprds=exprds,
    region = all_regs,
    adjPval_meanLogFC = adj_pvalfc[all_regs],
    adjPval_meanCorr = adj_pvalcorr[all_regs],
    stringsAsFactors = FALSE
  )
  stopifnot(!is.na(outdt))
  outdt
} # end-foreach iterating ds

save(final_dt, file="final_dt.Rdata", version=2)
save(all_pval_dt, file="all_pval_dt.Rdata", version=2)
final_allpval_dt <- merge(final_dt, all_pval_dt, by=c("hicds", "exprds", "region"), all.x=TRUE, all.y=FALSE)
stopifnot(nrow(final_allpval_dt) == nrow(final_dt))
stopifnot(!is.na(final_allpval_dt))


# outFile <- file.path(outFolder, paste0("meanLogFC_meanCorr_selectedImmune.", plotType))
# do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
# par(bty="l")
# plot(
#   meanLogFC ~ meanCorr, 
#   data =final_allpval_dt,
#   pch=16,
#   cex=0.7,
#   cex.lab=plotCex,
#   cex.axis=plotCex,
#   col="grey"
# )
# points(x = final_allpval_dt$meanCorr[matching_symbols],
#        y = final_allpval_dt$meanLogFC[matching_symbols],
#        cex =0.7,
#        pch=16,
#      col="red"
#        )
# 
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))    

final_allpval_dt$adjPval_meanLogFC_log10 <- -log10(final_allpval_dt$adjPval_meanLogFC)
final_allpval_dt$adjPval_meanCorr_log10 <- -log10(final_allpval_dt$adjPval_meanCorr)


all_y <- c("adjPval_meanLogFC", "adjPval_meanLogFC_log10", "meanLogFC")
all_x <- c("adjPval_meanCorr", "adjPval_meanCorr_log10", "meanCorr")


for(i in 1:length(all_x)) {
  
  myx <- all_x[i]
  myy <- all_y[i]
  
  
  outFile <- file.path(outFolder, paste0(paste0(immune_genes, collapse="_"), "_", myy, "_vs_", myx,"_selectedImmune.", plotType))
  do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
  
  par(bty="l")
  
  plot(
    as.formula(paste0(myy, " ~ ", myx)),
    data =final_allpval_dt,
    pch=16,
    main = paste0("all TADs, all datasets"),
    cex=0.7,
    col="grey",
    cex.lab=plotCex,
    cex.axis=plotCex
  )
  points(x = final_allpval_dt[, paste0(myx)][matching_symbols],
         y = final_allpval_dt[,paste0(myy)][matching_symbols],
         cex =0.7,
         pch=16,
         col="red"
  )
  mtext(side=3, text = paste0("TADs with genes matching: ", paste0(immune_genes, collapse = ",")))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))    
  
}

# outFile <- file.path(outFolder, paste0("adjPvalMeanLogFC_adjPvalMeanCorr_selectedImmune.", plotType))
# do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
# 
# par(bty="l")
# 
# plot(
#   adjPval_meanLogFC ~ adjPval_meanCorr, 
#   data =final_allpval_dt,
#   pch=16,
#   cex=0.7,
#   col="grey",
#   cex.lab=plotCex,
#   cex.axis=plotCex
# )
# points(x = final_allpval_dt$adjPval_meanCorr[matching_symbols],
#        y = final_allpval_dt$adjPval_meanLogFC[matching_symbols],
#        cex =0.7,
#        pch=16,
#        col="red"
# )
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))    
# 
# 
# 
# final_allpval_dt$adjPval_meanLogFC_log10 <- -log10(final_allpval_dt$adjPval_meanLogFC)
# final_allpval_dt$adjPval_meanCorr_log10 <- -log10(final_allpval_dt$adjPval_meanCorr)
# 
# outFile <- file.path(outFolder, paste0("adjPvalMeanLogFC_log10_adjPvalMeanCorr_log10_selectedImmune.", plotType))
# do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
# 
# par(bty="l")
# 
# plot(
#   adjPval_meanLogFC_log10 ~ adjPval_meanCorr_log10, 
#   data =final_allpval_dt,
#   pch=16,
#   cex=0.7,
#   col="grey",
#   cex.lab=plotCex,
#   cex.axis=plotCex
# )
# points(x = final_allpval_dt$adjPval_meanCorr[matching_symbols],
#        y = final_allpval_dt$adjPval_meanLogFC[matching_symbols],
#        cex =0.7,
#        pch=16,
#        col="red"
# )
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))    
# 
# 
# 
# 
# 


