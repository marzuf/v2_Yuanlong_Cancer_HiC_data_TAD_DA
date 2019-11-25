########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript gene_TADpval_corrExprPurity.R\n"))

script_name <- "gene_TADpval_corrExprPurity.R"

# Rscript  gene_TADpval_corrExprPurity.R
# Rscript  gene_TADpval_corrExprPurity.R EPIC

myHeight <- 400
myWidth <- 400
plotType <- "png"
plotCex <- 1.4


args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  purity_ds <- "" 
  file_suffix <- ""
} else {
  stopifnot(length(args) == 1)
  purity_ds <- args[1]
  file_suffix <- paste0("_", purity_ds)
}


outFolder <- file.path(paste0("GENE_TADPVAL_CORREXPRPURITY", file_suffix))
dir.create(outFolder, recursive = TRUE)

corr_expr_purity_dt <- get(load(file.path(paste0("CORR_EXPR_AND_PURITY", file_suffix),"corr_expr_purity_dt.Rdata")))
gene_tad_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

all_dt <- merge(corr_expr_purity_dt, gene_tad_dt, by=c("hicds", "exprds", "entrezID"))
all_dt$tad_adjCombPval_log10 <- -log10(all_dt$tad_adjCombPval)

myxlab <- "TAD adj. comb. pval [log10]"
myylab <- "corr. gene expression and purity"

outFile <- file.path(outFolder, paste0("gene_corrExprPurity_tadPval.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  all_corr_gene_purity ~ tad_adjCombPval_log10,
  data=all_dt,
  xlab=myxlab,
  ylab=myylab,
  pch=16,
  cex=0.7,
  cex.axis=plotCex,
  cex.lab=plotCex
)
mtext(side=3, text = paste0(purity_ds))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))    


outFile <- file.path(outFolder, paste0("gene_corrExprPurity_tadPval_onlySignifTADs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  all_corr_gene_purity ~ tad_adjCombPval_log10,
  data=all_dt[all_dt$tad_adjCombPval <= 0.01,],
  xlab=myxlab,
  ylab=myylab,
  pch=16,
  cex=0.7,
  cex.axis=plotCex,
  cex.lab=plotCex
)
mtext(side=3, text = paste0(purity_ds))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))    




##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
