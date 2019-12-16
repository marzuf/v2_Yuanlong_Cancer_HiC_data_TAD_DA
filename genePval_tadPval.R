
# Rscript genePval_tadPval.R

inDT <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

outFolder <- "GENEPVAL_TADPVAL"
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- 400
myWidth <- myHeight

geneSignifThresh <- 0.05
tadSignifThresh <- 0.01

geneTopRank <- 100
tadTopRank <- 50

outFile <- file.path(outFolder, paste0("gene_pval_vs_tad_pval_grey_red.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  main = paste0("Gene p-val vs. TAD p-val"),
  xlab = c("TAD adj. combined p-val [-log10]"),
  ylab = c("Gene adj. p-val [-log10]"),
  x = -log10(inDT$tad_adjCombPval),
  y = -log10(inDT$adj.P.Val),
  pch = 16,
  col="grey",
  cex=0.7
)
points(
  x = -log10(inDT$tad_adjCombPval[inDT$tad_adjCombPval <= tadSignifThresh & inDT$adj.P.Val > geneSignifThresh]),
  y =  -log10(inDT$adj.P.Val[inDT$tad_adjCombPval <= tadSignifThresh & inDT$adj.P.Val > geneSignifThresh]),
  col="red",
  pch=16,
  cex=0.7
)
legend(
  "topright",
  legend = c(paste0("n = ",nrow(inDT)) ,
             paste0("gene p-val > ",geneSignifThresh, " &\nTAD p-val <= ", tadSignifThresh, "\nn=", sum(inDT$tad_adjCombPval <= tadSignifThresh & inDT$adj.P.Val > geneSignifThresh))),
  col = c("grey", "red"),
  pch = 16,
  cex = 0.7,
  bty="n"
  )

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("gene_rank_vs_tad_rank_grey_red.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  main = paste0("Gene rank vs. TAD rank"),
  xlab = c("TAD rank"),
  ylab = c("Gene rank"),
  x = inDT$tad_rank,
  y = inDT$gene_rank,
  pch = 16,
  col="grey",
  cex=0.7
)
points(
  x = inDT$tad_rank[inDT$tadTopRank <= tadSignifThresh & inDT$gene_rank > geneTopRank],
  y =  inDT$gene_rank[inDT$tadTopRank <= tadSignifThresh & inDT$gene_rank > geneTopRank],
  col="red",
  pch=16,
  cex=0.7
)
legend(
  "topright",
  legend = c(paste0("n = ",nrow(inDT)) ,
             paste0("gene rank > ",geneTopRank, " &\nTAD rank <= ", tadTopRank, "\nn=", 
                    sum(inDT$tad_rank <= tadTopRank & inDT$gene_rank > geneTopRank))),
  col = c("grey", "red"),
  pch = 16,
  cex = 0.7,
  bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

