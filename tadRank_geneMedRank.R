
startTime <- Sys.time()

# Rscript tadRank_geneMedRank.R
# Rscript tadRank_geneMedRank.R LG1_40kb TCGAluad_norm_luad

outFolder <- file.path("TADRANK_GENEMEDRANK")
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2 | length(args) == 0 )
hicds <- args[1]
exprds <- args[2]

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- 400
myWidth <- 400

plotCex <- 1.4

inDT <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

nDS <- length(unique(paste0(inDT$hicds, inDT$exprds)))

signifThresh <- 0.01

if(length(args) == 2) {
  inDT <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]
  outprefix <- paste0(hicds, "_", exprds)
  subTit <- paste0(hicds, " - ", exprds)
} else {
  outprefix <- "all_ds"
  subTit <- paste0("all DS (n=", nDS, ")")
}

medRank_dt <- aggregate(gene_rank ~ hicds + exprds + region, data = inDT, FUN=median)
colnames(medRank_dt)[colnames(medRank_dt) == "gene_rank"] <- "med_gene_rank"

tadRank_dt <- inDT[,c("hicds", "exprds", "region", "tad_rank", "tad_adjCombPval")]
tadRank_dt <- unique(tadRank_dt)

tadRank_geneMedRank_dt <- merge(tadRank_dt, medRank_dt, all.x=TRUE, all.y=TRUE, by=c("hicds", "exprds", "region"))
stopifnot(!is.na(tadRank_geneMedRank_dt))

outFile <- file.path(outFolder, paste0(outprefix, "_medGeneRank_tadRank_allTADs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = tadRank_geneMedRank_dt$tad_rank,
  y = tadRank_geneMedRank_dt$med_gene_rank,
  xlab = paste0("TAD rank"),
  ylab = paste0("med. TAD gene ranks"),
  cex.lab = plotCex,
  cex.axis = plotCex,
  main=paste0("med. TAD gene ranks vs. TAD rank")
)
legend("topleft", legend = paste0("n = ", nrow(tadRank_geneMedRank_dt)), bty="n")
mtext(side = 3, text=paste0(subTit, " - all TADS"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

tadRank_geneMedRank_dt_signif <- tadRank_geneMedRank_dt[tadRank_geneMedRank_dt$tad_adjCombPval <= signifThresh ,]

outFile <- file.path(outFolder, paste0(outprefix, "_medGeneRank_tadRank_signifTADs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = tadRank_geneMedRank_dt_signif$tad_rank,
  y = tadRank_geneMedRank_dt_signif$med_gene_rank,
  xlab = paste0("TAD rank"),
  ylab = paste0("med. TAD gene ranks"),
  cex.lab = plotCex,
  cex.axis = plotCex,
  main=paste0("med. TAD gene ranks vs. TAD rank")
)
legend("topleft", legend = paste0("n = ", nrow(tadRank_geneMedRank_dt_signif)), bty="n")
mtext(side = 3, text=paste0(subTit, " - signif. TADS"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




cat(paste0("*** DONE \n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))