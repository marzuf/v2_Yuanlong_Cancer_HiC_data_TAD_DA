


require(foreach)
require(doMC)
registerDoMC(40)

# Rscript nbr_signif_all_random.R

outFolder <- "NBR_SIGNIF_ALL_RANDOM"
dir.create(outFolder)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <-  "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4

qt_dt<- read.delim("hicds_sparsity.csv", sep="\t", header = TRUE)

pipOutFolder <- "PIPELINE/OUTPUT_FOLDER"

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[! (grepl("PERMUTG2T", all_hicds) | grepl("RANDOMNBR", all_hicds) | grepl( "RANDOMSHIFT", all_hicds))]
hicds = all_hicds[1]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

tad_signifThresh <- 0.01

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    
    combPval_file <- file.path(pipOutFolder, hicds, exprds, "11random_runEmpPvalCombined", "emp_pval_combined.Rdata")
    stopifnot(file.exists(combPval_file))
    pval <- get(load(combPval_file))
    adjPval <- p.adjust(pval, method="BH")
    
    data.frame(
      hicds=hicds,
      exprds=exprds,
      nSignif=sum(adjPval <= tad_signifThresh),
      nTot = length(adjPval),
      stringsAsFactors = FALSE
    )
  }
}


all_dt$hicds_lab <- gsub(".+_(.+?)_40kb","\\1",  all_dt$hicds)
all_dt$hicds_lab <- ifelse(grepl("RANDOM",all_dt$hicds_lab) | grepl("PERMUT", all_dt$hicds_lab),all_dt$hicds_lab,  "OBSERVED" )

stopifnot(diff(table(all_dt$hicds_lab)) == 0)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path(outFolder, paste0("nbrSignifTADs_densityplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  split(all_dt$nSignif, all_dt$hicds_lab),
  plotTit = paste0("# signif. TADs")
)
mtext(side=3, text = paste0("(sameNbr) adj. comb. p-val <= ", tad_signifThresh), font = 3)
foo <- dev.off()

outFile <- file.path(outFolder, paste0("nbrSignifTADs_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
par(mar = par()$mar + c(9,0,0,0))
boxplot(nSignif~hicds_lab, data = all_dt, main = "# signif. TADs", xlab="", ylab="# signif. TADs", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
mtext(side=3, text = paste0("(sameNbr) adj. comb. p-val <= ", tad_signifThresh), font = 3)
foo <- dev.off()

all_dt$ratioSignif <- all_dt$nSignif/all_dt$nTot
outFile <- file.path(outFolder, paste0("ratioSignifTADs_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
par(mar = par()$mar + c(9,0,0,0))
boxplot(ratioSignif~hicds_lab, data = all_dt, main = "Ratio signif. TADs", xlab="", ylab="Ratio signif. TADs", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
mtext(side=3, text = paste0("(sameNbr) adj. comb. p-val <= ", tad_signifThresh), font = 3)
foo <- dev.off()





