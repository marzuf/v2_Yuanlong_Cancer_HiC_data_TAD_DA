
# Rscript check_corrMeth.R

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "check_corrMeth.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script11a_file <- "11sameNbr_runEmpPvalCombined/emp_pval_combined.Rdata"
script11b_file <- "11sameNbrSpearman_runEmpPvalCombined/emp_pval_combined.Rdata"

require(foreach)
require(doMC)
registerDoMC(40)

plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.4

pipFolder<- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CHECK_CORRMETH"
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds)) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

hicds = all_hicds[1]
# all_hicds = all_hicds[1]
# all_hicds <-  all_hicds[grepl("460", all_hicds)]
all_pval_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]

  hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    a_pval <- get(load(file.path(pipOutFolder, hicds, exprds, script11a_file)))
    b_pval <- get(load(file.path(pipOutFolder, hicds, exprds, script11b_file)))
    # 
    a_adjPval <- p.adjust(a_pval, method="BH")
    b_adjPval <- p.adjust(b_pval, method="BH")
    stopifnot(setequal(names(a_adjPval), names(b_adjPval)))
    stopifnot(names(a_adjPval) == names(b_pval))
    
    data.frame(
      hicds=hicds,
      exprds=exprds,
      region=names(a_adjPval),
      PCC_adjPval=a_adjPval,
      SCC_adjPval=b_adjPval,
      stringsAsFactors = FALSE
    )
  }
  hicds_dt
}
outFile <- file.path(outFolder, "all_pval_dt.Rdata")
save(all_pval_dt, file=outFile, version=2)
cat(paste0("... written:", outFile, "\n"))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


my_x <- -log10(all_pval_dt$PCC_adjPval)
my_y <- -log10( all_pval_dt$SCC_adjPval)

x_lab <- "adj. comb. p-val - Pearson's corr. [-log10]"
y_lab <- "adj. comb. p-val - Spearman's corr. [-log10]"

nDS <- length(unique(file.path(all_pval_dt$hicds, all_pval_dt$exprds)))

outFile <- file.path(outFolder, paste0("SCC_vs_PCC_adjCombPval_log10_densplot.", plotType))

do.call(plotType, list(outFile,height=myHeight, width=myWidth))
densplot(
  x= my_x,
  y=my_y,
  xlab=x_lab,
  ylab=y_lab,
  cex=0.7,
  main=paste0("Spearman's vs. Pearson's"),
  cex.axis=plotCex,
  cex.main=plotCex,
  cex.lab=plotCex
)
mtext(side=3, text = paste0("all DS (n=", nDS, "); # TADs=", nrow(all_pval_dt)))
addCorr(x=my_x, y=my_y, bty="n", legPos="topleft")
foo <- dev.off()
cat(paste0("... written:", outFile, "\n"))

my_x <- all_pval_dt$PCC_adjPval
my_y <-  all_pval_dt$SCC_adjPval

x_lab <- "adj. comb. p-val - Pearson's corr."
y_lab <- "adj. comb. p-val - Spearman's corr."

outFile <- file.path(outFolder, paste0("SCC_vs_PCC_adjCombPval_densplot.", plotType))

do.call(plotType, list(outFile,height=myHeight, width=myWidth))
densplot(
  x= my_x,
  y=my_y,
  xlab=x_lab,
  ylab=y_lab,
  cex=0.7,
  main=paste0("Spearman's vs. Pearson's"),
  cex.axis=plotCex,
  cex.main=plotCex,
  cex.lab=plotCex
)
mtext(side=3, text = paste0("all DS (n=", nDS, "); # TADs=", nrow(all_pval_dt)))
addCorr(x=my_x, y=my_y, bty="n", legPos="topleft")
foo <- dev.off()
cat(paste0("... written:", outFile, "\n"))




