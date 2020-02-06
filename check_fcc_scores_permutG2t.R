

options(scipen=100)

# Rscript check_fcc_scores_permutG2t.R

script_name <- "check_fcc_scores_permutG2t.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildTable <- FALSE

require(flux)
require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)
registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("CHECK_FCC_SCORES_PERMUTG2T")
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.4
plotType <- "png"
myHeight <- 400
myWidth <- 400
myWidthLeg <- 500


hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"
# exprds <- "TCGAluad_mutKRAS_mutEGFR"

nCol <- 1000

ratioFC_permDT <- get(load(file.path("PERMG2T_TAD_FC_RATIO", hicds, exprds, "ratioFC_1000Permut_permDT.Rdata")))
ratioFC_permDT <- ratioFC_permDT[sort(rownames(ratioFC_permDT)),]
ratioDown_permDT <- get(load(paste0(hicds, "_", exprds, "_ratioDown_1000Permut_permDT.Rdata")))
fcc_permDT <- get(load(paste0(hicds, "_", exprds, "_fcc_1000Permut_permDT.Rdata")))

stopifnot( (2*ratioDown_permDT[1,1]-1)*(2*ratioFC_permDT[1,1]-1) == fcc_permDT[1,1])
stopifnot(rownames(ratioFC_permDT) == rownames(ratioDown_permDT))
stopifnot(rownames(ratioFC_permDT) == rownames(fcc_permDT))





  
  
  fcc_fract <- seq(from=-1, to=1, by=0.25)
  # fcc_fract_names <- paste0("FCC > ", fcc_fract[1:(length(fcc_fract)-1)], " and FCC <= ",fcc_fract[2:length(fcc_fract)])
  fcc_fract_names <- paste0("FCC \u2208 ]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
  fcc_fract_names <- paste0("]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
  fcc_fract_names[fcc_fract_names == "]-1, -0.75]"] <- "[-1, -0.75]"
  
  fccScores_fract <- sapply(as.numeric(fcc_permDT[,1:nCol]), function(x) which(hist(x, breaks=fcc_fract, plot=F)$counts == 1))
  
  require(ggsci)
  ggsci_pal <- "lancet"
  ggsci_subpal <- ""
  myPals <-  eval(parse(text=paste0("pal_", ggsci_pal, "(", ggsci_subpal, ")")))(length(unique(fcc_fract_names)))
  myPals <- rev(myPals)
  
  outFile <- file.path(outFolder, paste0(hicds, exprds, "_ratioDown_ratioFC_colFCCfract_", nCol, "col.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidthLeg))
  
  par(xpd = T, mar = par()$mar + c(0,0,0,6))
  
  plot(
    x = as.numeric(ratioFC_permDT[,1:nCol]),
    y = as.numeric(ratioDown_permDT[,1:nCol]),
    xlab=paste0("ratioNegFC_", nCol, "col"),
    ylab=paste0("ratioDown_", nCol, "col"),
    col = myPals[fccScores_fract],
    pch = 16,
    cex = 0.7,
    cex.lab = plotCex,
    cex.axis = plotCex,
    main = paste0(hicds, " - ", exprds)
  )
  mtext(side=3, text = paste0("# col = ", nCol, "; # TADs = ", nrow(ratioFC_permDT)))
  
  # legend("bottomright",
  legend(1.1,1,
         title="FCC score",
         legend= rev(fcc_fract_names),
         col =rev(myPals),
         pch=20,
         cex = 0.8,
         ncol=1,
         horiz = F,
         bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))












