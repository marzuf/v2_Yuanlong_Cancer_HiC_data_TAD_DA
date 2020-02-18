
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
# Rscript cmp_meanCorr_meanFC_permut.R

plotType <- "png"
myHeight <- myWidth <- 400
myHeightGG <- 7
myWidthGG <- 9
plotCex <- 1.4

outFolder <- "CMP_MEANCORR_MEANFC_PERMUT"
dir.create(outFolder, recursive = TRUE)
  
all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
hicds = all_hicds[1]
all_hicds <- all_hicds[grep("NCI-H460", all_hicds)]

myHicds <- "ENCSR489OCU_NCI-H460_"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")



tadSignifThresh <- 0.01

plot_dt  <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  mean_FC <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "3_runMeanTADLogFC", "all_meanLogFC_TAD.Rdata")))
  mean_corr <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata")))
  
  
  empPval_FC <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "9_runEmpPvalMeanTADLogFC", "emp_pval_meanLogFC.Rdata")))
  empPval_FC <- p.adjust(empPval_FC, method="BH")
  
  empPval_Corr <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "10sameNbr_runEmpPvalMeanTADCorr", "emp_pval_meanCorr.Rdata")))
  empPval_Corr <- p.adjust(empPval_Corr, method="BH")
  
  empPval <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
  empPval <- p.adjust(empPval, method="BH")
  
  
  data.frame(
    hicds= hicds,
    meanFC = as.numeric(mean_FC[names(mean_corr)]),
    meanCorr = as.numeric(mean_corr),
    adjCombPval = empPval[names(mean_corr)],
    adjFCPval = empPval_FC[names(mean_corr)],
    adjCorrPval = empPval_Corr[names(mean_corr)],
    region = names(mean_corr),
    stringsAsFactors = FALSE
  )
}
stopifnot(!is.na(plot_dt))
plot_dt$signif <- ifelse(plot_dt$adjCombPval <= tadSignifThresh, "signif.", "not signif.")
plot_dt$hicds_lab <- gsub(myHicds, "", plot_dt$hicds)
plot_dt$signif_col <- ifelse(plot_dt$signif == "signif.", "red", "darkgrey")

for(rd_type in unique(plot_dt$hicds_lab)) {
  
  my_x <- plot_dt$meanFC[plot_dt$hicds_lab == rd_type]
  my_y <- plot_dt$meanCorr[plot_dt$hicds_lab == rd_type]
  
  outFile <- file.path(outFolder, paste0("meanCorr_meanFC_", rd_type, "_dotplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  plot(my_x, 
       my_y,
       xlim = range(plot_dt$meanFC),
       ylim = range(plot_dt$meanCorr),
       main = paste0(rd_type),
       xlab = "TAD meanFC",
       ylab = "TAD meanCorr",
       cex.main = plotCex,
       cex.axis = plotCex,
       cex.lab = plotCex,
       pch = 16,
       cex = 0.7,
       cols = plot_dt$signif_col
       )
  points(
    x= plot_dt$meanFC[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ] ,
    y = plot_dt$meanCorr[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ],
    col="red",
    pch=16
  )
  mtext(side=3, paste0(myHicds, " - signif. thresh <= ",tadSignifThresh ))
  # addCorr(my_x,my_y,bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  my_x <- -log10(plot_dt$adjCombPval[plot_dt$hicds_lab == rd_type])
  my_y <- plot_dt$meanCorr[plot_dt$hicds_lab == rd_type]
  
  outFile <- file.path(outFolder, paste0("meanCorr_adjPvalComb_", rd_type, "_dotplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  plot(my_x, 
       my_y,
       xlim = range( -log10(plot_dt$adjCombPval)),
       ylim = range(plot_dt$meanCorr),
       main = paste0(rd_type),
       xlab = "TAD adj. comb. pval [-log10]",
       ylab = "TAD meanCorr",
       cex.main = plotCex,
       cex.axis = plotCex,
       cex.lab = plotCex,
       pch = 16,
       cex = 0.7,
       cols = plot_dt$signif_col
  )
  points(
    x= -log10(plot_dt$adjCombPval[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ] ),
    y = plot_dt$meanCorr[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ],
    col="red",
    pch=16
  )
  mtext(side=3, paste0(myHicds, " - signif. thresh <= ",tadSignifThresh ))
  # addCorr(my_x,my_y,bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  my_x <- -log10(plot_dt$adjCombPval[plot_dt$hicds_lab == rd_type])
  my_y <- plot_dt$meanFC[plot_dt$hicds_lab == rd_type]
  
  outFile <- file.path(outFolder, paste0("meanFC_adjPvalComb_", rd_type, "_dotplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  plot(my_x, 
       my_y,
       xlim = range( -log10(plot_dt$adjCombPval)),
       ylim = range(plot_dt$meanFC),
       main = paste0(rd_type),
       xlab = "TAD adj. comb. pval [-log10]",
       ylab = "TAD meanFC",
       cex.main = plotCex,
       cex.axis = plotCex,
       cex.lab = plotCex,
       pch = 16,
       cex = 0.7,
       cols = plot_dt$signif_col
  )
  points(
    x=  -log10(plot_dt$adjCombPval[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ]) ,
    y = plot_dt$meanFC[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ],
    col="red",
    pch=16
  )
  mtext(side=3, paste0(myHicds, " - signif. thresh <= ",tadSignifThresh ))
  # addCorr(my_x,my_y,bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  my_x <- -log10(plot_dt$adjCombPval[plot_dt$hicds_lab == rd_type])
  my_y <- -log10(plot_dt$adjFCPval[plot_dt$hicds_lab == rd_type])
  
  outFile <- file.path(outFolder, paste0("adjPvalFC_adjPvalComb_", rd_type, "_dotplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  plot(my_x, 
       my_y,
       xlim = range( -log10(plot_dt$adjCombPval)),
       ylim = range(-log10(plot_dt$adjFCPval)),
       main = paste0(rd_type),
       xlab = "TAD adj. comb. pval [-log10]",
       ylab = "TAD adj. FC pval [-log10]",
       cex.main = plotCex,
       cex.axis = plotCex,
       cex.lab = plotCex,
       pch = 16,
       cex = 0.7,
       cols = plot_dt$signif_col
  )
  points(
    x=  -log10(plot_dt$adjCombPval[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ]) ,
    y = -log10(plot_dt$adjFCPval[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ]),
    col="red",
    pch=16
  )
  mtext(side=3, paste0(myHicds, " - signif. thresh <= ",tadSignifThresh ))
  # addCorr(my_x,my_y,bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  my_x <- -log10(plot_dt$adjCombPval[plot_dt$hicds_lab == rd_type])
  my_y <- -log10(plot_dt$adjCorrPval[plot_dt$hicds_lab == rd_type])
  
  outFile <- file.path(outFolder, paste0("adjPvalCorr_adjPvalComb_", rd_type, "_dotplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  plot(my_x, 
       my_y,
       xlim = range( -log10(plot_dt$adjCombPval)),
       ylim = range(-log10(plot_dt$adjCorrPval)),
       main = paste0(rd_type),
       xlab = "TAD adj. comb. pval [-log10]",
       ylab = "TAD adj. Corr pval [-log10]",
       cex.main = plotCex,
       cex.axis = plotCex,
       cex.lab = plotCex,
       pch = 16,
       cex = 0.7,
       cols = plot_dt$signif_col
  )
  points(
    x=  -log10(plot_dt$adjCombPval[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ]) ,
    y = -log10(plot_dt$adjCorrPval[plot_dt$hicds_lab == rd_type & plot_dt$signif=="signif." ]),
    col="red",
    pch=16
  )
  mtext(side=3, paste0(myHicds, " - signif. thresh <= ",tadSignifThresh ))
  # addCorr(my_x,my_y,bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
}











