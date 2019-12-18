
startTime <- Sys.time()

# Rscript meanCorr_with_nTFs.R

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

plotType <- "png"
myHeight <- 400
myWidth <- 400
plotCex <- 1.4



all_db <- c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg")
dsIn="crisp"
foo <- foreach(dsIn = all_db) %dopar% {
  
  outFolder <- paste0("MEANCORR_WITH_NTFS")
  dir.create(outFolder, recursive = TRUE)
  
  tf_dt <- get(load(file.path(paste0("TFS_BY_TADS_", toupper(dsIn)), "all_dt.Rdata")))
  colnames(tf_dt)[colnames(tf_dt) == "targetRegion"] <- "region"
  
  tf_corr_dt <- merge(final_dt, tf_dt, by=c("hicds", "exprds", "region"), all.x=TRUE, all.y=TRUE)
  stopifnot(na.omit(tf_corr_dt$regSymbol) > 0)
  tf_corr_dt$regSymbol[is.na(tf_corr_dt$regSymbol)] <- 0
  
  plotTit <- paste0("all DS (", length(unique(file.path(tf_corr_dt$hicds, tf_corr_dt$exprds))),")")
  subTit <- paste0(dsIn)

  outFile <- file.path(outFolder, paste0(dsIn, "_meanCorr_nbrTF_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    y = tf_corr_dt$meanCorr,
    x = tf_corr_dt$regSymbol,
    cex = 0.7,
    ylab = "TAD meanCorr",
    xlab = "TAD # reg. elements",
    main = plotTit,
    cex.axis=plotCex,
    cex.lab = plotCex
  )
  mtext(side=3, text = subTit)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  tf_corr_dt$adjPvalComb.y <- NULL
  # stopifnot(round(tf_corr_dt$adjPvalComb.x,3) == round(tf_corr_dt$adjPvalComb.y,3)) // there are some NA for those that did not have reg symbols
  
  
  outFile <- file.path(outFolder, paste0(dsIn, "_adjPvalComb_nbrTF_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    y = -log10(tf_corr_dt$adjPvalComb.x),
    x = tf_corr_dt$regSymbol,
    cex = 0.7,
    ylab = "TAD adj. comb. p-val [-log10]",
    xlab = "TAD # reg. elements",
    main = plotTit,
    cex.axis=plotCex,
    cex.lab = plotCex
  )
  mtext(side=3, text = paste0(subTit, " (n=", nrow(tf_corr_dt), ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  signif_tf_corr_dt <- tf_corr_dt[tf_corr_dt$adjPvalComb.x <= 0.01,]
  
  outFile <- file.path(outFolder, paste0(dsIn, "_adjPvalComb_nbrTF_densplot_signif.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    y = -log10(signif_tf_corr_dt$adjPvalComb.x),
    x = signif_tf_corr_dt$regSymbol,
    cex = 0.7,
    ylab = "TAD adj. comb. p-val [-log10]",
    xlab = "TAD # reg. elements",
    main = plotTit,
    cex.axis=plotCex,
    cex.lab = plotCex
  )
  mtext(side=3, text = paste0(subTit, " - signif. TADs (n=", nrow(signif_tf_corr_dt), ")"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  tf_corr_dt$regSymbolCat <- ifelse(tf_corr_dt$regSymbol == 0, "0 TF", "> 0 TF(s)")
  
  outFile <- file.path(outFolder, paste0(dsIn, "_meanCorr_nbrTF_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(
    meanCorr ~ regSymbolCat,
    data = tf_corr_dt,
    main  = paste0(plotTit, " - ", dsIn),
    xlab = "",
    ylab = paste0("TAD meanCorr"),
    cex.axis=plotCex,
    cex.lab = plotCex
  )
  mtext(side=3, text = paste0("# 0 TF = ", sum(tf_corr_dt$regSymbolCat ==  "0 TF"), 
                              "; # > 0 TF(s) = ",  sum(tf_corr_dt$regSymbolCat ==  "> 0 TF(s)")))  
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  tf_corr_dt$regSymbolCat <- ifelse(tf_corr_dt$regSymbol == 0, "0 TF", 
                                    ifelse(tf_corr_dt$regSymbol == 1, "1 TF", "> 1 TFs"))
  

    
  tf_corr_dt$regSymbolCat <- factor(tf_corr_dt$regSymbolCat, levels=c("0 TF", "1 TF", "> 1 TFs"))
  
  outFile <- file.path(outFolder, paste0(dsIn, "_meanCorr_nbrTF_with1_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(
    meanCorr ~ regSymbolCat,
    data = tf_corr_dt,
    main  = paste0(plotTit, " - ", dsIn),
    xlab = "",
    ylab = paste0("TAD meanCorr"),
    cex.axis=plotCex,
    cex.lab = plotCex
  )
  mtext(side=3, text = paste0("# 0 TF = ", sum(tf_corr_dt$regSymbolCat ==  "0 TF"), 
                              "# 1 TF = ", sum(tf_corr_dt$regSymbolCat ==  "1 TF"), 
                              "; # > 1 TFs = ",  sum(tf_corr_dt$regSymbolCat ==  "> 1 TFs")))  
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  
  tf_corr_dt$pvalSignif <- tf_corr_dt$adjPvalComb.x <= 0.01
  
  
  tf_corr_dt$tfAndSignif <- tf_corr_dt$pvalSignif & tf_corr_dt$regSymbol>0
  
  tf_corr_dt$noTfAndSignif <- tf_corr_dt$pvalSignif & tf_corr_dt$regSymbol==0
  
  nTFandSignif <- aggregate( tfAndSignif ~ hicds+exprds, data=tf_corr_dt, FUN=sum)
  nNoTFandSignif <- aggregate( noTfAndSignif ~ hicds+exprds, data=tf_corr_dt, FUN=sum)
  
  tfAndSignif_dt <- merge(nTFandSignif, nNoTFandSignif,by=c("hicds", "exprds"), all=TRUE)
  
  
  outFile <- file.path(outFolder, paste0(dsIn, "_nbrTADs_signif_TF_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(tfAndSignif_dt[,-c(1,2)],
    # meanCorr ~ regSymbolCat,
    # data = tf_corr_dt,
    main  = paste0(plotTit, " - ", dsIn),
    xlab = "",
    ylab = paste0("# TADs"),
    cex.axis=plotCex,
    cex.lab = plotCex
  )
  mtext(side=3, text = paste0("signif: adj. comb. p-val <= 0.01")) 
  
  # mtext(side=3, text = paste0("# 0 TF = ", sum(tf_corr_dt$regSymbolCat ==  "0 TF"), 
  #                             "; # > 0 TF(s) = ",  sum(tf_corr_dt$regSymbolCat ==  "> 0 TF(s)")))  
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}



cat("***DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))