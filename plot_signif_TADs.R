options(scipen=100)

setDir=""

# Rscript plot_signif_TADs.R <FDR_threshold> <p_value_threshold> <hicds> <exprds>
# Rscript plot_signif_TADs.R <FDR_threshold> <p_value_threshold> # to run all datasets in one shot
# Rscript plot_signif_TADs.R 0.1 0.01 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_signif_TADs.R 0.2 0.01 

script_name <- "plot_signif_TADs.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

buildData <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(ggplot2)
require(reshape2)


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("plot_lolliTAD_funct.R")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.4

nToPlot <- 10

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script8c_name <- "8cOnlyRatioDownFastSave_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script19_name <- "19onlyFC_SAM_emp_measurement"
script19sameNbr_name <- "19sameNbr_SAM_emp_measurement"
script11same_name <- "11sameNbr_runEmpPvalCombined"

pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))


outFolder <- "PLOT_SIGNIF_TADS"
dir.create(outFolder, recursive=TRUE)

twoSidedStouffer <- FALSE

FDRthresh=0.1
pvalThresh=0.01
hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2 | length(args) == 4)
FDRthresh <- as.numeric(args[1])
pvalThresh <- as.numeric(args[2])
stopifnot(!is.na(FDRthresh))
stopifnot(!is.na(pvalThresh))
stopifnot(FDRthresh >= 0 & FDRthresh <=1 )
stopifnot(pvalThresh >= 0 & pvalThresh <=1 )
hicds <- args[3]
exprds <- args[4]

if(length(args) == 2) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}

### BUILD SIGNIF ALONG FDR THRESH
cat("... start retrieving FDR signif. TADs\n")

allDS_signifFDR_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_signifFDR_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    ##########> FDR SIGNIF TADS
    # RETRIEVE FDR DATA FOR MEAN LOGFC
    cat("... ", hicds , " - ", exprds, ": retrieved signif. TADs, FDR thresh. = ", FDRthresh, "\n")
    logFC_FDR_file <- file.path(pipOutFolder, hicds, exprds, script19_name, "empFDR_list.Rdata")
    stopifnot(file.exists(logFC_FDR_file))
    all_FDR <- eval(parse(text = load(logFC_FDR_file)))
    logFC_FDR <- all_FDR[["empFDR_logFC"]]  # the names is the FC threshold, the value is the FDR
    stopifnot(length(logFC_FDR) > 0)
    # RETRIEVE FDR DATA FOR MEAN CORR




    # RETRIEVE FDR DATA FOR MEAN CORR
    # the same for meanCorr
    meanCorr_FDR_file <-  file.path(pipOutFolder, hicds, exprds, script19sameNbr_name, "meanCorr_empFDR.Rdata")
    stopifnot(file.exists(meanCorr_FDR_file))
    all_corr_FDR <- eval(parse(text = load(meanCorr_FDR_file)))
    meanCorr_FDR <- all_corr_FDR[["empFDR"]]  # the names is the meanCorr threshold, the value is the FDR
    stopifnot(length(meanCorr_FDR) > 0)
      



    # PREPARE logFC and meanCorr observed data
    fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
    corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
    stopifnot(file.exists(fc_file))
    stopifnot(file.exists(corr_file))  
    tad_logFC <- eval(parse(text = load(fc_file)))
    tad_meanCorr <- eval(parse(text = load(corr_file)))
    all_regs <- names(tad_logFC)
    stopifnot(setequal(all_regs, names(tad_meanCorr)))
    # => REORDER TO PLOT THE "TOP RANKING"
    tad_logFC_rank <- rank(-abs(tad_logFC), ties.method = "min")
    stopifnot(abs(tad_logFC[names(tad_logFC_rank[which(tad_logFC_rank ==1)])]) == max(abs(tad_logFC)))
    tad_meanCorr_rank <- rank(-tad_meanCorr, ties.method = "min")
    stopifnot(tad_meanCorr[names(tad_meanCorr_rank[which(tad_meanCorr_rank ==1)])] == max(tad_meanCorr))
    tad_avgRank <- (tad_meanCorr_rank+tad_logFC_rank)/2
    tad_avgRank_sorted <- sort(tad_avgRank, decreasing=FALSE) # smaller rank -> better
    stopifnot(setequal(names(tad_avgRank_sorted), names(tad_logFC)))
    stopifnot(setequal(names(tad_avgRank_sorted), names(tad_meanCorr)))
    # reorder to have top ranking first
    tad_logFC <- tad_logFC[names(tad_avgRank_sorted)]
    tad_meanCorr <- tad_meanCorr[names(tad_avgRank_sorted)]
    
    plotDT <- data.frame(regions = names(tad_logFC), meanFC = tad_logFC, meanCorr = tad_meanCorr)
    
    
    # => SIGNIF TADs FOR THE DESIRED FDR THRESHOLD
    logFC_cut_off <- min(as.numeric(as.character(na.omit(names(logFC_FDR)[logFC_FDR <= FDRthresh]))))  # the smallest FC cut-off that leads to desired FDR; if not returns Inf
    meanCorr_cut_off <- min(as.numeric(as.character(na.omit(names(meanCorr_FDR)[meanCorr_FDR <= FDRthresh]))))
    
    FDR_signifTADs <- names(tad_logFC)[( abs(tad_logFC) >= logFC_cut_off & tad_meanCorr >= meanCorr_cut_off)]
    
    if(length(FDR_signifTADs) > 0) {
      nPlotted <- min(c(nToPlot, length(FDR_signifTADs)))
      plot_FDR_signifTADs <-   FDR_signifTADs[1:nPlotted]
      cat("... ", hicds , " - ", exprds, ": signif. TADs, FDR thresh - start plotting (max. ", nToPlot, "; found: ", length(FDR_signifTADs), ")\n")  
      plotList <- list()
      for(i_tad in 1:length(plot_FDR_signifTADs)) {
        plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                          hicds = hicds,
                                          all_TADs = plot_FDR_signifTADs[i_tad],
                                          orderByLolli = "startPos")
      } # end-for iterating over TADs to plot
      outFile <- file.path(outFolder, paste0(hicds, exprds, "_FDRsignifTADs",  "_nToPlot", nToPlot, ".", plotType ))
      mytit <- paste0(hicds, " - ", exprds, " - top ", nToPlot, "\n(FDR <=", FDRthresh, "; sorted avgRank)")
      all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nToPlot == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      outHeightGG <- min(c(7 * nPlotted/2, 49))
      outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
      outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
      
      ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
      cat("... written: ", outFile, "\n")
      
      
      stopifnot(names(tad_logFC) == names(tad_meanCorr))
      
      outFileSuffix <- file.path(outFolder, paste0(hicds, "_", exprds, "_meanCorr_vs_meanFC_FDRthresh", FDRthresh))
      
      plot_meanFC_meanCorr_FDRthresh(
        dataDT=plotDT,
        var1 = "meanFC",
        var2 = "meanCorr",
        annotCol = "regions",
        var1_cutoff = logFC_cut_off,
        var2_cutoff = meanCorr_cut_off,
        abs_var1 = TRUE,
        plotTit = paste0(hicds, " - ", exprds),
        plotSub = paste0("# TADs=", nrow(plotDT), " (# signif.=", length(FDR_signifTADs), ")"),
        fileSuffix = outFileSuffix,
        plotType= plotType,
        legSuppTxt = paste0("FDR thresh.: ",FDRthresh),
        mTextLine = 0,
        toPlot = c("grey","density"),
        fullName = FALSE,
        twoInOne = TRUE
      )
    }else{
      # NO SIGNIF TADs FOUND
      cat("... ", hicds , " - ", exprds, ": no signif. FDR TADs found \n")  
    }
    
     #stop("--ok\n")
    
    ##########> ADJ. COMB PVAL SIGNIF TADS



    # RETRIEVE COMBINED EMP PVAL      
    comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
    stopifnot(file.exists(comb_empPval_file))
    comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
    stopifnot(setequal(all_regs, names(comb_empPval)))
    comb_empPval <- comb_empPval[all_regs]
    # ADJUST THE PVAL
    adj_empPval_comb <- p.adjust(comb_empPval, method="BH")
    stopifnot(names(adj_empPval_comb) == all_regs)

    adj_empPval_comb <- sort(adj_empPval_comb, decreasing=FALSE) # => in order to plot the top-ranking; smaller pval = better
    # => SIGNIF TADs FOR THE DESIRED PVAL THRESHOLD
    adjCombPval_signifTADs <- names(adj_empPval_comb)[(adj_empPval_comb <= pvalThresh)]
    if(length(adjCombPval_signifTADs) > 0) {
      nPlotted <- min(c(nToPlot, length(adjCombPval_signifTADs)))
      plot_adjCombPval_signifTADs <-   adjCombPval_signifTADs[1:nPlotted]
      cat("... ", hicds , " - ", exprds, ": signif. TADs, adj. pval. comb. thresh - start plotting (max. ", nToPlot, "; found: ", length(adjCombPval_signifTADs), ")\n")
      plotList <- list()
      for(i_tad in 1:length(plot_adjCombPval_signifTADs)) {
        plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                          hicds = hicds,
                                          all_TADs = plot_adjCombPval_signifTADs[i_tad],
                                          orderByLolli = "startPos")
      } # end-for iterating over TADs to plot
      outFile <- file.path(outFolder, paste0(hicds, exprds, "_adjCombPvalSignifTADs",  "_nToPlot", nToPlot, ".", plotType ))
      mytit <- paste0(hicds, " - ", exprds, " - top ", nToPlot, "\n(adj. comb. pval. <= ", pvalThresh, "; sorted adjPvalComb)")
      all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nToPlot == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      outHeightGG <- min(c(7 * nPlotted/2, 49))
      outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
      ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
      cat("... written: ", outFile, "\n")
      
      outFileSuffix <- file.path(outFolder, paste0(hicds, "_", exprds, "_meanCorr_vs_meanFC_pvalThresh", pvalThresh))
      

      plot_meanFC_meanCorr_signifTADs(
        dataDT=plotDT,
        var1 = "meanFC",
        var2 = "meanCorr",
        annotCol = "regions",
        signifCol = "regions",
        signifTADs = adjCombPval_signifTADs,
        abs_var1 = FALSE,
        plotTit = paste0(hicds, " - ", exprds),
        plotSub = paste0("# TADs=", nrow(plotDT), " (# signif.=", length(adjCombPval_signifTADs), ")"),
        fileSuffix = outFileSuffix,
        plotType= plotType,
        mTextLine = 0,
        toPlot = c("grey","density"),
        fullName = FALSE,
        twoInOne = TRUE,
        legTxt = paste0("Pval. thresh.:", pvalThresh)
      )
    }else{
      # NO SIGNIF TADs FOUND
      cat("... ", hicds , " - ", exprds, ": no signif. adj. pval. comb. TADs found \n")  
    }
  } # end-foreach iterating over exprds
} # end-foreach iterating over hicds





##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

