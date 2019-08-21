options(scipen=100)

setDir=""

# Rscript intersect_FDR_combPval_lolliPlot.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript intersect_FDR_combPval_lolliPlot.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich
# Rscript intersect_FDR_combPval_lolliPlot.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"
hicds="ENCSR079VIJ_G401_40kb"
exprds="TCGAkich_norm_kich"


script_name <- "intersect_FDR_combPval_lolliPlot.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

buildData <- TRUE
separateHeatmap <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

require(reshape2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("plot_lolliTAD_funct.R")

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script19_name <- "19onlyFC_SAM_emp_measurement"
script19sameNbr_name <- "19sameNbr_SAM_emp_measurement"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "INTERSECT_FDR_COMBPVAL_LOLLIPLOT"
dir.create(outFolder, recursive=TRUE)

FDRthresh_seq <- seq(from=0.1, to=0.5, by=0.1)
pvalThresh_seq <- seq(from=0.01, to=0.05, by = 0.01)

myHeightGG <- length(pvalThresh_seq)*1.2
myWidthGG <- length(FDRthresh_seq)*1.2

twoSidedStouffer <- FALSE

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2 | length(args) == 0)
hicds <- args[1]
exprds <- args[2]

if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}

if(buildData) {
  allDS_intersect_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_intersect_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat("... start ", hicds, " - " , exprds, "\n")
      # PREPARE logFC and meanCorr observed data
      fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
      corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(fc_file))
      stopifnot(file.exists(corr_file))  
      
      tad_logFC <- eval(parse(text = load(fc_file)))
      tad_meanCorr <- eval(parse(text = load(corr_file)))
      
      all_regs <- names(tad_logFC)
      stopifnot(setequal(all_regs, names(tad_meanCorr)))
      tad_logFC <- tad_logFC[all_regs]
      tad_meanCorr <- tad_meanCorr[all_regs]

      # RETRIEVE COMBINED EMP PVAL      
      comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      stopifnot(setequal(all_regs, names(comb_empPval)))
      comb_empPval <- comb_empPval[all_regs]
      # ADJUST THE PVAL
      adj_empPval_comb <- p.adjust(comb_empPval, method="BH")
      stopifnot(names(adj_empPval_comb) == all_regs)
      
      # => SIGNIF TADs FOR VARIOUS PVAL THRESHOLD
      # for each of the threshold of empPvals -> retrieve signif TADs
      pval_thresh=pvalThresh_seq[1]
      adjCombPval_signifTADs <- lapply(pvalThresh_seq, function(pval_thresh) {
        names(adj_empPval_comb)[adj_empPval_comb <= pval_thresh]
      })
      names(adjCombPval_signifTADs) <- paste0(pvalThresh_seq)
      
      # RETRIEVE FDR DATA FOR LOGFC
      # for each of the FDR threshold -> get FC cut-off and meanCorr cut-off => signif TADs those with abs(logFC) >= cut-off & meanCorr >= cut-off
      logFC_FDR_file <- file.path(pipOutFolder, hicds, exprds, script19_name, "empFDR_list.Rdata")
      stopifnot(file.exists(logFC_FDR_file))
      all_FDR <- eval(parse(text = load(logFC_FDR_file)))
      logFC_FDR <- all_FDR[["empFDR_logFC"]]  # the names is the FC threshold, the value is the FDR
      stopifnot(length(logFC_FDR) > 0)
      
      # RETRIEVE FDR DATA FOR MEAN CORR
      # the same for meanCorr
      meanCorr_FDR_file <-  file.path(pipOutFolder, hicds, exprds, script19sameNbr_name, "meanCorr_empFDR.Rdata")
      stopifnot(file.exists(meanCorr_FDR_file))
      all_corr_FDR <- eval(parse(text = load(meanCorr_FDR_file)))
      meanCorr_FDR <- all_corr_FDR[["empFDR"]]  # the names is the meanCorr threshold, the value is the FDR
      stopifnot(length(meanCorr_FDR) > 0)
      
      # => SIGNIF TADs FOR VARIOUS FDR THRESHOLD
      # for each of the FDR threshold => FC cut-off, meanCorr cut-off => signif TADs
      cutoff_fdr <- FDRthresh_seq[1]
      FDR_signifTADs <- lapply(FDRthresh_seq , function(cutoff_fdr) {
        logFC_cut_off <- min(as.numeric(as.character(na.omit(names(logFC_FDR)[logFC_FDR <= cutoff_fdr]))))  # the smallest FC cut-off that leads to desired FDR; if not returns Inf
        meanCorr_cut_off <- min(as.numeric(as.character(na.omit(names(meanCorr_FDR)[meanCorr_FDR <= cutoff_fdr]))))
        stopifnot(names(tad_logFC) == names(tad_meanCorr))
        names(tad_logFC)[ abs(tad_logFC) >= logFC_cut_off & tad_meanCorr >= meanCorr_cut_off]
      })
      names(FDR_signifTADs) <- paste0(FDRthresh_seq)
      
      fdr=FDRthresh_seq[1]
      # => ITERATE OVER VARIOUS FDR AND PVAL THRESH TO GET THE INTERSECT OF SIGNIF TADS
      all_intersectTADs <- foreach(fdr = FDRthresh_seq, .combine='rbind') %do% {
        stopifnot(paste0(fdr) %in% names(FDR_signifTADs))
        pval=pvalThresh_seq[1]
        fdr_tads <- FDR_signifTADs[[paste0(fdr)]]
        stopifnot(fdr_tads %in% all_regs)
        intersectTADs_pval_dt <- foreach(pval = pvalThresh_seq, .combine='rbind') %do% {
          stopifnot(paste0(pval) %in% names(adjCombPval_signifTADs))
          pval_tads <- adjCombPval_signifTADs[[paste0(pval)]]
          stopifnot(pval_tads %in% all_regs)
          intersect_tads_not_sorted <- intersect(pval_tads, fdr_tads)
          
          # => REORDER TO PLOT SORTING BY AVG RANK
          tad_logFC_rank <- rank(-abs(tad_logFC), ties.method = "min")
          stopifnot(abs(tad_logFC[names(tad_logFC_rank[which(tad_logFC_rank ==1)])]) == max(abs(tad_logFC)))
          tad_meanCorr_rank <- rank(-tad_meanCorr, ties.method = "min")
          stopifnot(tad_meanCorr[names(tad_meanCorr_rank[which(tad_meanCorr_rank ==1)])] == max(tad_meanCorr))
          tad_avgRank <- (tad_meanCorr_rank+tad_logFC_rank)/2
          tad_avgRank_sorted <- sort(tad_avgRank, decreasing=FALSE) # smaller rank -> better
          stopifnot(setequal(names(tad_avgRank_sorted), names(tad_logFC)))
          stopifnot(setequal(names(tad_avgRank_sorted), names(tad_meanCorr)))
          
          nPlotted <- length(intersect_tads_not_sorted)
          stopifnot(intersect_tads_not_sorted %in% names(tad_avgRank_sorted))
          intersect_tads <- names(tad_avgRank_sorted)[names(tad_avgRank_sorted) %in% intersect_tads_not_sorted]
          
          if(nPlotted > 0 & nPlotted <= 10) {
            
            stopifnot(diff(tad_avgRank_sorted[intersect_tads]) >= 0)
            stopifnot(length(intersect_tads) == nPlotted)
            
            
            
            plotList <- list()
            for(i_tad in 1:length(intersect_tads)) {
              plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                                    hicds = hicds,
                                                    all_TADs = intersect_tads[i_tad],
                                                    orderByLolli = "startPos")
            } # end-for iterating over TADs to plot
            outFile <- file.path(outFolder, paste0(hicds, exprds, "_intersect_FDR", fdr, "_adjPvalComb", pval,"_signifTADs", ".", plotType ))
            mytit <- paste0(hicds, " - ", exprds, "\n(FDR<=", fdr, "; pval<=", pval, "; sorted avgRank)")
            all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
            outHeightGG <- min(c(7 * nPlotted/2, 49))
            outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
            outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
            
            ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
            cat("... written: ", outFile, "\n")
            
        } # end if enough and not too much plots
          FDR_nSignifTADs <- length(fdr_tads)
          adjCombPval_nSignifTADs <- length(fdr_tads)
          intersect_nSignifTADs <- length(intersect_tads)
          data.frame(
            hicds = hicds,
            exprds = exprds,
            FDR_threshold = fdr,
            adjPval_threshold = pval,
            FDR_signifTADs = paste0(fdr_tads, collapse=","),
            adjCombPval_signifTADs = paste0(pval_tads, collapse=","),
            intersect_signifTADs = paste0(intersect_tads, collapse=","),
            FDR_nSignifTADs = length(fdr_tads),
            adjCombPval_nSignifTADs = length(pval_tads),
            intersect_nSignifTADs = length(intersect_tads),
            stringsAsFactors = FALSE
          )
        } # end-foreach iterating over pval cutoffs
        intersectTADs_pval_dt
      } # end-foreach iterating over FDR cutoffs
      all_intersectTADs
    } # end-foreach iterating over exprds
    exprds_intersect_DT
  }# end-foreach iterating over hicds => allDS_intersect_DT
  outFile <- file.path(outFolder, "allDS_intersect_DT.Rdata")  
  save(allDS_intersect_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else { # end-if build data
  outFile <- file.path(outFolder, "allDS_intersect_DT.Rdata")  
  allDS_intersect_DT <- eval(parse(text = load(outFile)))
}
allDS_intersect_DT$hicds <- as.character(allDS_intersect_DT$hicds)
allDS_intersect_DT$exprds <- as.character(allDS_intersect_DT$exprds)
allDS_intersect_DT$dataset <- paste0(allDS_intersect_DT$hicds, " - ", allDS_intersect_DT$exprds)



##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

