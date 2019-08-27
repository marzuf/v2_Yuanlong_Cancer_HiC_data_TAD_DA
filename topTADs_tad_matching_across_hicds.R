options(scipen=100)

setDir="/media/electron"
setDir=""

SSHFS=TRUE
SSHFS=FALSE

# Rscript topTADs_tad_matching_across_hicds.R 0.2 0.01
# Rscript topTADs_tad_matching_across_hicds.R 0.2 0.01 K562_40kb TCGAlaml_wt_mutFLT3

script_name <- "topTADs_tad_matching_across_hicds.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("plot_lolliTAD_funct.R")

buildTable <- FALSE

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"
hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAluad_norm_luad"


outFolder <- "TOP_TADS_TAD_MATCHING_ACROSS_HICDS"
dir.create(outFolder, recursive=TRUE)

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script11_name <- "11sameNbr_runEmpPvalCombined"
script19_name <- "19onlyFC_SAM_emp_measurement"
script19sameNbr_name <- "19sameNbr_SAM_emp_measurement"

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.2

inFolder <- "TAD_MATCHING_ACROSS_HICDS"
inFile <- file.path(inFolder, "all_matching_dt.Rdata")
stopifnot(file.exists(inFile))

all_matching_dt <- get(load(inFile))

all_matching_dt$overlapBp_ratio <- all_matching_dt$overlapBP/all_matching_dt$ref_totBp
all_matching_dt$overlapGenes_ratio <- all_matching_dt$nOverlapGenes/all_matching_dt$ref_nGenes
stopifnot(na.omit(all_matching_dt$overlapBp_ratio) >= 0 & na.omit(all_matching_dt$overlapBp_ratio) <= 1)
stopifnot(na.omit(all_matching_dt$overlapGenes_ratio) >= 0 & na.omit(all_matching_dt$overlapGenes_ratio <= 1))

matching_ratio_thresh <- 0.8
conservQt <- 0.95
pvalQt <- 0.95


FDRthresh=0.2
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

cat(paste0("> FDRthresh\t=\t", FDRthresh, "\n"))
cat(paste0("> pvalThresh\t=\t", pvalThresh, "\n"))

if(length(args) == 2) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}

# 
# hicds = all_hicds[7]
# exprds = all_exprds[[paste0(hicds)]][4]

if(buildTable) {
  
  all_data <- foreach(hicds = all_hicds) %do% {
    
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_data <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      regionfile <- file.path(pipOutFolder,  hicds, exprds, script0_name, "pipeline_regionList.Rdata")
      stopifnot(file.exists(regionfile))
      regionList <- get(load(regionfile))
      
      cat("... ", hicds , " - ", exprds, ": retrieved signif. TADs, comb. emp. pval = ", pvalThresh, "\n")
      empPval_file <- file.path(pipOutFolder, hicds, exprds, script11_name, "emp_pval_combined.Rdata")
      stopifnot(file.exists(empPval_file))
      empPval <- get(load(empPval_file))
      adj_empPval <- p.adjust(empPval, method="BH")
      adjEmpPval_signifTADs <- names(adj_empPval)[adj_empPval <= pvalThresh]
      
      
      cat("... ", hicds , " - ", exprds, ": retrieved signif. TADs, FDR thresh. = ", FDRthresh, "\n")
      # PREPARE logFC and meanCorr observed data
      fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
      corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(fc_file))
      stopifnot(file.exists(corr_file))  
      tad_logFC <- eval(parse(text = load(fc_file)))
      tad_meanCorr <- eval(parse(text = load(corr_file)))
      all_regs <- names(tad_logFC)
      stopifnot(setequal(all_regs, names(tad_meanCorr)))
      stopifnot(setequal(all_regs, names(adj_empPval)))
      stopifnot(length(all_regs) == length(regionList))
      # => REORDER TO PLOT THE "TOP RANKING"
      tad_logFC <- tad_logFC[all_regs]
      tad_meanCorr <- tad_meanCorr[all_regs]
      
      
      ## => retrieve the FC and corr cut-off yielding the desired FDR
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
      
      
      # => SIGNIF TADs FOR THE DESIRED FDR THRESHOLD
      logFC_cut_off <- min(as.numeric(as.character(na.omit(names(logFC_FDR)[logFC_FDR <= FDRthresh]))))  # the smallest FC cut-off that leads to desired FDR; if not returns Inf
      meanCorr_cut_off <- min(as.numeric(as.character(na.omit(names(meanCorr_FDR)[meanCorr_FDR <= FDRthresh]))))
      stopifnot(names(tad_logFC) == names(tad_meanCorr))
      FDR_signifTADs <- names(tad_logFC)[( abs(tad_logFC) >= logFC_cut_off & tad_meanCorr >= meanCorr_cut_off)]
      
      stopifnot(file.path(hicds, exprds) %in% all_matching_dt$ref_dataset)
      ds_matching_dt <- all_matching_dt[all_matching_dt$ref_dataset == file.path(hicds, exprds) ,]  
      stopifnot(nrow(ds_matching_dt) > 0)
      stopifnot(c(FDR_signifTADs, adjEmpPval_signifTADs) %in% ds_matching_dt$refID)
      
      stopifnot(length(unique(ds_matching_dt$refID)) == length(all_regs))
      
      ds_matching_dt$isSignifFDR <- ds_matching_dt$refID %in% FDR_signifTADs
      ds_matching_dt$isSignifEmpPval <- ds_matching_dt$refID %in% adjEmpPval_signifTADs
      
      bp_nOverThresh_FDR <- sapply(FDR_signifTADs, function(fdr_signif_tad) {
        sum(na.omit(ds_matching_dt$overlapBp_ratio[ds_matching_dt$refID == fdr_signif_tad]) >= matching_ratio_thresh)  
      })
      bp_nOverThresh_empPval <- sapply(adjEmpPval_signifTADs, function(ep_signif_tad) {
        sum(na.omit(ds_matching_dt$overlapBp_ratio[ds_matching_dt$refID == ep_signif_tad]) >= matching_ratio_thresh)  
      })
      genes_nOverThresh_FDR <- sapply(FDR_signifTADs, function(fdr_signif_tad) {
        sum(na.omit(ds_matching_dt$overlapGenes_ratio[ds_matching_dt$refID == fdr_signif_tad]) >= matching_ratio_thresh)  
      })
      genes_nOverThresh_empPval <- sapply(adjEmpPval_signifTADs, function(ep_signif_tad) {
        sum(na.omit(ds_matching_dt$overlapGenes_ratio[ds_matching_dt$refID == ep_signif_tad]) >= matching_ratio_thresh)  
      })
      
      emptyDT <- data.frame(region = character(0),
                            nOverThresh = numeric(0),
                            isSignifFDR = logical(0),
                            isSignifEmpPval = logical(0),
                            matching = character(0),
                            stringsAsFactors = FALSE      )
      
      if(length(FDR_signifTADs) > 0) {
        stopifnot( names(bp_nOverThresh_FDR) %in% names(adj_empPval))
        bp_dt1 <- data.frame(
          region = names(bp_nOverThresh_FDR),
          adjEmpPval = adj_empPval[names(bp_nOverThresh_FDR)],
          nOverThresh = bp_nOverThresh_FDR,
          isSignifFDR = names(bp_nOverThresh_FDR) %in% FDR_signifTADs,
          isSignifEmpPval = names(bp_nOverThresh_FDR) %in% adjEmpPval_signifTADs,
          matching = "overlapBp",
          stringsAsFactors = FALSE
        )
      } else {
        bp_dt1 <- emptyDT
      }
      if(length(adjEmpPval_signifTADs) > 0) {
        stopifnot( names(bp_nOverThresh_empPval) %in% names(adj_empPval))
        bp_dt2 <- data.frame(
          region = names(bp_nOverThresh_empPval),
          adjEmpPval = adj_empPval[names(bp_nOverThresh_empPval)],
          nOverThresh = bp_nOverThresh_empPval,
          isSignifFDR = names(bp_nOverThresh_empPval) %in% FDR_signifTADs,
          isSignifEmpPval = names(bp_nOverThresh_empPval) %in% adjEmpPval_signifTADs,
          matching = "overlapBp",
          stringsAsFactors = FALSE
        )
      } else{
        bp_dt2 <- emptyDT
      }
      bp_dt <- rbind(bp_dt1, bp_dt2)
      
      if(length(FDR_signifTADs) > 0) {
        stopifnot( names(genes_nOverThresh_FDR) %in% names(adj_empPval))
        ng_dt1 <- data.frame(
          region = names(genes_nOverThresh_FDR),
          adjEmpPval = adj_empPval[names(genes_nOverThresh_FDR)],
          nOverThresh = genes_nOverThresh_FDR,
          isSignifFDR = names(genes_nOverThresh_FDR) %in% FDR_signifTADs,
          isSignifEmpPval = names(genes_nOverThresh_FDR) %in% adjEmpPval_signifTADs,
          matching = "overlapGenes",
          stringsAsFactors = FALSE
        )
      } else {
        ng_dt1 <- emptyDT
      }
      if(length(adjEmpPval_signifTADs) > 0) {
        stopifnot( names(genes_nOverThresh_empPval) %in% names(adj_empPval))
        ng_dt2 <- data.frame(
          region = names(genes_nOverThresh_empPval),
          adjEmpPval = adj_empPval[names(genes_nOverThresh_empPval)],
          nOverThresh = genes_nOverThresh_empPval,
          isSignifFDR = names(genes_nOverThresh_empPval) %in% FDR_signifTADs,
          isSignifEmpPval = names(genes_nOverThresh_empPval) %in% adjEmpPval_signifTADs,
          matching = "overlapGenes",
          stringsAsFactors = FALSE
        )
      } else {
        ng_dt2 <- emptyDT
      }
      ng_dt <- rbind(ng_dt1, ng_dt2)
      
      nMatch_ratio_dt <- rbind(bp_dt, ng_dt)
      nMatch_ratio_dt$hicds <- hicds
      nMatch_ratio_dt$exprds <- exprds
      
      nTotTADs <- length(all_regs)
      nMatchDS <- length(unique(ds_matching_dt$matching_dataset))
      nSignifFDR <- length(FDR_signifTADs)
      nSignifEmpPval <- length(adjEmpPval_signifTADs)
      
      subTxtFDR <- paste0("nMatchDS=", nMatchDS, "; nSignifFDR=", nSignifFDR, "/", nTotTADs)
      subTxtPval <- paste0("nMatchDS=", nMatchDS, "; nSignifEmpPval=", nSignifEmpPval, "/", nTotTADs)
      subTxt <- paste0("nMatchDS=", nTotTADs, "; nSignifFDR=", nSignifFDR, "/", nTotTADs, "; nSignifEmpPval=", nSignifEmpPval, "/", nTotTADs)
      
      ot="overlapBp"
      st="isSignifFDR"
      subTxt1 <- paste0("pval.<=", pvalThresh, "; FDR <=", FDRthresh)
      subTxt2 <- paste0("pval.<=", pvalThresh, "; FDR <=", FDRthresh, "; match>=", matching_ratio_thresh)
      for(ot in c("overlapBp", "overlapGenes")) {
        for(st in c("isSignifFDR", "isSignifEmpPval")) {
          # if(st == "isSignifFDR") {
          #   subTxt <- subTxtFDR
          # } else if(st == "isSignifEmpPval"){
          #   subTxt <- subTxtPval
          # } else {
          #   stop("error")
          # }
          outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", ot, "Ratio_", st, ".", plotType))
          do.call(plotType, list(outFile, height=myHeight, width=myWidth))
          boxplot(as.formula(paste0(ot, "_ratio ~ ", st)), data=ds_matching_dt,
                  ylab = paste0(ot, "_ratio"),
                  xlab = paste0(st),
                  cex.lab = plotCex,
                  cex.axis  = plotCex,
                  main = paste0(hicds, " - ", exprds),
                  sub=subTxt1
          )
          mtext(side=3, text=subTxt)
          foo <- dev.off()
          cat(paste0("... written: ", outFile, "\n"))
          
          outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_nOverThresh_", ot, "_", st, ".", plotType))
          do.call(plotType, list(outFile, height=myHeight, width=myWidth))
          boxplot(as.formula(paste0("nOverThresh ~ ", st)), data = nMatch_ratio_dt[nMatch_ratio_dt$matching == paste0(ot),],
                  ylab = paste0("nOverThresh (", ot, ")"),
                  xlab = paste0(st),
                  cex.lab = plotCex,
                  cex.axis  = plotCex,
                  main = paste0(hicds, " - ", exprds),
                  sub=subTxt2
          )
          mtext(side=3, text=subTxt)        
          foo <- dev.off()
          cat(paste0("... written: ", outFile, "\n"))
        }
      }
      rownames(ds_matching_dt) <- NULL
      rownames(nMatch_ratio_dt) <- NULL
      list(
        # ds_matching_dt=ds_matching_dt[,c("overlapBp_ratio", "overlapGenes_ratio", "isSignifFDR", "isSignifEmpPval")],
        ds_matching_dt=ds_matching_dt,
        nMatch_ratio_dt=nMatch_ratio_dt
      )
      
    } # end-foreach iterating over exprds
    names(hicds_data) <-  all_exprds[[paste0(hicds)]]
    hicds_data
  } # end-foreach iterating over hicds  
  names(all_data) <- all_hicds
  
  outFile <- file.path(outFolder, "all_data.Rdata")
  save(all_data, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  subTxt1 <- paste0("pval.<=", pvalThresh, "; FDR <=", FDRthresh)
  subTxt2 <- paste0("pval.<=", pvalThresh, "; FDR <=", FDRthresh, "; match>=", matching_ratio_thresh)
  
  outFile <- file.path(outFolder, "all_data.Rdata")
  all_data <- get(load(outFile))
  
}
all_data <- unlist(all_data, recursive = FALSE)

all_ds_matching_dt <- do.call(rbind, lapply(all_data, function(x) x[["ds_matching_dt"]]))
all_nMatch_ratio_dt <- do.call(rbind, lapply(all_data, function(x) x[["nMatch_ratio_dt"]]))
rownames(all_ds_matching_dt) <- NULL
rownames(all_nMatch_ratio_dt) <- NULL

stopifnot(nrow(all_ds_matching_dt) == nrow(all_matching_dt))

totTADs <- nrow(all_ds_matching_dt)

nTotDS <- length(unique(all_matching_dt$ref_dataset))
stopifnot(length(all_data) == nTotDS)

subTxtAll <- paste0("nDS = ", nTotDS)

ot="overlapBp"
st="isSignifFDR"

for(ot in c("overlapBp", "overlapGenes")) {
  for(st in c("isSignifFDR", "isSignifEmpPval")) {
    
    nSignif <- sum(all_ds_matching_dt[,paste0(st)])
    
    subTxtAll_curr <- paste0(subTxtAll, "; signif. TADs = ", nSignif, "/", totTADs)
    
    outFile <- file.path(outFolder, paste0("allDS", "_", ot, "Ratio_", st, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    boxplot(as.formula(paste0(ot, "_ratio ~ ", st)), data=all_ds_matching_dt,
            ylab = paste0(ot, "_ratio"),
            xlab = paste0(st),
            cex.lab = plotCex,
            cex.axis  = plotCex,
            main = paste0("all datasets"),
            sub=subTxt1
    )
    mtext(side=3, text=subTxtAll_curr)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFolder, paste0("allDS", "_nOverThresh_", ot, "_", st, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    boxplot(as.formula(paste0("nOverThresh ~ ", st)), data = all_nMatch_ratio_dt[all_nMatch_ratio_dt$matching == paste0(ot),],
            ylab = paste0("nOverThresh (", ot, ")"),
            xlab = paste0(st),
            cex.lab = plotCex,
            cex.axis  = plotCex,
            main = paste0("all datasets"),
            sub=subTxt2
    )
    mtext(side=3, text=subTxtAll_curr)        
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
  
  
# ot = "overlapBp"
# for(ot in c("overlapBp", "overlapGenes")) {
  
  
  sub_dt <- all_nMatch_ratio_dt[all_nMatch_ratio_dt$matching == ot, c("region", "adjEmpPval", "nOverThresh", "hicds", "exprds")]
  sub_dt <- unique(sub_dt)
  stopifnot(nrow(sub_dt) > 0)
  all_ids <- paste0(sub_dt$hicds, sub_dt$exprds, sub_dt$region)
  stopifnot(!duplicated(all_ids))

  sub_dt$adjEmpPval_log10 <- -log10(sub_dt$adjEmpPval)
  
  x_var <- "adjEmpPval_log10"
  y_var <- "nOverThresh"
  myx <- sub_dt[,paste0(x_var)]
  myy <- sub_dt[,paste0(y_var)]
  

  x_thresh <- quantile(x=myx, pvalQt)
  y_thresh <- quantile(x=myy, conservQt)
  
  best_conserved_and_signif_dt <- sub_dt[which(sub_dt[,paste0(x_var)] >= x_thresh & sub_dt[,paste0(y_var)] >= y_thresh),]
  nPlotted <- nrow(best_conserved_and_signif_dt)
  nToPlot <- nPlotted
  
  outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_", ot, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x=(myx),
    y=(myy),
    xlab = paste0(x_var, " [-log10]"),
    ylab = paste0(y_var, ""),
    main = paste0(ot, " - ", y_var, " vs. ", x_var, ""),
    sub = paste0("# conserv. thresh.=", round(y_thresh, 2), " (", conservQt ," -qt); adj. emp. p-val. thresh.=", round(x_thresh, 2), " (", pvalQt, "-qt); n=", nPlotted),
    cex.lab=plotCex,
    cex.axis = plotCex
  )
  addCorr(x = myx, y = myy, bty="n", legPos = "bottomleft")
  mtext(side=3, text = paste0("nDS=", length(all_ids), "; match. thresh=", matching_ratio_thresh))
  
  text(x = best_conserved_and_signif_dt[,paste0(x_var)],
       y = best_conserved_and_signif_dt[,paste0(y_var)],
       labels = paste0(best_conserved_and_signif_dt[,paste0("hicds")],"\n", best_conserved_and_signif_dt[,paste0("exprds")],"-",best_conserved_and_signif_dt[,paste0("region")]),
       cex=0.4
       )

  if(nPlotted > 0) {
    
    plotList <- list()
    for(i_tad in 1:nPlotted) {
      plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = best_conserved_and_signif_dt$exprds[i_tad],
                                            hicds = best_conserved_and_signif_dt$hicds[i_tad],
                                            all_TADs = best_conserved_and_signif_dt$region[i_tad],
                                            orderByLolli = "startPos")
    } # end-for iterating over TADs to plot
    outFile <- file.path(outFolder, paste0("signif_and_conserv_", ot, "_TADs_lolliplot.", plotType ))
    mytit <- paste0(ot, " - # conserv. thresh.=", round(y_thresh, 2), " (", conservQt ," -qt); adj. emp. p-val. thresh.=", round(x_thresh, 2), " (", pvalQt, "-qt)")
    all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nToPlot == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
    outHeightGG <- min(c(7 * nPlotted/2, 49))
    outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
    outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
    
    ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
    cat("... written: ", outFile, "\n")
    
  
  
  }
  
  
  
  
  
  
  
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
    
}
  
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

