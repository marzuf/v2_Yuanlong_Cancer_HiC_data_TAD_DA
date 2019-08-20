# Rscript pval_biases.R

buildTable <- TRUE

script_name <- "pval_biases.R"

cat("> START: ", script_name, "\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")


library(foreach)
library(doMC)


source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("cancer_cols.R")


startTime <- Sys.time()

pointPch <- 16
pointCex <- 1
cexAxis <- 1.2
cexLab <- 1.2
plotCex <- 1.2


settingFolder <- file.path("PIPELINE", "INPUT_FILES")
pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(settingFolder))
stopifnot(dir.exists(pipOutFolder))

pvalThresh <- 0.05


script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10sameNbr_name <- "10sameNbr_runEmpPvalMeanTADCorr"
script11same_name <- "11sameNbr_runEmpPvalCombined"

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



outFolder <- file.path("PVAL_BIASES")
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)


if(buildTable) {
  
  
  all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      stopifnot(exprds %in% names(all_cmps))
      stopifnot(hicds %in% names(cl_cancer_annot))
      
      cat("... start building DT for :", hicds, " - ", exprds, "\n")
      
      ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
      script0_name <- "0_prepGeneData"
      stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
      geneListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneListFile))
      pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
      
      fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
      stopifnot(file.exists(fc_file))
      tad_fc <- eval(parse(text = load(fc_file)))
      all_regs <- names(tad_fc)
      
      
      corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(corr_file))  
      tad_corr <- eval(parse(text = load(corr_file)))
      stopifnot(setequal(all_regs, names(tad_corr)))
      
      fc_empPval_file <- file.path(pipOutFolder, hicds, exprds, script9_name, "emp_pval_meanLogFC.Rdata" )
      stopifnot(file.exists(fc_empPval_file))
      fc_empPval <- eval(parse(text = load(paste0(fc_empPval_file))))
      stopifnot(setequal(all_regs, names(fc_empPval)))
      fc_empPval <- fc_empPval[all_regs]
      # ADJUST THE PVAL
      tad_adjFCpval <- p.adjust(fc_empPval, method="BH")
      stopifnot(names(tad_adjFCpval) == all_regs)
      
      corr_empPval_file <- file.path(pipOutFolder, hicds, exprds, script10sameNbr_name, "emp_pval_meanCorr.Rdata" )
      stopifnot(file.exists(corr_empPval_file))
      corr_empPval <- eval(parse(text = load(paste0(corr_empPval_file))))
      stopifnot(setequal(all_regs, names(corr_empPval)))
      corr_empPval <- corr_empPval[all_regs]
      # ADJUST THE PVAL
      tad_adjCorrPval <- p.adjust(corr_empPval, method="BH")
      stopifnot(names(tad_adjCorrPval) == all_regs)
      
      comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      stopifnot(setequal(all_regs, names(comb_empPval)))
      comb_empPval <- comb_empPval[all_regs]
      # ADJUST THE PVAL
      tad_adjCombPval <- p.adjust(comb_empPval, method="BH")
      stopifnot(names(tad_adjCombPval) == all_regs)
      
      # PIPELINE/INPUT_FILES/Panc1_rep12_40kb/run_settings_TCGApaad_wt_mutKRAS.R
      settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
      stopifnot(file.exists(settingFile))
      cat("... source settingFile ", basename(settingFile), "\n")
      source(settingFile)
      cat("... load samp1\n")
      samp1 <- eval(parse(text = load(file.path(setDir, sample1_file))))
      cat("... load samp2\n")
      samp2 <- eval(parse(text = load(file.path(setDir, sample2_file))))
      
      
      data.frame(
        hicds=hicds,
        exprds=exprds,
        
        hicds_cancer = as.numeric(cl_cancer_annot[paste0(hicds)]),
        exprds_type = as.character(all_cmps[paste0(exprds)]),
        
        nGenes = length(pipeline_geneList),
        nTADs = length(all_regs),
        
        # nSamp1 = length(samp1),
        # nSamp2 = length(samp2),
        
        totSamp = length(samp1)+length(samp2),
        
        ratioSamp = length(samp1)/length(samp2),
        
        meanFC = mean(tad_fc),
        medianFC = median(tad_fc),
        meanCorr = mean(tad_corr),
        medianCorr = median(tad_corr),
        
        nSignif_meanCorrEmpPval = sum(tad_adjCorrPval <= pvalThresh),
        nSignif_meanFCempPval = sum(tad_adjFCpval <= pvalThresh),
        nSignif_combEmpPval = sum(tad_adjCombPval <= pvalThresh),
        
        stringsAsFactors = FALSE
      )
      
      
      
      
    } # end-foreach iterating exprds      
  } # end-foreach iterating hicds


  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  save(all_result_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,  "\n"))


} else {
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  all_result_dt <- eval(parse(text = load(outFile)))
}  
  
all_x <- c("totSamp", "ratioSamp", "nGenes", "nTADs")
stopifnot(all_x %in% colnames(all_result_dt))

all_y <- colnames(all_result_dt)[grepl("mean|median|nSignif", colnames(all_result_dt))]
stopifnot(length(all_y) > 0)

all_result_dt$exprds_type_col <- all_cols[all_result_dt$exprds_type]
all_result_dt$hicds_cancer_col <- ifelse(all_result_dt$hicds_cancer == 1, "red", 
                                         ifelse(all_result_dt$hicds_cancer == 0, "black", NA))
stopifnot(!is.na(all_result_dt$hicds_cancer_col))
stopifnot(!is.na(all_result_dt$exprds_type_col))

all_result_dt$dataset <- paste0(all_result_dt$hicds, "\n", all_result_dt$exprds)

x_var=all_x[1]
y_var=all_y[1]

foo <- foreach(y_var = all_y) %dopar% {
  
  myy <- all_result_dt[,paste0(y_var)]
  
  for(x_var in all_x) {
    
  
    myx <- all_result_dt[,paste0(x_var)]  

    outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_colByExprds_and_colByHicds.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))    
    par(mfrow=c(1,2))
    # outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_colByExprds.", plotType))
    # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(
      x=myx,
      y=myy,
      xlab=paste0(x_var),
      ylab=paste0(y_var),
      main=paste0(y_var, " vs. ", x_var),
      pch = 16,
      cex=0.7,
      cex.lab=plotCex,
      cex.axis=plotCex,
      col=all_result_dt$exprds_type_col
    )
    text(x=myx, y=myy, cex=0.5, labels=all_result_dt$dataset, col=all_result_dt$exprds_type_col)
    addCorr(x=myx, y=myy, bty="n")
    addSubtypeLeg(bty="n", pch=16)
    # foo <- dev.off()
    # cat(paste0("... written: ", outFile, "\n"))
    
    # outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_colByHicds.", plotType))
    # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(
      x=myx,
      y=myy,
      xlab=paste0(x_var),
      ylab=paste0(y_var),
      main=paste0(y_var, " vs. ", x_var),
      pch = 16,
      cex=0.7,
      cex.lab=plotCex,
      cex.axis=plotCex,
      col=all_result_dt$hicds_cancer_col
    )
    addCorr(x=myx, y=myy, bty="n")
    text(x=myx, y=myy, cex=0.5, labels=all_result_dt$dataset, col=all_result_dt$hicds_cancer_coll)
    legend("bottomleft", c("cancer", "not cancer"), 
           pch=16, col=c("red", "black"), bty="n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  } # end-for iterating all_x
} # end-foreach iterating all_y







#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))



