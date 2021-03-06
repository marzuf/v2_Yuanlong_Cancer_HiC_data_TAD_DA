# Rscript pval_biases_TADlevel.R

buildTable <- FALSE

script_name <- "pval_biases_TADlevel.R"

cat("> START: ", script_name, "\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

library(foreach)
library(doMC)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("cancer_cols.R")


startTime <- Sys.time()

pointPch <- 16
pointCex <- 1
cexAxis <- 1.2
cexLab <- 1.2
plotCex <- 1.2

pipFolder <- "."
settingFolder <- file.path("PIPELINE", "INPUT_FILES")
pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(settingFolder))
stopifnot(dir.exists(pipOutFolder))

pvalThresh <- 0.05


script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script8cRatioDown_name <- "8cOnlyRatioDownFastSave_runAllDown"
script8cFCC_name <- "8cOnlyFCC_runAllDown" 
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

outFolder <- file.path("PVAL_BIASES_TADLEVEL")
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)


"/LOG2FPKM/"

varFile <- file.path("EXPR_VARIANCE_BYTAD/LOG2FPKM/all_ds_geneVarDT.Rdata")
stopifnot(file.exists(varFile))
varData <- get(load(varFile))
var_variable <- "tadMeanVar"

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
      
      
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
      g2t_DT <- g2t_DT[ g2t_DT$entrezID %in% pipeline_geneList,]
      gt <- table(g2t_DT$region)
      tad_genes <- setNames(as.numeric(gt), names(gt))
      
      
      
      
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
      
      
      rd_file <-  file.path(pipOutFolder, hicds, exprds, script8cRatioDown_name, "all_obs_ratioDown.Rdata")
      stopifnot(file.exists(rd_file))
      tad_rd <- eval(parse(text = load(rd_file)))
      stopifnot(names(tad_rd) == all_regs)
      
      
      fcc_file <-  file.path(pipOutFolder, hicds, exprds, script8cFCC_name, "all_obs_prodSignedRatio.Rdata")
      stopifnot(file.exists(fcc_file))
      tad_fcc <- eval(parse(text = load(fcc_file)))
      stopifnot(names(tad_fcc) == all_regs)
      
      stopifnot(paste0(hicds, "_", exprds) %in% names(varData))
      curr_varData <- varData[[paste0(hicds, "_", exprds)]]
      stopifnot(var_variable %in% names(curr_varData))
      tad_var <- curr_varData[[paste0(var_variable)]]
      stopifnot(setequal(names(tad_var), all_regs))
      
      # PIPELINE/INPUT_FILES/Panc1_rep12_40kb/run_settings_TCGApaad_wt_mutKRAS.R
      settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
      stopifnot(file.exists(settingFile))
      cat("... source settingFile ", basename(settingFile), "\n")
      source(settingFile)
      cat("... load samp1\n")
      samp1 <- eval(parse(text = load(file.path(setDir, sample1_file))))
      cat("... load samp2\n")
      samp2 <- eval(parse(text = load(file.path(setDir, sample2_file))))
      

      stopifnot(names(tad_genes) == all_regs)
      
      data.frame(
        hicds=hicds,
        exprds=exprds,
        
        hicds_cancer = as.numeric(cl_cancer_annot[paste0(hicds)]),
        exprds_type = as.character(all_cmps[paste0(exprds)]),
        
        
        nGenes = tad_genes[all_regs],
        
        # nSamp1 = length(samp1),
        # nSamp2 = length(samp2),
        
        totSamp = length(samp1)+length(samp2),
        
        ratioSamp = length(samp1)/length(samp2),
        
        
        logFC = as.numeric(tad_fc[all_regs]),
        
        meanCorr = as.numeric(tad_corr[all_regs]),
        
        
        meanVar = as.numeric(tad_var[all_regs]),
        
        adj_empPval_meanCorr = as.numeric(tad_adjCorrPval[all_regs]),
        adj_empPval_meanFC = as.numeric(tad_adjFCpval[all_regs]),
        adj_combEmpPval = as.numeric(tad_adjCombPval[all_regs]),
        
        
        
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
  



all_result_dt$meanVar_log10 <- log10(all_result_dt$meanVar)
x_var <- "meanCorr"
y_var <- "meanVar_log10"
 
all_result_dt$exprds_type_col <- all_cols[all_result_dt$exprds_type]
all_result_dt$hicds_cancer_col <- ifelse(all_result_dt$hicds_cancer == 1, "red",
                                         ifelse(all_result_dt$hicds_cancer == 0, "black", NA))
myy <- all_result_dt[,paste0(y_var)]
myx <- all_result_dt[,paste0(x_var)]  

outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_colByExprds_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x=myx,
  y=myy,
  xlab=paste0(x_var, "[log10]"),
  ylab=paste0(y_var, "[log10]"),
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
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



x_var <- "logFC"
y_var <- "meanVar_log10"

all_result_dt$exprds_type_col <- all_cols[all_result_dt$exprds_type]
all_result_dt$hicds_cancer_col <- ifelse(all_result_dt$hicds_cancer == 1, "red",
                                         ifelse(all_result_dt$hicds_cancer == 0, "black", NA))
myy <- all_result_dt[,paste0(y_var)]
myx <- all_result_dt[,paste0(x_var)]  

outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_colByExprds_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x=myx,
  y=myy,
  xlab=paste0(x_var, "[log10]"),
  ylab=paste0(y_var, "[log10]"),
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
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
stop("-ok")

































all_vars <- colnames(all_result_dt)[! colnames(all_result_dt) %in% c( "hicds", "exprds", "hicds_cancer", "exprds_type")]

all_cmbs <- combn(all_vars, 2)

for(i in 1:ncol(all_cmbs)){
  
  x_var <- all_cmbs[1,i]
  y_var <- all_cmbs[2,i]


# all_y <- colnames(all_result_dt)[grepl("adj_", colnames(all_result_dt))]
# stopifnot(length(all_y) > 0)
# 
# all_x <- colnames(all_result_dt)[! colnames(all_result_dt) %in% c(all_y, "hicds", "exprds", "hicds_cancer", "exprds_type")]
# 
all_result_dt$exprds_type_col <- all_cols[all_result_dt$exprds_type]
all_result_dt$hicds_cancer_col <- ifelse(all_result_dt$hicds_cancer == 1, "red",
                                         ifelse(all_result_dt$hicds_cancer == 0, "black", NA))
# stopifnot(!is.na(all_result_dt$hicds_cancer_col))
# stopifnot(!is.na(all_result_dt$exprds_type_col))
# 
# all_result_dt$dataset <- paste0(all_result_dt$hicds, "\n", all_result_dt$exprds)
# 
# x_var=all_x[1]
# y_var=all_y[1]

# foo <- foreach(y_var = all_y) %dopar% {
#   
#   
#   for(x_var in all_x) {
#     
    myy <- all_result_dt[,paste0(y_var)]
    
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
    text(x=myx, y=myy, cex=0.5, labels=all_result_dt$dataset, col=all_result_dt$hicds_cancer_col)
    legend("bottomleft", c("cancer", "not cancer"), 
           pch=16, col=c("red", "black"), bty="n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
}
























#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))



