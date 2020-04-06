# Rscript cmp_pval_byCmp.R

outFolder <- "CMP_PVAL_BYCMP"
dir.create(outFolder)


plotType <- "png"
myHeightGG <- 7
myWidthGG <- 9
plotCex <- 1.4

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")


require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(ggpubr)

buildTable <- TRUE

script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"

script9_name <- "9_runEmpPvalMeanTADLogFC"
script9sameNbr_name <- "9sameNbr_runEmpPvalMeanTADLogFC"
script9random_name <- "9random_runEmpPvalMeanTADLogFC"
script9randomResc_name <- "9randomRescaled_runEmpPvalMeanTADLogFC"
script10_name <- "10sameNbr_runEmpPvalMeanTADCorr"
script10random_name <- "10random_runEmpPvalMeanTADCorr"

pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

hicds="LG1_40kb"
exprds="TCGAluad_norm_luad"

var_dt <- get(load("../OLDER_v2_Yuanlong_Cancer_HiC_data_TAD_DA/EXPR_VARIANCE/LOG2FPKM/all_ds_geneVarDT.Rdata"))

if(buildTable){
  
  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      meanFC <- get(load(file.path(pipFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")))
      meanCorr <- get(load(file.path(pipFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")))
      
      logfc_pval_permG2T <- get(load(file.path(pipFolder, hicds, exprds, script9_name, "emp_pval_meanLogFC.Rdata")))
      logfc_pval_sameNbr <- get(load(file.path(pipFolder, hicds, exprds, script9sameNbr_name, "emp_pval_meanLogFC.Rdata")))
      logfc_pval_random <- get(load(file.path(pipFolder, hicds, exprds, script9random_name, "emp_pval_meanLogFC.Rdata")))
      logfc_pval_randomResc <- get(load(file.path(pipFolder, hicds, exprds, script9randomResc_name, "emp_pval_meanLogFC.Rdata")))
      
      corr_pval_sameNbr <- get(load(file.path(pipFolder, hicds, exprds, script10_name, "emp_pval_meanCorr.Rdata")))
      corr_pval_random <- get(load(file.path(pipFolder, hicds, exprds, script10random_name, "emp_pval_meanCorr.Rdata")))
      
      stopifnot(length(logfc_pval_permG2T) == length(corr_pval_sameNbr))
      
      medianMostVar <- var_dt$medianMostVar[var_dt$hicds == hicds & var_dt$exprds == exprds]
      
      data.frame(
        hicds = hicds,
        exprds=exprds,
        meanFC = meanFC,
        meanCorr = meanCorr,
        
        logfc_pval_permG2T = logfc_pval_permG2T,
        logfc_pval_sameNbr = logfc_pval_sameNbr,
        logfc_pval_random = logfc_pval_random,
        logfc_pval_randomResc = logfc_pval_randomResc,
        
        meanCorr_pval_sameNbr = corr_pval_sameNbr,
        meanCorr_pval_random = corr_pval_random,
        
        stringsAsFactors = FALSE
      )
      
    }
    exprds_dt
  }
  outFile <- file.path(outFolder, "all_dt.Rdata")
  save(all_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_dt.Rdata")
  all_dt <- get(load(outFile))
}

curr_var="medianMostVar"
var_dt$cmp <- all_cmps[var_dt$exprds]
stopifnot(!is.na(var_dt$cmp))
p <- ggboxplot(
  var_dt, x = "cmp", y = paste0(curr_var),
  xlab = "", ylab=curr_var, title = paste0(curr_var)
) + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(size=16, face = "bold", hjust=0.5)  )
outFile <- file.path(outFolder, paste0(curr_var, "_by_cmpType.", plotType))
ggsave(p, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



load("CMP_PVAL_BYCMP/all_dt.Rdata")

all_vars <- colnames(all_dt)
all_vars <- all_vars[! all_vars %in% c("hicds", "exprds")]

all_dt$cmp <- all_cmps[all_dt$exprds]
stopifnot(!is.na(all_dt$cmp))

curr_var="logfc_pval_permG2T"

for(curr_var in all_vars) {
  
  
  p <- ggboxplot(
    all_dt, x = "cmp", y = paste0(curr_var),
    xlab = "", ylab=curr_var, title = paste0(curr_var)
  )+ 
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16, face = "bold", hjust=0.5)
    )
  outFile <- file.path(outFolder, paste0(curr_var, "_by_cmpType.", plotType))
  ggsave(p, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  if(grepl("pval", curr_var)) {
    
    curr_var_log10 <- paste0(curr_var, "_log10")
    all_dt[,curr_var_log10] <- -log10(all_dt[, paste0(curr_var)])
    
    p_log10 <- ggboxplot(
      all_dt, x = "cmp", y = paste0(curr_var_log10),
      xlab = "", ylab = curr_var_log10, title = paste0(curr_var_log10)
    )+ 
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=14),
            plot.title = element_text(size=16, face = "bold", hjust=0.5)      )
    
    outFile <- file.path(outFolder, paste0(curr_var_log10, "_by_cmpType.", plotType))
    ggsave(p_log10, file=outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  }

}
