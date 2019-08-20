options(scipen=100)

setDir <- ""

# Rscript cmp_same_cancer_notCancer.R 

script_name <- "cmp_same_cancer_notCancer.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(reshape2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("cancer_cols.R")

script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight*1.2  # for density plot
axisCex <- 1.4
labCex <- 1.4

tiesMeth <- "min"


pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CMP_SAME_CANCER_NOTCANCER"
dir.create(outFolder, recursive=TRUE)

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


# retrieve the "duplicated" all_exprds
all_exprds_vect <- unlist(all_exprds)
dup_exprds <- unique(all_exprds_vect[duplicated(all_exprds_vect)])

expr_hicds <- sapply(dup_exprds, function(exprds) {
  names(all_exprds)[sapply(all_exprds, function(x) exprds %in% x )]
})

stopifnot(dup_exprds %in% names(expr_hicds))

stopifnot(unlist(expr_hicds) %in% names(cl_cancer_annot))

### BUILD SIGNIF ALONG FDR THRESH
exprds = dup_exprds[1]
foo <- foreach(exprds = dup_exprds) %dopar% {
  
  all_hicds <- expr_hicds[[paste0(exprds)]]

  if(all(cl_cancer_annot[all_hicds] == 1) | all(cl_cancer_annot[all_hicds] == 0)) return(NULL)
  
  
  sameExpr_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    
    pvalFile <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
    stopifnot(file.exists(pvalFile))
    combPvals <- eval(parse(text = load(pvalFile)))
    adj_combPvals <- p.adjust(combPvals, method="BH")
    adj_combPvals_log10 <- -log10(adj_combPvals)
    
    fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
    corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
    stopifnot(file.exists(fc_file))
    stopifnot(file.exists(corr_file))  
    
    tad_logFC <- get(load(fc_file))
    tad_meanCorr <- get(load(corr_file))
    all_regs <- names(tad_logFC)
    
    stopifnot(length(adj_combPvals_log10) == length(tad_logFC))
    stopifnot(length(adj_combPvals_log10) == length(tad_meanCorr))
    
    stopifnot(setequal(all_regs, names(tad_meanCorr)))
    stopifnot(setequal(all_regs, names(adj_combPvals_log10)))
    
    data.frame(
      hicds = hicds,
      exprds=exprds,
      isCancer=as.numeric(cl_cancer_annot[hicds]),
      adj_combPvals_log10=adj_combPvals_log10[all_regs],
      logFC = tad_logFC[all_regs],
      meanCorr = tad_meanCorr[all_regs],
      stringsAsFactors = FALSE
    )
    
  } # end-foreach iterating over hicds of current exprds
  all_vars <- c("logFC", "meanCorr", "adj_combPvals_log10")
  stopifnot(all_vars %in% colnames(sameExpr_dt))
  for(plot_var in all_vars) {
    
    outFile <- file.path(outFolder, paste0(exprds, "_", plot_var, "_cancer_notCancer_mutlidens.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens(
      split(sameExpr_dt[, paste0(plot_var)], sameExpr_dt$isCancer),
      plotTit = paste0(exprds, " - ", plot_var)
    )
    mtext(side=3, text=paste0(all_hicds, collapse=","))
    legend("topleft", c("0 = not cancer", "1 = cancer" ), bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
  } # end-for iterating over plotting variables
  
} # end-foreach iterating over exprds


##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



