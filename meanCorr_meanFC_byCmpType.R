# compare the FDR~corr thresholds and FDR~FC thresholds

setDir=""

options(scipen=100)

# Rscript meanCorr_meanFC_byCmpType.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript meanCorr_meanFC_byCmpType.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "meanCorr_meanFC_byCmpType.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
axisCex <- 1.4

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))


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

outFolder <- "MEANCORR_MEANFC_BYCMPTYPE"
dir.create(outFolder, recursive=TRUE)

buildTable <- TRUE

if(buildTable) {

  all_data_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  
  exprds_data_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
    corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
if(!file.exists(fc_file))cat(fc_file,"\n")
if(!file.exists(corr_file))cat(corr_file,"\n")
    stopifnot(file.exists(fc_file))
    stopifnot(file.exists(corr_file))  
    
    tad_fc <- eval(parse(text = load(fc_file)))
    tad_corr <- eval(parse(text = load(corr_file)))
    
    all_regs <- names(tad_fc)
    stopifnot(setequal(all_regs, names(tad_corr)))
    
    stopifnot(length(tad_fc) == length(all_regs))
    stopifnot(length(tad_corr) == length(all_regs))
    
    plotDT <- data.frame(
      hicds = hicds,
      exprds=exprds,
      region = all_regs,
      meanFC = as.numeric(tad_fc[all_regs]),
      meanCorr = as.numeric(tad_corr[all_regs]),
      stringsAsFactors = FALSE
    )
    plotDT
    
  } # end iterating over exprds
  exprds_data_dt
  } # end iterating over hicds
  
  
  outFile <- file.path(outFolder, "all_data_dt.Rdata")
  save(all_data_dt, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))


} else {  #end-if buildTable
  outFile <- file.path(outFolder, "all_data_dt.Rdata")
  all_data_dt <- eval(parse(text = load(outFile)))
  
}


all_data_dt$subtype_col <- all_cols[paste0(all_cmps[paste0(all_data_dt$exprds)])]
all_data_dt$dataset <- paste0( all_data_dt$hicds, " - ", all_data_dt$exprds)
all_datasets <- unique(all_data_dt$dataset) 
nDS <- length(unique( all_datasets ))

###################################
################################### PLOT EACH DATASET SEPARATELY
###################################

foo <- foreach(ds = all_datasets) %dopar% {
  
  sub_dt <- all_data_dt[all_data_dt$dataset == ds,]
  
  x_var = "meanFC"
  y_var = "meanCorr"
  myx <- sub_dt[,paste0(x_var)]
  myy <- sub_dt[,paste0(y_var)]
  outFile <- file.path(outFolder, paste0(gsub(" - ", "_", ds), "_", y_var, "_vs_", x_var, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = myx,
    y = myy,
    cex.axis = axisCex,
    cex.lab = axisCex,
    xlab = paste0(gsub("_", " ", x_var)),
    ylab = paste0(gsub("_", " " , y_var)),
    main = paste0(y_var, " vs. ", x_var)
  )
  mtext(side=3, text = paste0(ds, " (n=", nrow(sub_dt),")"), font=3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  
}


###################################
################################### PLOT FOR ALL DATASETS
###################################


x_var = "meanFC"
y_var = "meanCorr"
myx <- all_data_dt[,paste0(x_var)]
myy <- all_data_dt[,paste0(y_var)]
outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  cex.axis = axisCex,
  cex.lab = axisCex,
  xlab = paste0(gsub("_", " ", x_var)),
  ylab = paste0(gsub("_", " " , y_var)),
  main = paste0(y_var, " vs. ", x_var)
)
mtext(side=3, text = paste0(nDS, " DS (n=", nrow(all_data_dt),")"), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


x_var = "meanFC"
y_var = "meanCorr"
myx <- all_data_dt[,paste0(x_var)]
myy <- all_data_dt[,paste0(y_var)]
outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_plotColType.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = myx,
  y = myy,
  cex.axis = axisCex,
  cex.lab = axisCex,
  xlab = paste0(gsub("_", " ", x_var)),
  ylab = paste0(gsub("_", " " , y_var)),
  main = paste0(y_var, " vs. ", x_var),
  col = all_data_dt$subtype_col,
  pch=16,
  cex=0.7
)
addSubtypeLeg(mypos="bottomright", bty="n")
mtext(side=3, text = paste0(nDS, " DS (n=", nrow(all_data_dt),")"), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))








##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
