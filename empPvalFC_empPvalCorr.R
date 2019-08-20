# Rscript empPvalFC_empPvalCorr.R

script_name <- "empPvalFC_empPvalCorr.R"
cat("> Start ", script_name, "\n")
startTime <- Sys.time()

require(foreach)
require(doMC)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

registerDoMC(40)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight

cexPlot <- 1.2


outFold <- "EMPPVALFC_EMPPVALCORR"
dir.create(outFold)




pipFolder <- file.path(".")
pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))




#script9_name <- "910000_runEmpPvalMeanTADLogFC" # => EMPPVALFC_EMPPVALCORR_10000
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10sameNbr_runEmpPvalMeanTADCorr"

all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanLogFC.Rdata", full.names = FALSE)
all_fc_files <- all_fc_files[grep(script9_name, all_fc_files)]
stopifnot(length(all_fc_files) > 0)

all_meanCorr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_meanCorr.Rdata", full.names = FALSE)
all_meanCorr_files <- all_meanCorr_files[grep(script10_name, all_meanCorr_files)]
stopifnot(length(all_meanCorr_files) > 0)

stopifnot(length(all_meanCorr_files) == length(all_fc_files))




### BUILD THE LOGFC TABLE
cat("... start build fc_DT \n")
fc_file = all_fc_files[1]
fc_DT <- foreach(fc_file = all_fc_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, fc_file)
  stopifnot(file.exists(curr_file))
  tad_fc <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(fc_file))
  hicds <- dirname(dataset)
  exprds <- basename(dataset)
  data.frame(
    dataset = dataset,
    hicds = hicds,
    exprds=exprds,
    region = names(tad_fc),
    empPval_meanFC = as.numeric(tad_fc),
    stringsAsFactors = FALSE
  )
}
### BUILD THE MEANCORR TABLE
cat("... start build meanCorr_DT \n")
meanCorr_file = all_meanCorr_files[1]
meanCorr_DT <- foreach(meanCorr_file = all_meanCorr_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, meanCorr_file)
  stopifnot(file.exists(curr_file))
  tad_meanCorr <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(meanCorr_file))
  hicds <- dirname(dataset)
  exprds <- basename(dataset)
  data.frame(
    dataset = dataset,
    hicds = hicds,
    exprds=exprds,
    region = names(tad_meanCorr),
    empPval_meanCorr = as.numeric(tad_meanCorr),
    stringsAsFactors = FALSE
  )
}




stopifnot(nrow(meanCorr_DT) == nrow(fc_DT))




id_cols <-  c("dataset", "hicds", "exprds", "region")

all_dt <- merge(meanCorr_DT, fc_DT, by = id_cols, all = TRUE)
stopifnot(nrow(all_dt) == nrow(fc_DT))




outFile <- file.path(outFold, "all_dt.Rdata")
save(all_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


nDS <- length(unique(all_dt$dataset))
nTADs <- nrow(all_dt)



#############################################################################################################################
############################################################################################################################# empPvalMeanCorr vs. empPvalMeanFC
#############################################################################################################################

all_x <- c("empPval_meanFC")

all_y <- c("empPval_meanCorr")

for(x_var in all_x) {
  
  stopifnot(x_var %in% colnames(all_dt))
  
  myx <- all_dt[,paste0(x_var)]
  
  foo <- foreach(y_var = all_y) %dopar% {
    
    
    
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = myx,
      y = all_dt[,paste0(y_var)],
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var),
      ylab = paste0(y_var),
      main = paste0(y_var, " vs.\n", x_var)
    )
    mtext(side = 3, text = paste0("nDS = ", nDS, " - nTADs = ", nTADs), font = 3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    outFile <- file.path(outFold, paste0(y_var, "_vs_", x_var, "_log10.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = - log10(myx),
      y = - log10(all_dt[,paste0(y_var)]),
      cex.axis = cexPlot,
      cex.lab = cexPlot,
      xlab = paste0(x_var, " [-log10]"),
      ylab = paste0(y_var, " [-log10]"),
      main = paste0(y_var, " vs.\n", x_var)
    )
    mtext(side = 3, text = paste0("nDS = ", nDS, " - nTADs = ", nTADs), font = 3)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
  
  
}

#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))
