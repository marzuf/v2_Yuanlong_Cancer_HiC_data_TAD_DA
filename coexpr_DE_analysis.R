# Rscript coexpr_DE_analysis.R

options(save.defaults = list(version = 2))

script_name <- "coexpr_DE_analysis.R"

startTime <- Sys.time()

cat("> START coexpr_DE_analysis.R \n")

SSHFS <- FALSE

require(tools)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

options(scipen=100)


outFolder <- "COEXPR_DE_ANALYSIS"
dir.create(outFolder, recursive=TRUE)


pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script8cRatioDown_name <- "8cOnlyRatioDownFastSave_runAllDown"
script8cFCC_name <- "8cOnlyFCC_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script9v2_name <- "9v2_runEmpPvalWilcoxStat"
script10_name <- "10_runEmpPvalMeanTADCorr"
script10v2_name <- "10v2_runEmpPvalMeanTADCorr"
script10b_name <- "10b_runEmpPvalProdSignedRatio"
script11_name <- "11sameNbr_runEmpPvalCombined"


all_fcc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_obs_prodSignedRatio.Rdata", full.names = FALSE)
all_fcc_files <- all_fcc_files[grepl(script8cFCC_name, all_fcc_files)]
stopifnot(length(all_fcc_files) > 0)

all_corr_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanCorr_TAD.Rdata", full.names = FALSE)
all_corr_files <- all_corr_files[grepl(script4_name, all_corr_files)]
stopifnot(length(all_corr_files) > 0)


all_fc_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_meanLogFC_TAD.Rdata", full.names = FALSE)
all_fc_files <- all_fc_files[grepl(script3_name, all_fc_files)]
stopifnot(length(all_fc_files) > 0)

all_pvalcomb_files <- list.files(pipOutFolder, recursive = TRUE, pattern="emp_pval_combined.Rdata", full.names = FALSE)
all_pvalcomb_files <- all_pvalcomb_files[grepl(script11_name, all_pvalcomb_files)]
stopifnot(length(all_pvalcomb_files) > 0)

all_ratioDown_files <- list.files(pipOutFolder, recursive = TRUE, pattern="all_obs_ratioDown.Rdata", full.names = FALSE)
all_ratioDown_files <- all_ratioDown_files[grepl(script8cRatioDown_name, all_ratioDown_files)]
stopifnot(length(all_ratioDown_files) > 0)

stopifnot(length(all_ratioDown_files) == length(all_fc_files))
stopifnot(length(all_fc_files) == length(all_pvalcomb_files) )
stopifnot(length(all_fc_files) == length(all_corr_files) )
stopifnot(length(all_fc_files) == length(all_fcc_files) )
  

### BUILD THE ratio down TABLE
rd_file = all_ratioDown_files[1]
rD_DT <- foreach(rd_file = all_ratioDown_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder,rd_file)
  stopifnot(file.exists(curr_file))
  tad_rd <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(rd_file))
  data.frame(
    dataset = dataset,
    region = names(tad_rd),
    ratioDown = as.numeric(tad_rd),
    stringsAsFactors = FALSE
  )
}

#### BUILD THE corr TABLE
corr_file = all_corr_files[1]
corr_DT <- foreach(corr_file = all_corr_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, corr_file)
  stopifnot(file.exists(curr_file))
  tad_corr <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(corr_file))
  data.frame(
    dataset = dataset,
    region = names(tad_corr),
    meanCorr = as.numeric(tad_corr),
    stringsAsFactors = FALSE
  )
}

### BUILD THE FCC TABLE

fcc_file = all_fcc_files[1]
fcc_DT <- foreach(fcc_file = all_fcc_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder,fcc_file)
  stopifnot(file.exists(curr_file))
  tad_fcc <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(fcc_file))
  data.frame(
    dataset = dataset,
    region = names(tad_fcc),
    FCC = as.numeric(tad_fcc),
    stringsAsFactors = FALSE
  )
}



### BUILD THE LOGFC TABLE
fc_file = all_fc_files[1]
fc_DT <- foreach(fc_file = all_fc_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, fc_file)
  stopifnot(file.exists(curr_file))
  tad_fc <- eval(parse(text = load(curr_file)))
  dataset <- dirname(dirname(fc_file))
  data.frame(
    dataset = dataset,
    region = names(tad_fc),
    meanFC = as.numeric(tad_fc),
    stringsAsFactors = FALSE
  )
}


### BUILD THE LOGFC TABLE
pvalcomb_file = all_pvalcomb_files[1]
pvalcomb_DT <- foreach(pvalcomb_file = all_pvalcomb_files, .combine = 'rbind') %dopar% {
  curr_file <- file.path(pipOutFolder, pvalcomb_file)
  stopifnot(file.exists(curr_file))
  tad_pvalcomb <- eval(parse(text = load(curr_file)))
  # adj pval comb
  adj_tad_pvalcomb <- p.adjust(tad_pvalcomb, method="BH")
  dataset <- dirname(dirname(pvalcomb_file))
  data.frame(
    dataset = dataset,
    region = names(adj_tad_pvalcomb),
    adjPvalComb = as.numeric(adj_tad_pvalcomb),
    stringsAsFactors = FALSE
  )
}


stopifnot(nrow(pvalcomb_DT) == nrow(fc_DT) )
stopifnot(nrow(pvalcomb_DT) == nrow(rD_DT) )
stopifnot(nrow(pvalcomb_DT) == nrow(corr_DT) )
stopifnot(nrow(pvalcomb_DT) == nrow(fcc_DT) )



############################################################################## MERGE THE TABLES


tad_coexpr_fc_DT <- merge(fcc_DT, 
                        merge(corr_DT, 
                            merge(rD_DT, 
                                merge(pvalcomb_DT, fc_DT, by=c("dataset", "region")),
                            by=c("dataset", "region")), 
                        by=c("dataset", "region")), 
                    by=c("dataset", "region"))

cat(nrow(pvalcomb_DT), "\n")
cat(nrow(tad_coexpr_fc_DT), "\n")

stopifnot(nrow(pvalcomb_DT) == nrow(tad_coexpr_fc_DT))
stopifnot(!is.na(tad_coexpr_fc_DT))


tad_coexpr_fc_DT$adjPvalComb_log10 <- -log10(tad_coexpr_fc_DT$adjPvalComb)


all_vars <- colnames(tad_coexpr_fc_DT)
all_vars <- all_vars[! all_vars %in% c("adjPvalComb", "dataset", "region")]

all_comb <- combn(all_vars, 2)


# for each, plot i) densplot; ii) plot with color for subtypes
tad_coexpr_fc_DT$cmps <- basename(tad_coexpr_fc_DT$dataset)

colDT <- data.frame(
  cmps = names(all_cmps),
  cmpType = all_cmps,
stringsAsFactors = FALSE
)
colDT$cmpCol <- all_cols[colDT$cmpType]
stopifnot(!is.na(colDT))

stopifnot(tad_coexpr_fc_DT$cmps %in% colDT$cmps)

tad_coexpr_fc_DT <- merge(tad_coexpr_fc_DT, colDT, by = "cmps", all.x = TRUE, all.y = FALSE)

stopifnot(!is.na(tad_coexpr_fc_DT$cmpCol))


outFile <- file.path(outFolder, "tad_coexpr_fc_DT.Rdata")
save(tad_coexpr_fc_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


myplot_densplot <- function(xvar, yvar, addCurve=FALSE, dt = tad_coexpr_fc_DT, outPrefix="", saveFile=TRUE) {
  
  stopifnot(xvar %in% colnames(dt))
  myx <- dt[,xvar]
  stopifnot(yvar %in% colnames(dt))
  myy <- dt[,yvar]
  mycols <- dt$cmpCol
  
  outFile <- file.path(outFolder, paste0(outPrefix, "", yvar, "_vs_", xvar, "_densplot.", plotType))
  if(saveFile) do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(y=myy,
           x=myx,
           pch=16,
           cex=0.7,
           cex.axis = myCexAxis,
           cex.lab = myCexLab,
           xlab=xvar,
           ylab=yvar,
           main = paste0(yvar, " vs. ", xvar)
  )
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  if(addCurve) {
    curve(1*x, lty=2, col="grey", add = TRUE)
  }
  if(saveFile) foo <- dev.off()
  if(saveFile) cat(paste0("... written: ", outFile, "\n"))
}

myplot_colplot <- function(xvar, yvar, mycols, addCurve = FALSE, dt = tad_coexpr_fc_DT, outPrefix="",saveFile=TRUE) {
  
  stopifnot(xvar %in% colnames(dt))
  myx <- dt[,xvar]
  stopifnot(yvar %in% colnames(dt))
  myy <- dt[,yvar]
  
  outFile <- file.path(outFolder, paste0(outPrefix, "", yvar, "_vs_", xvar, "_colplot.", plotType))
  if(saveFile) do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(y=myy,
           x=myx,
           pch=16,
           cex=0.7,
           cex.axis = myCexAxis,
           cex.lab = myCexLab,
           xlab=xvar,
           ylab=yvar,
       col=mycols,
           main = paste0(yvar, " vs. ", xvar)
  )
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  if(addCurve) {
    curve(1*x, lty=2, col="grey",add=TRUE)
  }
  addSubtypeLeg(bty="n")
  if(saveFile) foo <- dev.off()
  if(saveFile) cat(paste0("... written: ", outFile, "\n"))
}

mycols <- tad_coexpr_fc_DT$cmpCol

##########################
### plot all possible pairs of variables
##########################


i_comp=1
for(i_comp in 1:ncol(all_comb)) {

  xvar <- all_comb[1, i_comp]
  yvar <- all_comb[2, i_comp]

  myplot_densplot(xvar,yvar, saveFile=F)
  myplot_colplot(xvar,yvar,mycols, saveFile=F)
  myplot_densplot(xvar,yvar)
  myplot_colplot(xvar,yvar,mycols)


}







# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





