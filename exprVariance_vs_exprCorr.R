# Rscript exprVariance_vs_exprCorr.R 

startTime <- Sys.time()
options(scipen=100)

cat("> START: exprVariance_vs_exprCorr.R\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

library(foreach)
library(doMC)
library(tools)
require(dplyr)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("subtype_cols.R")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4

myCexAxis <- 1.2
myCexLab <- 1.2

outFolder <- "EXPRVARIANCE_VS_EXPRCORR"
dir.create(outFolder, recursive=TRUE)

coexprFile <- file.path("COEXPR_DE_ANALYSIS2", "tad_coexpr_fc_DT.Rdata")
stopifnot(file.exists(coexprFile))
coexprDT <- eval(parse(text = load(coexprFile)))


varFile <- file.path("EXPR_VARIANCE_BYTAD/LOG2FPKM/all_ds_geneVarDT.Rdata")
stopifnot(file.exists(varFile))


toKeep <- c("hicds", "exprds", "tadMeanVar", "tadMeanVar_cond1", "tadMeanVar_cond2")
varData <- eval(parse(text = load(varFile)))
varData2 <- lapply(varData, function(x) x[toKeep])
varDT <- do.call(rbind, lapply(varData2, data.frame)) # otherwise the columns remain as lists !

varDT3 <- do.call(rbind, lapply(seq_len(length(varData2)), function(i) {
  x <- varData2[[i]]
  rL <- names(x[["tadMeanVar"]])
  tmp2 <- data.frame(x)
  tmp2$region <- rL
  tmp2
  }
  )) # otherwise the columns remain as lists !

varDT <- varDT3

varDT$dataset <- paste0(varDT$hicds, "/", varDT$exprds)

stopifnot(nrow(coexprDT) == nrow(varDT))


varDT_short <- varDT[,c("dataset", "region", "tadMeanVar", "tadMeanVar_cond1", "tadMeanVar_cond2")]
coexprDT_short <- coexprDT[,c("dataset", "region", "cmpCol",
                              "withinCoexpr", "withinCoexpr_cond1", "withinCoexpr_cond2")]
                              # "betweenNbrCoexpr", "betweenNbrCoexpr_cond1", "betweenNbrCoexpr_cond2")]
rownames(varDT_short) <- NULL
rownames(coexprDT_short) <- NULL
cat("... start merging:\n")
cat(paste0(Sys.time(), "\n"))
all_DT <- left_join(varDT_short, coexprDT_short, by=c("dataset", "region"))
cat(paste0(Sys.time(), "\n"))
# all_DTb <- merge(varDT_short, coexprDT_short, by = "dataset")
# cat(paste0(Sys.time(), "\n"))

stopifnot(nrow(coexprDT) == nrow(varDT))
stopifnot(nrow(coexprDT) == nrow(all_DT))

outFile <- file.path(outFolder, "all_DT.Rdata")
save(all_DT, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

all_DT$tadMeanVar_log10 <- log10(all_DT$tadMeanVar)
all_DT$tadMeanVar_cond1_log10 <- log10(all_DT$tadMeanVar_cond1)
all_DT$tadMeanVar_cond2_log10 <- log10(all_DT$tadMeanVar_cond2)

#################################################################################################################
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
myplot_densplot <- function(xvar, yvar, addCurve=FALSE, savePlot=TRUE, mysubtit="") {
  stopifnot(xvar %in% colnames(all_DT))
  myx <- all_DT[,xvar]
  stopifnot(yvar %in% colnames(all_DT))
  myy <- all_DT[,yvar]
  mycols <- all_DT$cmpCol
  if(savePlot) {
    outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "_densplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  }
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
  mtext(side=3, text=mysubtit)
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  if(addCurve) {
    curve(1*x, lty=2, col="grey", add = TRUE)
  }
  if(savePlot){
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}

myplot_colplot <- function(xvar, yvar, mycols, addCurve = FALSE, savePlot=TRUE, mysubtit="") {
  stopifnot(xvar %in% colnames(all_DT))
  myx <- all_DT[,xvar]
  stopifnot(yvar %in% colnames(all_DT))
  myy <- all_DT[,yvar]
  if(savePlot){
    outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "_colplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  }
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
  mtext(side=3, text=mysubtit)
  addCorr(x=myx,y=myy,legPos="topleft", bty='n')
  if(addCurve) {
    curve(1*x, lty=2, col="grey",add=TRUE)
  }
  addSubtypeLeg(bty="n")
  if(savePlot){
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}

mycols <- all_DT$cmpCol
#################################################################################################################

nDS <- length(unique(all_DT$dataset))
nTADs <- nrow(all_DT)

subTit <- paste0("(nDS = ", nDS, " - nTADs = ", nTADs, ")")

############################
############################ tadMeanVar vs. withinCoexpr
############################

myx <- "withinCoexpr"
myy <- "tadMeanVar"


myplot_densplot(xvar = myx, 
                yvar = myy, 
                addCurve = FALSE, 
                savePlot = !SSHFS,
                mysubtit = subTit)

myplot_colplot(xvar = myx, 
                yvar = myy, 
               mycols = mycols,
                addCurve = FALSE, 
                savePlot = !SSHFS,
               mysubtit = subTit)



############################
############################ tadMeanVar (log10) vs. withinCoexpr
############################

myx <- "withinCoexpr"
myy <- "tadMeanVar_log10"


myplot_densplot(xvar = myx, 
                yvar = myy, 
                addCurve = FALSE, 
                savePlot = !SSHFS,
                mysubtit = subTit)

myplot_colplot(xvar = myx, 
               yvar = myy, 
               mycols = mycols,
               addCurve = FALSE, 
               savePlot = !SSHFS,
               mysubtit = subTit)


############################
############################ varDiff vs. withinDiff
############################

all_DT$withinCoexprDiff <- all_DT$withinCoexpr_cond2 - all_DT$withinCoexpr_cond1
all_DT$varDiff <- all_DT$tadMeanVar_cond2 - all_DT$tadMeanVar_cond1


myx <- "withinCoexprDiff"
myy <- "varDiff"


myplot_densplot(xvar = myx, 
                yvar = myy, 
                addCurve = FALSE, 
                savePlot = !SSHFS,
                mysubtit = subTit)

myplot_colplot(xvar = myx, 
               yvar = myy, 
               mycols = mycols,
               addCurve = FALSE, 
               savePlot = !SSHFS,
               mysubtit = subTit)


# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



