# Rscript expr_variance.R log2fpkm


options(save.defaults = list(version = 2))


buildTable <- TRUE

cat("> START: expr_variance.R\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")


library(foreach)
library(doMC)
library(tools)
# library(DESeq2)

startTime <- Sys.time()

pointPch <- 16
pointCex <- 1
cexAxis <- 1.2
cexLab <- 1.2

rangeOffset <- 0.15
source(file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/ezh2_utils_fct.R"))
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R"))# cancer_subColors from here
source( file.path("../Cancer_HiC_data_TAD_DA/colors_utils.R"))
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

require(ggplot2)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R")) # cancer_subColors from here
source( file.path("../Cancer_HiC_data_TAD_DA/colors_utils.R")) # !! redefine cancer_subAnnot !!!
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Cancer_HiC_data_TAD_DA/utils_plot_fcts.R")
source("../Cancer_HiC_data_TAD_DA/plot_lolliTAD_funct.R")

cancerDS <- score_DT$dataset[score_DT$process_short == "cancer"]

registerDoMC(ifelse(SSHFS,2,40))

famType1 <- "hgnc"
famType2 <- "hgnc_family_short"

script17_name <- "170revision2EZH2_score_auc_pval_permGenes"

correctTCGA_factor <- 10^6
# I used scaled estimates -> similar to FPKM, but not multiplied by 10^6

args <- commandArgs(trailingOnly = TRUE)
# stopifnot(length(args) > 0)
exprType <- "log2fpkm"
nTopLast <- 1000
if(!is.na(args[1]))
  exprType <- args[1]
if(!is.na(args[2]))
  nTopLast <- as.numeric(args[2])
stopifnot(is.numeric(nTopLast))
stopifnot(exprType %in% c("fpkm", "log2fpkm", "NBvar", "voom"))

exprTypeName <- ifelse(exprType == "fpkm", "FPKM", 
				  ifelse(exprType == "log2fpkm", "log2 FPKM",
					ifelse(exprType == "NBvar", "negative binomial var.",
						ifelse(exprType == "voom", "voom transf.", NA))))
stopifnot(!is.na(exprTypeName))

cat("... START with exprType =\t", exprType, "\n")
cat("... START with nTopLast =\t", nTopLast, "\n")

settingFolder <- file.path("PIPELINE", "INPUT_FILES")
# PIPELINE/INPUT_FILES/NCI-H460_40kb/run_settings_TCGAluad_mutKRAS_mutEGFR.R

dsFold <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_setting_files <- list.files(settingFolder, full.names=T, recursive = TRUE)
all_setting_files <- all_setting_files[grepl("^run_settings_.+\\.R$", basename(all_setting_files))]

stopifnot(length(all_setting_files) > 0)


#auc_coexprdist_fold <- file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP") #>>>>>>>>>>>>>><<<<< to uncomment when will be done !!! 20.08.19
#stopifnot(dir.exists(auc_coexprdist_fold))

outFold <- file.path(paste0("EXPR_VARIANCE"), toupper(exprType))
dir.create(outFold, recursive = TRUE)

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 400, 7)

logFile <- file.path(outFold, "gene_variance_TCGA_logFile.txt")
if(SSHFS) logFile <- ""
if(!SSHFS) file.remove(logFile)

txt <- paste0("... exprType\t=\t", exprType, "\n")
printAndLog(txt, logFile)
txt <- paste0("... nTopLast\t=\t", nTopLast, "\n")
printAndLog(txt, logFile)
txt <- paste0("... correctTCGA_factor\t=\t", correctTCGA_factor, "\n")
printAndLog(txt, logFile)

ds_file = all_setting_files[1]
# all_setting_files <- all_setting_files[22:22]

# all_setting_files <- all_setting_files[1:3]

if(buildTable) {
  all_ds_geneVarDT <- foreach(ds_file = all_setting_files, .combine="rbind") %dopar% {
    
    hicds <- basename(dirname(ds_file))
    exprds <- gsub("run_settings_(.+)\\.R", "\\1", basename(ds_file))
    
    curr_ds_name <- paste0(hicds, "_", exprds)
    curr_ds_path <- file.path(hicds, exprds)
    ds_pipFolder <- file.path(dsFold, hicds, exprds)
    
    cat(file.path(dsFold, hicds, exprds))
    stopifnot(dir.exists(file.path(dsFold, hicds, exprds)))
      
    cat("... source settingFile", basename(ds_file), "\n")
    source(ds_file)
    cat("... load samp1\n")
    samp1 <- eval(parse(text = load(file.path(setDir, sample1_file))))
    cat("... load samp2\n")
    samp2 <- eval(parse(text = load(file.path(setDir, sample2_file))))
    
    cat("... load geneList\n")
    geneList <- eval(parse(text = load(
      file.path(ds_pipFolder, "0_prepGeneData", "pipeline_geneList.Rdata")
    )))
    
    cat("... load exprDT\n")
    
    if(exprType == "fpkm" | exprType == "log2fpkm" | exprType == "NBvar") {
      exprDT <- eval(parse(text = load(
        file.path(ds_pipFolder, "0_prepGeneData", paste0("rna_", "fpkm", "DT.Rdata"))
      )))
    } else if(exprType == "voom") {
      exprDT <- eval(parse(text = load(
        file.path(ds_pipFolder, "1_runGeneDE", paste0(exprType, "_lmFitInputDT.Rdata"))
      )))
    } else{
      stop(paste0("!!! unimplemented exprType", exprType, "\n"))
    }
    
    stopifnot(samp1 %in% colnames(exprDT))
    stopifnot(samp2 %in% colnames(exprDT))
    stopifnot(names(geneList) %in% rownames(exprDT))
    
    curr_exprDT <- exprDT[rownames(exprDT) %in% names(geneList), c(samp1,samp2)]  
    stopifnot(is.numeric(curr_exprDT[1,1]))
    
    if(grepl("TCGA", curr_ds_name)) {
      cat("!!! For TCGA data, correctTCGA_factor = ", correctTCGA_factor, "\n")
      curr_exprDT <- curr_exprDT * correctTCGA_factor
    } else {
      stop("error")
    }

    if(exprType == "log2fpkm") {
      curr_exprDT <- log2(curr_exprDT + 1)
      stopifnot(is.numeric(curr_exprDT[1,1]))
    } else if(exprType == "NBvar") {
        cat("!!! WARNING to use \"DESeq2:::varianceStabilizingTransformation\" the expression data are rounded !!!\n")
        curr_exprDT <- try(DESeq2:::varianceStabilizingTransformation(round(as.matrix(curr_exprDT))))
        if(class(curr_exprDT) == "try-error"){
          return(data.frame(
            nTopLast = nTopLast,
            data_type = exprType,
            hicds = hicds,
            exprds=exprds,
            dataset = curr_ds_name,
            meanMostVar = NA,
            meanLeastVar = NA,
            stringsAsFactors = FALSE
          ))
        }
          
        stopifnot(is.numeric(curr_exprDT[1,1]))
      }
    
    geneVar <- apply(curr_exprDT, 1,  var, na.rm=T)
    geneVar <- sort(geneVar, decreasing = TRUE)
      
    stopifnot(length(geneVar) >= nTopLast)
  
    mostVariant <- geneVar[1:nTopLast]
    leastVariant <- geneVar[(length(geneVar)-nTopLast+1):length(geneVar)]
    stopifnot(length(mostVariant) == nTopLast)
    stopifnot(length(leastVariant) == nTopLast)
    
    meanMostVar <- mean(mostVariant)
    meanLeastVar <- mean(leastVariant)
    
    medianMostVar <- median(mostVariant)
    medianLeastVar <- median(leastVariant)
    
    ################ COND1 ONLY
    stopifnot(samp1 %in% colnames(curr_exprDT))
    
    curr_exprDT_cond1 <- curr_exprDT[, samp1]
    stopifnot(dim(curr_exprDT_cond1) == c(nrow(curr_exprDT), length(samp1)))
    
    geneVar_cond1 <- apply(curr_exprDT_cond1, 1,  var, na.rm=T)
    geneVar_cond1 <- sort(geneVar_cond1, decreasing = TRUE)
    
    stopifnot(length(geneVar_cond1) >= nTopLast)
    
    mostVariant_cond1 <- geneVar_cond1[1:nTopLast]
    leastVariant_cond1 <- geneVar_cond1[(length(geneVar_cond1)-nTopLast+1):length(geneVar_cond1)]
    stopifnot(length(mostVariant_cond1) == nTopLast)
    stopifnot(length(leastVariant_cond1) == nTopLast)
    
    meanMostVar_cond1 <- mean(mostVariant_cond1)
    meanLeastVar_cond1 <- mean(leastVariant_cond1)
    
    medianMostVar_cond1 <- median(mostVariant_cond1)
    medianLeastVar_cond1 <- median(leastVariant_cond1)
    
    ################ COND2 ONLY
    curr_exprDT_cond2 <- curr_exprDT[, samp2]
    stopifnot(dim(curr_exprDT_cond2) == c(nrow(curr_exprDT), length(samp2)))
    
    geneVar_cond2 <- apply(curr_exprDT_cond2, 1,  var, na.rm=T)
    geneVar_cond2 <- sort(geneVar_cond2, decreasing = TRUE)
    
    stopifnot(length(geneVar_cond2) >= nTopLast)
    
    mostVariant_cond2 <- geneVar_cond2[1:nTopLast]
    leastVariant_cond2 <- geneVar_cond2[(length(geneVar_cond2)-nTopLast+1):length(geneVar_cond2)]
    stopifnot(length(mostVariant_cond2) == nTopLast)
    stopifnot(length(leastVariant_cond2) == nTopLast)
    
    meanMostVar_cond2 <- mean(mostVariant_cond2)
    meanLeastVar_cond2 <- mean(leastVariant_cond2)
    
    medianMostVar_cond2 <- median(mostVariant_cond2)
    medianLeastVar_cond2 <- median(leastVariant_cond2)
    
    
    data.frame(
      nTopLast = nTopLast,
      data_type = exprType,
      hicds = hicds,
      exprds=exprds,
      dataset = curr_ds_name,
      meanMostVar = meanMostVar,
      meanLeastVar = meanLeastVar,
      medianMostVar = medianMostVar,
      medianLeastVar = medianLeastVar,
      
      meanMostVar_cond1 = meanMostVar_cond1,
      meanLeastVar_cond1 = meanLeastVar_cond1,
      medianMostVar_cond1 = medianMostVar_cond1,
      medianLeastVar_cond1 = medianLeastVar_cond1,
      
      meanMostVar_cond2 = meanMostVar_cond2,
      meanLeastVar_cond2 = meanLeastVar_cond2,
      medianMostVar_cond2 = medianMostVar_cond2,
      medianLeastVar_cond2 = medianLeastVar_cond2,
      
      
      stringsAsFactors = FALSE
    )
    
  }
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  save(all_ds_geneVarDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else{
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  all_ds_geneVarDT <- eval(parse(text = load(outFile)))
}

stop("--ok")
# load("GENE_VARIANCE/LOG2FPKM/all_ds_geneVarDT.Rdata")

########################################################################################## RETRIEVE FCC AND COEXPRDIST
##########################################################################################

all_ds <- unique(all_ds_geneVarDT$dataset)

curr_ds <- all_ds[1]

### RETRIEVE FCC
aucFCC <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  
  hicds <- all_ds_geneVarDT$hicds[all_ds_geneVarDT$dataset == curr_ds]
  stopifnot(length(hicds) == 1)
  exprds <- all_ds_geneVarDT$exprds[all_ds_geneVarDT$dataset == curr_ds]
  stopifnot(length(exprds) == 1)
  
  step17_fold <- file.path(dsFold, hicds, exprds, script17_name)
  stopifnot(dir.exists(step17_fold))
  
  aucFCC_file <- file.path(step17_fold, "auc_ratios.Rdata")
  stopifnot(file.exists(aucFCC_file))
  
  # AUC_COEXPRDIST_WITHFAM_SORTNODUP/MCF-7_40kb/TCGAbrca_lum_bas_hgnc/hgnc_family_short/auc_values.Rdata
  aucCoexprDist_file <- file.path(auc_coexprdist_fold, hicds, paste0(exprds, "_", famType1), famType2, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  
  all_ratios <- eval(parse(text = load(aucFCC_file)))
  aucFCC <- as.numeric(all_ratios["prodSignedRatio_auc_permGenes"])
  stopifnot(!is.na(aucFCC))
  aucFCC
}
names(aucFCC) <- all_ds
outFile <- file.path(outFold, "aucFCC.Rdata")
save(aucFCC, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# load("GENE_VARIANCE/LOG2FPKM/aucFCC.Rdata")
name_aucFCC_expr <- gsub(".+(TCGA.+)$", "\\1", names(aucFCC))

stopifnot(name_aucFCC_expr %in% names(dataset_proc_colors) )
stopifnot(cancer_subAnnot %in% names(cancer_subColors))

curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[name_aucFCC_expr])])
stopifnot(!is.na(curr_colors))

aucCoexprDist <- foreach(curr_ds = all_ds, .combine='c') %dopar% {
  
  hicds <- all_ds_geneVarDT$hicds[all_ds_geneVarDT$dataset == curr_ds]
  stopifnot(length(hicds) == 1)
  exprds <- all_ds_geneVarDT$exprds[all_ds_geneVarDT$dataset == curr_ds]
  stopifnot(length(exprds) == 1)
  
  ### RETRIEVE FCC
  step17_fold <- file.path(dsFold, hicds, exprds, script17_name)
  
  aucCoexprDist_file <- file.path(auc_coexprdist_fold, hicds, paste0(exprds, "_", famType1), famType2, "auc_values.Rdata")
  stopifnot(file.exists(aucCoexprDist_file))
  
  all_aucDist <- eval(parse(text = load(aucCoexprDist_file)))
  aucCoexprDist <- as.numeric(all_aucDist["auc_ratio_same_over_diff_distVect"])
  stopifnot(!is.na(aucCoexprDist))
  aucCoexprDist
}
names(aucCoexprDist) <- all_ds
outFile <- file.path(outFold, "aucCoexprDist.Rdata")
save(aucCoexprDist, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

# load("GENE_VARIANCE/LOG2FPKM/all_ds_geneVarDT.Rdata")
# load("GENE_VARIANCE/LOG2FPKM/aucFCC.Rdata")
# load("GENE_VARIANCE/LOG2FPKM/aucCoexprDist.Rdata")

subTypeDT <- data.frame(exprds=names(cancer_subAnnot), subtype=cancer_subAnnot, stringsAsFactors = FALSE)
colorDT <- data.frame(exprds=names(dataset_proc_colors), color=dataset_proc_colors, stringsAsFactors = FALSE)
stopifnot(all_ds_geneVarDT$exprds %in% subTypeDT$exprds)
stopifnot(all_ds_geneVarDT$exprds %in% colorDT$exprds)

subtype_data_DT <- unique(merge(merge(all_ds_geneVarDT, subTypeDT, by="exprds", all.x=T, all.y=F), 
                                colorDT, by="exprds", all.x=T, all.y=F))

stopifnot(!is.na(subtype_data_DT))

stopifnot(nrow(all_ds_geneVarDT) == nrow(subtype_data_DT))

dsByType <- table(unique(subtype_data_DT[, c("exprds", "subtype")])$subtype)
nDSbyType <- setNames(as.numeric(dsByType), names(dsByType))
subTit <- paste0("# DS: ", paste0(names(nDSbyType), "=", as.numeric(nDSbyType), collapse = " - "))

##########################################################################################
##########################################################################################
for(mytransf in c("normal", "log10")){
  
  mysuffix <- ifelse(mytransf=="normal", "", "_log10")
  mysuffixLab <- ifelse(mytransf=="normal", "", " [log10] ")
  for(myvar in c("meanMostVar", "medianMostVar", "meanLeastVar", "medianLeastVar")){
    
    myTit <- ifelse(mytransf == "normal", 
                    paste0(myvar, " by categories"),
                    paste0(myvar, " (log10) by categories"))
    
    outFile <- file.path(outFold, paste0("multidens_bySubtypes_", myvar, mysuffix, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens_setcols(
      lapply(split(subtype_data_DT, subtype_data_DT$subtype), function(x) {
        tmp <- x[[myvar]]
        if(mytransf == "log10") tmp <- log10(tmp)
        tmp
        }),
      my_cols = cancer_subColors[as.character(levels(as.factor(subtype_data_DT$subtype)))],
      my_xlab = paste0(myvar),
      plotTit = myTit
    )
    mtext(side=3, text=subTit)
    foo <- dev.off() 
    cat(paste0("... written: ", outFile, "\n"))
  } 
}

##########################################################################################
##########################################################################################

myTit <- paste0("mean variance top and least ", nTopLast, " genes (", exprType, " data)")

mySubTit <- paste0("all datasets (n=", nrow(all_ds_geneVarDT), ")")

if(exprType == "NBvar") {
  all_ds_geneVarDT <- all_ds_geneVarDT[!is.na(all_ds_geneVarDT$meanMostVar) & !is.na(all_ds_geneVarDT$meanMostVar),]
}

####### scatter plot most var, least var, and ratio

all_auc <- c("FCC", "CoexprDist")

auc_type ="FCC"

for(auc_type in all_auc) {
  
  cat("... start", auc_type, "\n")
  
  curr_auc <- eval(parse(text = paste0("auc", auc_type)))
  curr_auc <- (curr_auc-1)*100
  stopifnot(all_ds_geneVarDT$dataset %in% names(curr_auc))
  
  all_ds_geneVarDT[, auc_type] <- curr_auc[all_ds_geneVarDT$dataset]
  
  myTit <- paste0(auc_type, ": mean variance most ", nTopLast, " genes (", exprType, " data)")
  outFile <- file.path(outFold, paste0(auc_type, "auc_vs_mostVar_log10.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  plot(
       x =  log10(all_ds_geneVarDT$meanMostVar),
       y =  all_ds_geneVarDT[, auc_type],
       xlim =  range(log10(all_ds_geneVarDT$meanMostVar)) + c(-rangeOffset,rangeOffset),
       ylim =  range(all_ds_geneVarDT[, auc_type])  + c(-rangeOffset,rangeOffset),
       pch = 16, cex = 0.7,
       ylab = paste0("% AUC increase - ", auc_type),
       # xlab = "mean most var (log10)",
       xlab = paste0("Mean most var. [log10] top ", nTopLast, " most variant genes\n(", exprTypeName, ")"),
       main = myTit
  )
  text(x = log10(all_ds_geneVarDT$meanMostVar),
       y = all_ds_geneVarDT[, auc_type],
       labels = all_ds_geneVarDT$dataset,
       pos=3, cex = 0.7)
  add_legend_narm(x = log10(all_ds_geneVarDT$meanMostVar), y = all_ds_geneVarDT[, auc_type], mypos="bottomright", mymet="Pearson")  
  
  mtext(text = mySubTit, side = 3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

add_curv_fit <- function(x, y, withR2 = TRUE, R2shiftX = 0, R2shiftY = 0,...) {
  mymodel <- lm(y~x)
  abline(mymodel, ...)
  if(withR2) {
    r2Txt <- paste0("adj. R2 = ", sprintf("%.2f", summary(mymodel)$adj.r.squared))
    r2X <- x[which.min(x)] + R2shiftX
    r2Y <- fitted(mymodel)[which.min(x)]  + R2shiftY
    text(x = r2X, y = r2Y, 
         labels = r2Txt, 
         font=3,
         adj=c(1,0),
         pos=3,
         cex = 0.7)
  }
}

for(auc_type in all_auc) {
  
  cat("... start", auc_type, "\n")
  curr_auc <- eval(parse(text = paste0("auc", auc_type)))
  curr_auc <- (curr_auc-1)*100

  stopifnot(all_ds_geneVarDT$dataset %in% names(curr_auc))
  all_ds_geneVarDT[, auc_type] <- curr_auc[all_ds_geneVarDT$dataset]
  
  myTit <- paste0(auc_type, ": mean variance most ", nTopLast, " genes (", exprType, " data)")
  outFile <- file.path(outFold, paste0(auc_type, "auc_vs_mostVar_log10_withFit_noLab.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  xvar <- log10(all_ds_geneVarDT$meanMostVar)
  yvar <- all_ds_geneVarDT[, auc_type]
  
  plot(y =  yvar,
       x =  xvar,
       xlim = range(xvar) + c(-rangeOffset, rangeOffset),
       ylim = range(yvar) + c(-rangeOffset, rangeOffset),
       pch = pointPch, cex = pointCex,
       ylab = paste0("% AUC increase - ", auc_type),
       #xlab = "mean most var (log10)",
	   xlab = paste0("Mean var. [log10] top ", nTopLast, " most variant genes\n(", exprTypeName, ")"),
       main = myTit,
cex.axis = cexAxis, cex.lab = cexLab
  )
  # text(x = log10(all_ds_geneVarDT$meanMostVar),
  #      y = all_ds_geneVarDT[, auc_type],
  #      labels = all_ds_geneVarDT$dataset,
  #      pos=3, cex = 0.7)
  add_legend_narm(x = xvar, y = yvar, mypos="bottomright", mymet="Pearson")  
  add_curv_fit(x = xvar, y = yvar, 
               R2shiftX = 0, 
               R2shiftY = 0,
               withR2 = TRUE,
               col ="gray80",
               lty=2)
  
  mtext(text = mySubTit, side = 3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  myTit <- paste0(auc_type, ": mean variance most ", nTopLast, " genes (", exprType, " data)")
  outFile <- file.path(outFold, paste0(auc_type, "auc_vs_mostVar_log10_withFit.", plotType))
  do.call(plotType, list(outFile, height = myHeight, width = myWidth))
  
  xvar <- log10(all_ds_geneVarDT$meanMostVar)
  yvar <- all_ds_geneVarDT[, auc_type]
  
  plot(y =  yvar,
       x =  xvar,
       xlim = range(xvar) + c(-rangeOffset, rangeOffset),
       ylim = range(yvar) + c(-rangeOffset, rangeOffset),
       pch = pointPch, cex = pointCex,col = curr_colors,
       ylab = paste0("% AUC increase - ", auc_type),
       #xlab = "mean most var (log10)",
       xlab = paste0("Mean var. [log10] top ", nTopLast, " most variant genes\n(", exprTypeName, ")"),
       main = myTit,
       cex.axis = cexAxis, cex.lab = cexLab
  )
  text(x = log10(all_ds_geneVarDT$meanMostVar),
       y = all_ds_geneVarDT[, auc_type],
       labels = all_ds_geneVarDT$dataset,
       col = curr_colors,
       pos=3, cex = 0.7)
  add_legend_narm(x = xvar, y = yvar, mypos="bottomright", mymet="Pearson")  
  add_curv_fit(x = xvar, y = yvar, 
               withR2 = TRUE,
               col ="gray80",
               R2shiftX = -0.1, 
               R2shiftY = -0.1,
               lty=2)
  mtext(text = mySubTit, side = 3)
  
  my_colors <- my_colors[my_colors %in% curr_colors]
  legend("topleft",
         legend=unique(cancer_subAnnot[names(aucFCC)]),
         lty=1,
         col = unique(curr_colors),
         lwd = 5,
         bty="n",
         cex = 0.7)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

# load("EXPR_VARIANCE/LOG2FPKM/all_ds_geneVarDT.Rdata")


all_vars <- c("meanMostVar", "medianMostVar")
for(var in all_vars) {
  myx <-  all_ds_geneVarDT[, paste0(var, "_cond1")]
  myy <-  all_ds_geneVarDT[, paste0(var, "_cond2")]
  plot(
    x =myx,
    y = myy,
    main = paste0(var),
    cex.axis = cexAxis, cex.lab = cexLab,
    xlab = "cond1",
    ylab= "cond2"
  )
  mtext(side = 3, text = paste0("cond2 vs. cond1"))
  curve(1*x, lty=2, col="grey")
  addCorr(x = myx ,
          y = myy,
          bty="n")
}



################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
