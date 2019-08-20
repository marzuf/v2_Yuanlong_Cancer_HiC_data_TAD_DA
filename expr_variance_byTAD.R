# Rscript expr_variance_byTAD.R log2fpkm

options(save.defaults = list(version = 2))


buildTable <- TRUE

cat("> START: expr_variance_byTAD.R\n")

script_name <- "expr_variance_byTAD.R"

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

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

script0_name <- "0_prepGeneData"

all_setting_files <- list.files(settingFolder, full.names=T, recursive = TRUE)
all_setting_files <- all_setting_files[grepl("^run_settings_.+\\.R$", basename(all_setting_files))]

stopifnot(length(all_setting_files) > 0)


outFold <- file.path(paste0("EXPR_VARIANCE_BYTAD"), toupper(exprType))
dir.create(outFold, recursive = TRUE)

plotType <- "png"
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

outNames <- all_setting_files
outSuffix <- gsub("^run_settings_", "", basename(outNames))
outSuffix <- gsub("\\.R$", "", outSuffix)
outPrefix <- basename(dirname(outNames))
outNames <- paste0(outPrefix, "_", outSuffix)

if(buildTable) {
  all_ds_geneVarDT <- foreach(ds_file = all_setting_files) %dopar% {
    
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
    
    ### VARIANCE BY TAD
    stopifnot( setequal(names(geneList), rownames(curr_exprDT)))
    stopifnot( setequal(names(geneList), rownames(curr_exprDT_cond1)))
    stopifnot( setequal(names(geneList), rownames(curr_exprDT_cond2)))
        
    
    cat("... load regionList\n")
    regionList <- eval(parse(text = load(
      file.path(ds_pipFolder, "0_prepGeneData", "pipeline_regionList.Rdata")
    )))
    
    ### RETRIEVE THE GENE2TAD ASSIGNMENT
    g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
    stopifnot(nrow(g2t_DT) > 0)
    stopifnot(geneList %in% g2t_DT$entrezID)
    reg=regionList[1]
    tadMeanVar <- foreach(reg = regionList, .combine='c') %dopar% {
      reg_genes <- g2t_DT$entrezID[g2t_DT$region==reg]
      stopifnot(reg_genes %in% geneList)
      expr_genes <- names(geneList)[geneList %in% reg_genes]
      stopifnot(expr_genes %in% rownames(curr_exprDT))
      tad_exprDT <- curr_exprDT[expr_genes,]
      geneVar <- apply(tad_exprDT, 1,  var, na.rm=T)
      mean(geneVar)
    }
    names(tadMeanVar) <- regionList 
    
    tadMeanVar_cond1 <- foreach(reg = regionList, .combine='c') %dopar% {
      reg_genes <- g2t_DT$entrezID[g2t_DT$region==reg]
      stopifnot(reg_genes %in% geneList)
      expr_genes <- names(geneList)[geneList %in% reg_genes]
      stopifnot(expr_genes %in% rownames(curr_exprDT_cond1))
      tad_exprDT_cond1 <- curr_exprDT_cond1[expr_genes,]
      geneVar_cond1 <- apply(tad_exprDT_cond1, 1,  var, na.rm=T)
      mean(geneVar_cond1)
    }
    names(tadMeanVar_cond1) <- regionList 
    
    
    tadMeanVar_cond2 <- foreach(reg = regionList, .combine='c') %dopar% {
      reg_genes <- g2t_DT$entrezID[g2t_DT$region==reg]
      stopifnot(reg_genes %in% geneList)
      expr_genes <- names(geneList)[geneList %in% reg_genes]
      stopifnot(expr_genes %in% rownames(curr_exprDT_cond2))
      tad_exprDT_cond2 <- curr_exprDT_cond2[expr_genes,]
      geneVar_cond2 <- apply(tad_exprDT_cond2, 1,  var, na.rm=T)
      mean(geneVar_cond2)
    }
    names(tadMeanVar_cond2) <- regionList 
    
        
    list(
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
      
      region = regionList,
      
      tadMeanVar = tadMeanVar,
      tadMeanVar_cond1 = tadMeanVar_cond1,
      tadMeanVar_cond2 = tadMeanVar_cond2
    )
    
  } # end-foreach iterating over datasets
  names(all_ds_geneVarDT) <- outNames
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  save(all_ds_geneVarDT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else{
  outFile <- file.path(outFold, "all_ds_geneVarDT.Rdata")
  all_ds_geneVarDT <- eval(parse(text = load(outFile)))
}
all_ds_geneVar_with_TAD_data <- all_ds_geneVarDT

all_ds_geneVarDT_noTAD <- lapply(all_ds_geneVarDT, function(x) {
  toKeep <- ! names(x) %in% c("tadMeanVar", "tadMeanVar_cond1", "tadMeanVar_cond2")
  x[toKeep]
})
names(all_ds_geneVarDT_noTAD) <- names(all_ds_geneVarDT)


# all_ds_geneVarDT <- data.frame(do.call(rbind, all_ds_geneVarDT_noTAD))
all_ds_geneVarDT <- do.call(rbind, lapply(all_ds_geneVarDT_noTAD, data.frame)) # otherwise the columns remain as lists !


all_ds_geneVarDT_onlyTAD <- lapply(all_ds_geneVar_with_TAD_data, function(x) {
  toKeep <-  names(x) %in% c("tadMeanVar", "tadMeanVar_cond1", "tadMeanVar_cond2")
  x[toKeep]
})
names(all_ds_geneVarDT_onlyTAD) <- names(all_ds_geneVar_with_TAD_data)


# all_ds_geneVarDT <- data.frame(do.call(rbind, all_ds_geneVarDT_noTAD))
all_ds_geneVar_onlyTAD_DT <- do.call(rbind, lapply(all_ds_geneVarDT_onlyTAD, data.frame)) # otherwise the columns remain as lists !

all_ds_geneVar_onlyTAD_DT$hicds <- gsub("^(.+)_TCGA.+", "\\1", rownames(all_ds_geneVar_onlyTAD_DT))
all_ds_geneVar_onlyTAD_DT$exprds <- gsub("^.+_(TCGA.+)\\.chr.+", "\\1", rownames(all_ds_geneVar_onlyTAD_DT))
all_ds_geneVar_onlyTAD_DT$region <- gsub("^.+_TCGA.+\\.(chr.+$)", "\\1", rownames(all_ds_geneVar_onlyTAD_DT))

# stop("--ok")
# load("GENE_VARIANCE/LOG2FPKM/all_ds_geneVarDT.Rdata")


subTypeDT <- data.frame(exprds=names(cancer_subAnnot), subtype=cancer_subAnnot, stringsAsFactors = FALSE)
colorDT <- data.frame(exprds=names(dataset_proc_colors), color=dataset_proc_colors, stringsAsFactors = FALSE)
stopifnot(all_ds_geneVar_onlyTAD_DT$exprds %in% subTypeDT$exprds)
stopifnot(all_ds_geneVar_onlyTAD_DT$exprds %in% colorDT$exprds)

subtype_data_DT <- unique(merge(merge(all_ds_geneVar_onlyTAD_DT, subTypeDT, by="exprds", all.x=T, all.y=F), 
                                colorDT, by="exprds", all.x=T, all.y=F))

stopifnot(!is.na(subtype_data_DT))

stopifnot(nrow(all_ds_geneVarDT) == nrow(subtype_data_DT))

dsByType <- table(unique(subtype_data_DT[, c("exprds", "hicds","subtype")])$subtype)
nDSbyType <- setNames(as.numeric(dsByType), names(dsByType))
subTit <- paste0("# DS: ", paste0(names(nDSbyType), "=", as.numeric(nDSbyType), collapse = " - "))


##########################################################################################
##########################################################################################

subtypeDT <- subtype_data_DT[subtype_data_DT$subtype == "subtypes",]
stopifnot(nrow(subtypeDT) > 0)

nds <- length(unique(paste0(subtypeDT$exprds, subtypeDT$hicds)))
myTit <- "TAD variance by condition"
mySub <- paste0("(cancer_subtypes datasets only; n=", nds, ")")
outFile <- file.path(outFold, paste0("multidens_subtypes_compCond", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens( 
  list(
    tadMeanVar = subtypeDT[,"tadMeanVar"],
    tadMeanVar_cond1 = subtypeDT[,"tadMeanVar_cond1"],
    tadMeanVar_cond2 = subtypeDT[,"tadMeanVar_cond2"]
  ),
  legPos="topleft",
  my_xlab = "TAD mean variance",
  plotTit = myTit
)
mtext(side=3, text=mySub)
foo <- dev.off() 
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("multidens_subtypes_compCond_log10", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens( 
  list(
    tadMeanVar = log10(subtypeDT[,"tadMeanVar"]),
    tadMeanVar_cond1 = log10(subtypeDT[,"tadMeanVar_cond1"]),
    tadMeanVar_cond2 = log10(subtypeDT[,"tadMeanVar_cond2"])
  ),
  legPos="topleft",
  my_xlab = "TAD mean variance",
  plotTit = myTit
)
mtext(side=3, text=mySub)
foo <- dev.off() 
cat(paste0("... written: ", outFile, "\n"))

##########################################################################################
##########################################################################################

mutDT <- subtype_data_DT[subtype_data_DT$subtype == "mutation",]
stopifnot(nrow(mutDT) > 0)

nds <- length(unique(paste0(mutDT$exprds, mutDT$hicds)))
myTit <- "TAD variance by condition"
mySub <- paste0("(wt_vs_mut datasets only; n=", nds, ")")
outFile <- file.path(outFold, paste0("multidens_mutation_compCond", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens( 
  list(
    tadMeanVar = mutDT[,"tadMeanVar"],
    tadMeanVar_cond1 = mutDT[,"tadMeanVar_cond1"],
    tadMeanVar_cond2 = mutDT[,"tadMeanVar_cond2"]
  ),
  legPos="topleft",
  my_xlab = "TAD mean variance",
  plotTit = myTit
)
mtext(side=3, text=mySub)
foo <- dev.off() 
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("multidens_mutation_compCond", "_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens( 
  list(
    tadMeanVar = log10(mutDT[,"tadMeanVar"]),
    tadMeanVar_cond1 = log10(mutDT[,"tadMeanVar_cond1"]),
    tadMeanVar_cond2 = log10(mutDT[,"tadMeanVar_cond2"])
  ),
  legPos="topleft",
  my_xlab = "TAD mean variance",
  plotTit = myTit
)
mtext(side=3, text=mySub)
foo <- dev.off() 
cat(paste0("... written: ", outFile, "\n"))



##########################################################################################
##########################################################################################
normalDT <- subtype_data_DT[subtype_data_DT$subtype == "vs_normal",]
stopifnot(nrow(normalDT) > 0)

nds <- length(unique(paste0(normalDT$exprds, normalDT$hicds)))
myTit <- "TAD variance by condition"
mySub <- paste0("(norm_vs_tumor datasets only; n=", nds, ")")
outFile <- file.path(outFold, paste0("multidens_vs_normal_compCond", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens( 
  list(
    tadMeanVar = normalDT[,"tadMeanVar"],
    tadMeanVar_cond1 = normalDT[,"tadMeanVar_cond1"],
    tadMeanVar_cond2 = normalDT[,"tadMeanVar_cond2"]
  ),
  legPos="topleft",
  my_xlab = "TAD mean variance",
  plotTit = myTit
)
mtext(side=3, text=mySub)
foo <- dev.off() 
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFold, paste0("multidens_vs_normal_compCond", "_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens( 
  list(
    tadMeanVar = log10(normalDT[,"tadMeanVar"]),
    tadMeanVar_cond1 = log10(normalDT[,"tadMeanVar_cond1"]),
    tadMeanVar_cond2 = log10(normalDT[,"tadMeanVar_cond2"])
  ),
  legPos="topleft",
  my_xlab = "TAD mean variance",
  plotTit = myTit
)
mtext(side=3, text=mySub)
foo <- dev.off() 
cat(paste0("... written: ", outFile, "\n"))

##########################################################################################
##########################################################################################
for(mytransf in c("normal", "log10")){
  
  mysuffix <- ifelse(mytransf=="normal", "", "_log10")
  mysuffixLab <- ifelse(mytransf=="normal", "", " [log10] ")
  for(myvar in c("tadMeanVar", "tadMeanVar_cond1", "tadMeanVar_cond2")){
    
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
      legPos="topleft",
      plotTit = myTit
    )
    mtext(side=3, text=subTit)
    foo <- dev.off() 
    cat(paste0("... written: ", outFile, "\n"))
  } 
}

##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
for(mytransf in c("normal", "log10")){
  
  mysuffix <- ifelse(mytransf=="normal", "", "_log10")
  mysuffixLab <- ifelse(mytransf=="normal", "", " [log10] ")
  for(myvar in c("tadMeanVar", "tadMeanVar_cond1", "tadMeanVar_cond2")){
    
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
      legPos="topleft",
      plotTit = myTit
    )
    mtext(side=3, text=subTit)
    foo <- dev.off() 
    cat(paste0("... written: ", outFile, "\n"))
  } 
}

##########################################################################################
##########################################################################################



################
cat("*** DONE - ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

stop("---ok\n")
## REMAINING NOT MEANINGFUL AT THE TAD LEVEL



################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
