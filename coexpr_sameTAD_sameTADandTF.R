startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)
require(dplyr)


# Rscript coexpr_sameTAD_sameTADandTF.R crisp 
# Rscript coexpr_sameTAD_sameTADandTF.R c3.mir
# Rscript coexpr_sameTAD_sameTADandTF.R c3.tft
# Rscript coexpr_sameTAD_sameTADandTF.R c3.all
# Rscript coexpr_sameTAD_sameTADandTF.R trrust
# Rscript coexpr_sameTAD_sameTADandTF.R tftg
# Rscript coexpr_sameTAD_sameTADandTF.R motifmap LG1_40kb
# Rscript coexpr_sameTAD_sameTADandTF.R motifmap 
# Rscript coexpr_sameTAD_sameTADandTF.R kegg LG1_40kb

plotType <- "png"
myHeight <- 400
myWidth <- 400
plotCex <- 1.4

corMet <- "pearson"

dsIn <- "crisp"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1 | length(args) == 3)
dsIn <- args[1]
if(length(args) == 3) {
  all_hicds <- args[2]
  all_exprds <- args[3]
} else {
  all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
}

stopifnot(dsIn %in% c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg"))

inFolder <- paste0("CREATE_SAMETF_", toupper(dsIn))

outFolder <- file.path(paste0("COEXPR_SAMETAD_SAMETADANDTF", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)

corMet <- "pearson"

buildData <- TRUE

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

if(length(args) == 1) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}


sameTF_dt <- get(load(file.path(inFolder, "sameTF_dt.Rdata")))
sameTF_dt$gene1 <- as.character(sameTF_dt$gene1)
sameTF_dt$gene2 <- as.character(sameTF_dt$gene2)
stopifnot(sameTF_dt$gene1 < sameTF_dt$gene2)

foo <- foreach(hicds = all_hicds) %dopar% {
  
  for(exprds in all_exprds[[paste0(hicds)]])  {
 
    coexprDT <- get(load(file.path("CREATE_COEXPR_SORTNODUP", hicds, exprds, corMet, "coexprDT.Rdata")))
    coexprDT$gene1 <- as.character(coexprDT$gene1)
    coexprDT$gene2 <- as.character(coexprDT$gene2)
    
    sameTAD_dt <- get(load(file.path("CREATE_SAME_TAD_SORTNODUP", hicds, "all_TAD_pairs.Rdata")))
    sameTAD_dt$gene1 <- as.character(sameTAD_dt$gene1)
    sameTAD_dt$gene2 <- as.character(sameTAD_dt$gene2)
    
    stopifnot(coexprDT$gene1 < coexprDT$gene2)
    
    stopifnot(sameTAD_dt$gene1 < sameTAD_dt$gene2)
    
    
    coexpr_sameTF_DT <- left_join(coexprDT, sameTF_dt, all.x=TRUE, by=c("gene1", "gene2") )
    coexpr_sameTAD_sameTF_DT <- left_join(coexpr_sameTF_DT, sameTAD_dt, all.x=TRUE, by=c("gene1", "gene2") )
    
    coexpr_sameTAD_sameTF_DT$sameTF[is.na(coexpr_sameTAD_sameTF_DT$sameTF)] <- 0
    coexpr_sameTAD_sameTF_DT$region[is.na(coexpr_sameTAD_sameTF_DT$region)] <- 0
    coexpr_sameTAD_sameTF_DT$region[! coexpr_sameTAD_sameTF_DT$region == 0] <- 1
    
    colnames(coexpr_sameTAD_sameTF_DT)[colnames(coexpr_sameTAD_sameTF_DT) == "region"] <- "sameTAD"
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_sameTAD.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    par(bty="l")
    boxplot(coexpr ~ sameTAD, data = coexpr_sameTAD_sameTF_DT,
            main = paste0(hicds, " - ", exprds))
    mtext(side=3, text = paste0("# 0 = ", sum(coexpr_sameTAD_sameTF_DT$sameTAD==0), "; # 1 = ", sum(coexpr_sameTAD_sameTF_DT$sameTAD==1)))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_sameTF.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    par(bty="l")
    boxplot(coexpr ~ sameTF, data = coexpr_sameTAD_sameTF_DT,
            main = paste0(hicds, " - ", exprds))
    mtext(side=3, text = paste0("# 0 = ", sum(coexpr_sameTAD_sameTF_DT$sameTF==0), "; # 1 = ", sum(coexpr_sameTAD_sameTF_DT$sameTF==1)))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    sameTAD_dt <- coexpr_sameTAD_sameTF_DT[paste0(coexpr_sameTAD_sameTF_DT$sameTAD) == 1,]
    
    sameTAD_dt$TAD_TF <- ifelse( paste0(sameTAD_dt$sameTF) == 1, "sameTAD_sameTF", "sameTAD" )
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_sameTAD_sameTF.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    par(bty="l")
    boxplot(coexpr ~ TAD_TF, data = sameTAD_dt,
            main = paste0(hicds, " - ", exprds))
    mtext(side=3, text = paste0("# sameTAD = ", sum(sameTAD_dt$TAD_TF=="sameTAD"), "; # sameTAD_sameTF = ", sum(sameTAD_dt$TAD_TF=="sameTAD_sameTF")))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
            
    
  }
}
    
