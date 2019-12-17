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
# Rscript coexpr_sameTAD_sameTADandTF.R motifmap
# Rscript coexpr_sameTAD_sameTADandTF.R kegg

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

outFolder <- file.path(paste0("COEXPR_SAMETAD_SAMETADANDTF_", toupper(dsIn)))
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

all_hicds <- "LG1_40kb"
all_exprds <- setNames("TCGAlusc_norm_lusc", "LG1_40kb")


if(buildData) {
  
  all_mean_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    
    ds_mean_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind')  %do% {
      
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
      
      stopifnot(coexpr_sameTAD_sameTF_DT$sameTF == 1 | coexpr_sameTAD_sameTF_DT$sameTF == 0)
      stopifnot(coexpr_sameTAD_sameTF_DT$sameTAD == 1 | coexpr_sameTAD_sameTF_DT$sameTAD == 0)
      
      coexpr_sameTAD_sameTF_DT$sameTAD_lab <- ifelse(coexpr_sameTAD_sameTF_DT$sameTAD == 0, "diffTAD", 
                                                    ifelse(coexpr_sameTAD_sameTF_DT$sameTAD == 1, "sameTAD", NA))
      
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_sameTAD.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l")
      boxplot(coexpr ~ sameTAD_lab, data = coexpr_sameTAD_sameTF_DT,
              xlab="",
              ylab="pairwise coexpr.",
              main = paste0(hicds, " - ", exprds, " (", dsIn, ")"))
      mtext(side=3, text = paste0("# diffTAD = ", sum(coexpr_sameTAD_sameTF_DT$sameTAD==0), "; # sameTAD = ", sum(coexpr_sameTAD_sameTF_DT$sameTAD==1)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      meanSameTAD <- mean(coexpr_sameTAD_sameTF_DT$coexpr[coexpr_sameTAD_sameTF_DT$sameTAD==1])
      meanDiffTAD <- mean(coexpr_sameTAD_sameTF_DT$coexpr[coexpr_sameTAD_sameTF_DT$sameTAD==0])
      
      coexpr_sameTAD_sameTF_DT$sameTF_lab <- ifelse(coexpr_sameTAD_sameTF_DT$sameTF == 0, "diffTF", 
                                                    ifelse(coexpr_sameTAD_sameTF_DT$sameTF == 1, "sameTF", NA))
      
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_sameTF.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l")
      boxplot(coexpr ~ sameTF_lab, data = coexpr_sameTAD_sameTF_DT,
              xlab="",
              ylab="pairwise coexpr.",
              main = paste0(hicds, " - ", exprds))
      mtext(side=3, text = paste0("# diffTF = ", sum(coexpr_sameTAD_sameTF_DT$sameTF==0), "; # sameTF = ", sum(coexpr_sameTAD_sameTF_DT$sameTF==1)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      sameTAD_dt <- coexpr_sameTAD_sameTF_DT[coexpr_sameTAD_sameTF_DT$sameTAD == 1,]
      sameTAD_dt$TAD_TF <- ifelse( paste0(sameTAD_dt$sameTF) == 1, "sameTAD_sameTF", "sameTAD" )
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_sameTAD_sameTF.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l")
      boxplot(coexpr ~ TAD_TF, data = sameTAD_dt,
              xlab="",
              ylab="pairwise coexpr.",
              main = paste0(hicds, " - ", exprds, " (", dsIn, ")"))
      mtext(side=3, text = paste0("# sameTAD = ", sum(sameTAD_dt$TAD_TF=="sameTAD"), "; # sameTAD_sameTF = ", sum(sameTAD_dt$TAD_TF=="sameTAD_sameTF")))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      meanSameTADNoTF <- mean(sameTAD_dt$coexpr[sameTAD_dt$TAD_TF=="sameTAD"])
      meanSameTADSameTF <- mean(sameTAD_dt$coexpr[sameTAD_dt$TAD_TF=="sameTAD_sameTF"])
      
      tmp1_dt <- coexpr_sameTAD_sameTF_DT
      tmp1_dt$sameTAD[tmp1_dt$sameTAD == 0] <- "diffTAD" 
      tmp1_dt$sameTAD[tmp1_dt$sameTAD == 1] <- "sameTAD" 
      colnames(tmp1_dt)[colnames(tmp1_dt) == "sameTAD"] <- "pairType"
      
      tmp2_dt <- sameTAD_dt
      tmp2_dt$TAD_TF[tmp2_dt$TAD_TF == "sameTAD"] <- "sameTAD_diffTF" 
      
      colnames(tmp2_dt)[colnames(tmp2_dt) == "TAD_TF"] <- "pairType"
      
      full_dt <- rbind(tmp1_dt[,c("gene1", "gene2", "coexpr", "pairType")],
                       tmp2_dt[,c("gene1", "gene2", "coexpr", "pairType")])
      
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_allTypes.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l")
      boxplot(coexpr ~ pairType, data = full_dt,
              xlab="",
              ylab="pairwise coexpr.",
              main = paste0(hicds, " - ", exprds, " (", dsIn, ")"))
      mtext(side=3, text = paste0("# sameTAD = ", sum(full_dt$pairType=="sameTAD"), "; # diffTAD = ", sum(full_dt$pairType=="diffTAD"),
                                  "\n# sameTADdiffTF = ", sum(full_dt$pairType=="sameTAD_diffTF"), "; # sameTAD_sameTF = ", sum(full_dt$pairType=="sameTAD_sameTF")))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))



outdt <- data.frame(
  hicds = hicds,
  exprds=exprds,
  meanSameTAD = meanSameTAD,
  meanDiffTAD = meanDiffTAD,
  meanSameTADNoTF = meanSameTADNoTF,
  meanSameTADSameTF = meanSameTADSameTF,
  stringsAsFactors = FALSE
)

save(outdt, file = file.path(outFolder, paste0(hicds, "_", exprds, "_", "outdt.Rdata")), version=2)

outdt
    }
    ds_mean_dt
  }
  
  if(length(all_hicds) == 1 | length(all_exprds) == 1)
    stop("--ok")
  
  outFile <- file.path(outFolder, "all_mean_dt.Rdata")
  save(all_mean_dt, file =  outFile, version=2)
  cat(paste0("... written:" , outFile, "\n"))
  
  
} else {
  outFile <- file.path(outFolder, "all_mean_dt.Rdata")
  all_mean_dt <- get(load(outFile))
}


all_mean_dt$dataset <- file.path(all_mean_dt$hicds, all_mean_dt$exprds)
rownames(all_mean_dt) <- all_mean_dt$dataset
all_mean_dt <- all_mean_dt[,c("meanDiffTAD", "meanSameTAD", "meanSameTADNoTF", "meanSameTADSameTF")]

outFile <- file.path(outFolder, paste0("allDS_allTypes.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="l")
boxplot(all_mean_dt, 
        xlab="",
        ylab = "mean pairwise coexpr.",
        main = paste0("all DS (n = ", nrow(all_mean_dt), ")")
)
mtext(side=3, text = paste0(""))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))





