
# Rscript tadSignif_cnv.R LG1_40kb TCGAluad_norm_luad
# Rscript tadSignif_cnv.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc
# Rscript tadSignif_cnv.R  ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad
# Rscript tadSignif_cnv.R <hicds> <exprds>

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(reshape2)
require(foreach)
require(doMC)
registerDoMC(40)


setDir <- "/media/electron"
setDir <- ""
load(file.path(setDir, "/mnt/ndata/marco/databank/TCGA/TCGA_PancanAtlas/cnv/pancan_cnv/data_CNA.RData"))

# > cna[1:5,1:5]
# TCGA-OR-A5J1-01 TCGA-OR-A5J2-01 TCGA-OR-A5J3-01 TCGA-OR-A5J4-01 TCGA-OR-A5J5-01
# ACAP3                 0               0               0               1               0
# ACTRT2                0               0               0               1               0
# 
# Hugo_Symbol Entrez_Gene_Id Cytoband
# ACAP3         ACAP3         116983  1p36.33
# ACTRT2       ACTRT2         140625  1p36.32
# AGRN           AGRN         375790  1p36.33

outFolder <- file.path("TADSIGNIF_CNV")
dir.create(outFolder, recursive = TRUE)


# args <- commandArgs(trailingOnly = TRUE)
# hicds <- args[1]
# exprds <- args[2]

plotType <- "png"
plotCex <- 1.4
myHeight <- myWidth <- 400

buildData <- TRUE

tad_signif_thresh <- 0.01
signif_col <- "red"
not_signif_col <- "grey"

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))




all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))

if(buildData){
  
  all_cnv_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
  # all_cnv_dt <- foreach(hicds = all_hicds[1], .combine='rbind') %dopar%{
    
    
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    # exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]][1], .combine='rbind') %do% {
  
      pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
      
      

      result_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds,]
      
      settingFile <- file.path("PIPELINE/INPUT_FILES/", hicds, paste0("run_settings_", exprds, ".R"))
      source(settingFile)
      samp1 <- get(load(file.path(setDir, sample1_file)))
      samp2 <- get(load(file.path(setDir, sample2_file)))
      
      g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
      g2t_dt <- read.delim(g2t_file, header=FALSE, stringsAsFactors = FALSE,
                           col.names=c("entrezID", "chromo", "start", "end","region"))
      g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
      
      pipeline_geneList <- get(load(file.path(pipFolder, "0_prepGeneData", "pipeline_geneList.Rdata")))
      
      hugo_cnv_pip <- cna_meta$Hugo_Symbol[cna_meta$Entrez_Gene_Id %in% pipeline_geneList]
      stopifnot(length(hugo_cnv_pip) > 0)
      
      stopifnot(hugo_cnv_pip %in% rownames(cna))
      
      tmp_dt <- cna_meta[cna_meta$Hugo_Symbol %in% hugo_cnv_pip,]
      stopifnot(!duplicated(tmp_dt$Hugo_Symbol))
      hugo2entrez <- setNames(tmp_dt$Entrez_Gene_Id, tmp_dt$Hugo_Symbol)
      
      cnv_samp1_dt <- cna[rownames(cna) %in% hugo_cnv_pip , colnames(cna) %in% samp1, drop=FALSE]
      cnv_samp2_dt <- cna[rownames(cna) %in% hugo_cnv_pip , colnames(cna) %in% samp2, drop=FALSE]
      stopifnot(rownames(cnv_samp1_dt) %in% names(hugo2entrez))
      stopifnot(rownames(cnv_samp2_dt) %in% names(hugo2entrez))
      
      cnv_samp1_dt <- as.data.frame(cnv_samp1_dt)
      cnv_samp1_dt$entrez <- hugo2entrez[paste0(rownames(cnv_samp1_dt))]
      
      cnv_samp2_dt <- as.data.frame(cnv_samp2_dt)
      cnv_samp2_dt$entrez <- hugo2entrez[paste0(rownames(cnv_samp2_dt))]
      
      nSamp1 <- ncol(cnv_samp1_dt) - 1
      nSamp2 <- ncol(cnv_samp2_dt) - 1
      
      
      if(ncol(cnv_samp1_dt) > 1) { # only entrez
        stopifnot(cnv_samp1_dt$entrez %in% g2t_dt$entrezID)
        long_dt_samp1 <- melt(cnv_samp1_dt, id="entrez")
        mean_dt_samp1 <- aggregate(value ~ entrez, data=long_dt_samp1, FUN=mean)
        stopifnot(mean_dt_samp1$entrez %in% pipeline_geneList)
        
      } else {
        mean_dt_samp1 <- cnv_samp1_dt
        mean_dt_samp1$value <- 0  
      }
      
      
      if(ncol(cnv_samp2_dt) > 1) {
        
        stopifnot(cnv_samp2_dt$entrez %in% g2t_dt$entrezID)
        long_dt_samp2 <- melt(cnv_samp2_dt, id="entrez")
        mean_dt_samp2 <- aggregate(value ~ entrez, data=long_dt_samp2, FUN=mean)
        stopifnot(mean_dt_samp2$entrez %in% pipeline_geneList)
        
      } else {
        mean_dt_samp2 <- cnv_samp2_dt
        mean_dt_samp2$value <- 0
      }
      rownames(mean_dt_samp1) <- NULL
      rownames(mean_dt_samp2) <- NULL
      
      colnames(mean_dt_samp1)[colnames(mean_dt_samp1) == "value"] <- "mean_samp1"
      colnames(mean_dt_samp2)[colnames(mean_dt_samp2) == "value"] <- "mean_samp2"
      
      mean_dt_samp12 <- merge(mean_dt_samp1, mean_dt_samp2, by="entrez", all=TRUE)
      stopifnot(!is.na(mean_dt_samp12))
      
      de_dt <- get(load(file.path(pipFolder, "1_runGeneDE", "DE_topTable.Rdata")))
      de_dt <- de_dt[de_dt$genes %in% names(pipeline_geneList),]
      de_dt$entrez <- pipeline_geneList[paste0(de_dt$genes)]
      stopifnot(!is.na(de_dt$entrez))
      
      mean_dt_samp12$entrez <- as.character(mean_dt_samp12$entrez)
      de_dt$entrez <- as.character(de_dt$entrez)
      
      stopifnot(mean_dt_samp12$entrez %in% de_dt$entrez)
      
      # mean_dt_samp12 <- mean_dt_samp12[mean_dt_samp12$entrez %in% de_dt$entrez,]
      
      mean_dt_samp12_logFC <- merge(mean_dt_samp12, de_dt[,c("adj.P.Val", "entrez", "logFC")], by="entrez", all.x=TRUE, all.y=FALSE)
      stopifnot(!is.na(mean_dt_samp12_logFC))
      
      mean_dt_samp12_logFC$diff_cnv_samp12 <- mean_dt_samp12_logFC$mean_samp1-mean_dt_samp12_logFC$mean_samp2
      
      subTit <- paste0("#samp1=", nSamp1, "/", length(samp1), "; #samp2=", nSamp2, "/", length(samp2))
      plotTit <- paste0(hicds, " - ", exprds)
      plotTit <- paste0(hicds, "\n", exprds)
      
      my_x <- mean_dt_samp12_logFC$diff_cnv_samp12
      my_y <- mean_dt_samp12_logFC$logFC
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds,"_diffCNV_geneLogFC.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l")
      densplot(
        x=my_x,
        y=my_y,
        cex=0.7,
        cex.lab=plotCex,
        cex.main=plotCex,
        cex.axis=plotCex,
        xlab = "Mean gene diff. CNV samp1-samp2",
        ylab = "Gene logFC",
        sub=subTit,
        main = plotTit
      )
      # mtext(side = 3, text = subTit)
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
      
      mtext(side = 3, text = paste0("n=", nrow(mean_dt_samp12_logFC)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      my_x <- mean_dt_samp12_logFC$diff_cnv_samp12
      my_y <- mean_dt_samp12_logFC$adj.P.Val
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds,"_diffCNV_geneAdjPval.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l")
      densplot(
        x=my_x,
        y=my_y,
        cex=0.7,
        cex.lab=plotCex,
        cex.main=plotCex,
        cex.axis=plotCex,
        xlab = "Mean gene diff. CNV samp1-samp2",
        ylab = "Gene adj.P.Val",
        sub=subTit,
        main = plotTit
      )
      # mtext(side = 3, text = subTit)
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
      
      mtext(side = 3, text = paste0("n=", nrow(mean_dt_samp12_logFC)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      g2t_dt$entrez <- as.character(g2t_dt$entrezID)
      
      mean_cnv2t_dt <- merge(mean_dt_samp12, g2t_dt[,c("entrez", "region")], by="entrez", all.x=TRUE, all.y=FALSE)
      
      save(g2t_dt, version=2, file="g2t_dt.Rdata")
      save(mean_dt_samp12, version=2, file="mean_dt_samp12.Rdata")
      
      
      stopifnot(!is.na(mean_cnv2t_dt))
      stopifnot(grepl("_TAD", mean_cnv2t_dt$region))
        
      meanCNV_samp1_reg_dt <- aggregate(mean_samp1 ~ region, data = mean_cnv2t_dt, FUN=mean)
      meanCNV_samp2_reg_dt <- aggregate(mean_samp2 ~ region, data = mean_cnv2t_dt, FUN=mean)
      meanCNV_samp12_reg_dt <- merge(meanCNV_samp1_reg_dt, meanCNV_samp2_reg_dt, by="region", all=TRUE)                           
      stopifnot(!is.na(meanCNV_samp12_reg_dt))
      
      stopifnot(!duplicated(result_dt$region))
      
      meanCNV_reg_signif_dt <- merge(meanCNV_samp12_reg_dt, result_dt[,c("region", "meanLogFC", "adjPvalComb")], by="region", all.x=TRUE, all.y=FALSE)
      stopifnot(!is.na(meanCNV_reg_signif_dt))                                  
      meanCNV_reg_signif_dt$diff_cnv_samp12 <- meanCNV_reg_signif_dt$mean_samp1 - meanCNV_reg_signif_dt$mean_samp2
      
      meanCNV_reg_signif_dt$dotCol <- ifelse(meanCNV_reg_signif_dt$adjPvalComb <= tad_signif_thresh, signif_col, not_signif_col )
      
      my_x <- meanCNV_reg_signif_dt$diff_cnv_samp12
      my_y <- meanCNV_reg_signif_dt$meanLogFC
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds,"_diffCNV_tadLogFC.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l")
      densplot(
        x=my_x,
        y=my_y,
        cex=0.7,
        cex.lab=plotCex,
        cex.main=plotCex,
        cex.axis=plotCex,
        xlab = "Mean TAD diff. CNV samp1-samp2",
        ylab = "TAD meanLogFC",
        sub=subTit,
        main = plotTit
      )
      # mtext(side = 3, text = subTit)
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
      mtext(side = 3, text = paste0("n=", nrow(meanCNV_reg_signif_dt)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      my_x <- meanCNV_reg_signif_dt$diff_cnv_samp12
      my_y <- meanCNV_reg_signif_dt$meanLogFC
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds,"_diffCNV_tadLogFC_signifCol.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l")
      plot(
        x=my_x,
        y=my_y,
        col = meanCNV_reg_signif_dt$dotCol,
        cex=0.7,
        pch = 16,
        cex.lab=plotCex,
        cex.main=plotCex,
        cex.axis=plotCex,
        xlab = "Mean TAD diff. CNV samp1-samp2",
        ylab = "TAD meanLogFC",
        sub=subTit,
        main = plotTit
      )
      points(
        x=my_x[meanCNV_reg_signif_dt$dotCol == signif_col],
        y=my_y[meanCNV_reg_signif_dt$dotCol == signif_col],
        col = meanCNV_reg_signif_dt$dotCol[meanCNV_reg_signif_dt$dotCol == signif_col],
        cex=1.2,
        pch = 16
      )
      legend("bottomleft",
             legend=c(
               paste0("signif. TADs (p-val<=", tad_signif_thresh, ")\nn=",sum(meanCNV_reg_signif_dt$dotCol == signif_col))
             ),
             pch = 16,
             # col = c(signif_col, not_signif_col),
             col = c(signif_col),
             bty="n"
             )
      # mtext(side = 3, text = subTit)
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
      mtext(side = 3, text = paste0("n=", nrow(meanCNV_reg_signif_dt)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      
      my_x <- meanCNV_reg_signif_dt$diff_cnv_samp12
      my_y <- meanCNV_reg_signif_dt$adjPvalComb
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds,"_diffCNV_tadAdjCombPval.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      
      par(bty="l")
      densplot(
        x=my_x,
        y=my_y,
        cex=0.7,
        cex.lab=plotCex,
        cex.main=plotCex,
        cex.axis=plotCex,
        xlab = "Mean TAD diff. CNV samp1-samp2",
        ylab = "TAD adj. comb. p-val.",
        sub=subTit,
        main = plotTit
      )
      # mtext(side = 3, text = subTit)
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
      
      mtext(side = 3, text = paste0("n=", nrow(meanCNV_reg_signif_dt)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds,"_diffCNV_tadAdjCombPval_signifCol.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))

      my_x <- meanCNV_reg_signif_dt$diff_cnv_samp12
      my_y <- meanCNV_reg_signif_dt$adjPvalComb
            
      par(bty="l")
      plot(
        x=my_x,
        y=my_y,
        col=meanCNV_reg_signif_dt$dotCol,
        cex=0.7,
        pch = 16,
        cex.lab=plotCex,
        cex.main=plotCex,
        cex.axis=plotCex,
        xlab = "Mean TAD diff. CNV samp1-samp2",
        ylab = "TAD adj. comb. p-val.",
        sub=subTit,
        main = plotTit
      )
      
      points(
        x=my_x[meanCNV_reg_signif_dt$dotCol == signif_col],
        y=my_y[meanCNV_reg_signif_dt$dotCol == signif_col],
        col=meanCNV_reg_signif_dt$dotCol[meanCNV_reg_signif_dt$dotCol == signif_col],
        cex=1.2,
        pch = 16
        
      )
      legend("bottomleft",
             legend=c(
               paste0("signif. TADs (p-val<=", tad_signif_thresh, ")\nn=",sum(meanCNV_reg_signif_dt$dotCol == signif_col))
             ),
             pch = 16,
             # col = c(signif_col, not_signif_col),
             col = c(signif_col),
             bty="n"
      )
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
      # mtext(side = 3, text = subTit)
      mtext(side = 3, text = paste0("n=", nrow(meanCNV_reg_signif_dt)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds,"_diffCNV_tadAdjCombPval_log10_signifCol.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      
      my_x <- meanCNV_reg_signif_dt$diff_cnv_samp12
      my_y <- -log10(meanCNV_reg_signif_dt$adjPvalComb)
      
      par(bty="l")
      plot(
        x=my_x,
        y=my_y,
        col=meanCNV_reg_signif_dt$dotCol,
        cex=0.7,
        pch = 16,
        cex.lab=plotCex,
        cex.main=plotCex,
        cex.axis=plotCex,
        xlab = "Mean TAD diff. CNV samp1-samp2",
        ylab = "TAD adj. comb. p-val. [-log10]",
        sub=subTit,
        main = plotTit
      )
      
      points(
        x=my_x[meanCNV_reg_signif_dt$dotCol == signif_col],
        y=my_y[meanCNV_reg_signif_dt$dotCol == signif_col],
        col=meanCNV_reg_signif_dt$dotCol[meanCNV_reg_signif_dt$dotCol == signif_col],
        cex=1.2,
        pch = 16
        
      )
      
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
      # mtext(side = 3, text = subTit)
      mtext(side = 3, text = paste0("n=", nrow(meanCNV_reg_signif_dt)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      


      meanCNV_reg_signif_dt$hicds <- hicds
      meanCNV_reg_signif_dt$exprds <- exprds
      meanCNV_reg_signif_dt
    }# end-for iterating over exprds
    exprds_dt
  } # end-for iterating over hicds
  
  outFile <- file.path(outFolder, "all_cnv_dt.Rdata")
  save(all_cnv_dt , file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_cnv_dt.Rdata")
  all_cnv_dt <- get(load(outFile))
}  

if(exists("meanCNV_reg_signif_dt")) rm("meanCNV_reg_signif_dt")

plotTit <- paste0("all datasets (n=", length(unique(file.path(all_cnv_dt$hicds, all_cnv_dt$exprds))), ")")
subTit <- ""

my_x <- all_cnv_dt$diff_cnv_samp12
my_y <- all_cnv_dt$meanLogFC

outFile <- file.path(outFolder, paste0("allDS_diffCNV_tadLogFC_signifCol.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="l")
plot(
  x=my_x,
  y=my_y,
  col = all_cnv_dt$dotCol,
  cex=0.7,
  pch = 16,
  cex.lab=plotCex,
  cex.main=plotCex,
  cex.axis=plotCex,
  xlab = "Mean TAD diff. CNV samp1-samp2",
  ylab = "TAD meanLogFC",
  sub=subTit,
  main = plotTit
)
points(
  x=my_x[all_cnv_dt$dotCol == signif_col],
  y=my_y[all_cnv_dt$dotCol == signif_col],
  col = all_cnv_dt$dotCol[all_cnv_dt$dotCol == signif_col],
  cex=1.2,
  pch = 16
)
legend("topright",
       legend=c(
         paste0("signif. TADs (p-val<=", tad_signif_thresh, ")")
       ),
       pch = 16,
       # col = c(signif_col, not_signif_col),
       col = c(signif_col),
       bty="n"
)
# mtext(side = 3, text = subTit)
abline(lm(my_y~my_x), lty=2, col="grey")
addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
mtext(side = 3, text = paste0("n=", nrow(all_cnv_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("allDS_diffCNV_tadAdjCombPval_signifCol.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

my_x <- all_cnv_dt$diff_cnv_samp12
my_y <- all_cnv_dt$adjPvalComb

par(bty="l")
plot(
  x=my_x,
  y=my_y,
  col=all_cnv_dt$dotCol,
  cex=0.7,
  pch = 16,
  cex.lab=plotCex,
  cex.main=plotCex,
  cex.axis=plotCex,
  xlab = "Mean TAD diff. CNV samp1-samp2",
  ylab = "TAD adj. comb. p-val.",
  sub=subTit,
  main = plotTit
)

points(
  x=my_x[all_cnv_dt$dotCol == signif_col],
  y=my_y[all_cnv_dt$dotCol == signif_col],
  col=all_cnv_dt$dotCol[all_cnv_dt$dotCol == signif_col],
  cex=1.2,
  pch = 16
  
)

abline(lm(my_y~my_x), lty=2, col="grey")
addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
# mtext(side = 3, text = subTit)
mtext(side = 3, text = paste0("n=", nrow(all_cnv_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("allDS_diffCNV_tadAdjCombPval_log10_signifCol.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

my_x <- all_cnv_dt$diff_cnv_samp12
my_y <- -log10(all_cnv_dt$adjPvalComb)

par(bty="l")
plot(
  x=my_x,
  y=my_y,
  col=all_cnv_dt$dotCol,
  cex=0.7,
  pch = 16,
  cex.lab=plotCex,
  cex.main=plotCex,
  cex.axis=plotCex,
  xlab = "Mean TAD diff. CNV samp1-samp2",
  ylab = "TAD adj. comb. p-val. [-log10]",
  sub=subTit,
  main = plotTit
)

points(
  x=my_x[all_cnv_dt$dotCol == signif_col],
  y=my_y[all_cnv_dt$dotCol == signif_col],
  col=all_cnv_dt$dotCol[all_cnv_dt$dotCol == signif_col],
  cex=1.2,
  pch = 16
  
)

abline(lm(my_y~my_x), lty=2, col="grey")
addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
# mtext(side = 3, text = subTit)
mtext(side = 3, text = paste0("n=", nrow(all_cnv_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))






