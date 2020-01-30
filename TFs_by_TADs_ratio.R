startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript TFs_by_TADs_ratio.R crisp
# Rscript TFs_by_TADs_ratio.R c3.mir
# Rscript TFs_by_TADs_ratio.R c3.tft
# Rscript TFs_by_TADs_ratio.R c3.all
# Rscript TFs_by_TADs_ratio.R trrust
# Rscript TFs_by_TADs_ratio.R tftg
# Rscript TFs_by_TADs_ratio.R motifmap
# Rscript TFs_by_TADs_ratio.R kegg
# 


plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- ifelse(plotType=="png", 400, 7)
plotCex <- 1.4

nTop <- 10

fontFamily <- "Hershey"

require(ggsci)
top_col <- pal_d3()(2)[1]
last_col <- pal_d3()(2)[2]
# yarrr::transparent("grey", trans.val = .6)
mid_col <- "#BEBEBE66"

tad_signif_thresh <- 0.01
signif_col <- "red"

x_qt_val <- 0.2
y_qt_val <- 0.95


dsIn <- "crisp"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1 | length(args) == 3)
dsIn <- args[1]
if(length(args) == 3) {
  all_hicds <- args[2]
  all_exprds <- args[3]

} else {
  all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))
}

stopifnot(dsIn %in% c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg"))

outFolder <- file.path(paste0("TFS_BY_TADS_RATIO_", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)

stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))


if(buildData){
  allQt_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
  
    cat(paste0("> START - ", hicds,"\n"))
  
    if(dsIn == "crisp") {
      reg_file <- file.path("gene_set_library_crisp_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
      # regSymbol targetSymbol targetEntrezID
      # 2    ARID3A         A1BG              1
      # 3    ARID3A     A1BG-AS1         503538
      # 4    ARID3A         A1CF          29974
      
    } else if(dsIn == "trrust"){
      reg_file <- file.path("trrust_rawdata.human.tsv")
      reg_dt <- read.delim(reg_file, sep="\t", header=FALSE, stringsAsFactors = FALSE,
                           col.names = c("regSymbol", "targetSymbol", "direction", "ID"))
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    } else if(dsIn == "tftg") {
      reg_file <- file.path("tftg_db_all_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, 
                           col.names=c("regSymbol", "targetEntrezID"))
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    } else if(dsIn == "motifmap"){
      reg_file <- file.path("MOTIFMAP_ALLGENES/overlapDT_bp.Rdata")
      reg_dt <- get(load(reg_file))
      colnames(reg_dt)[colnames(reg_dt)=="entrezID"] <- "targetEntrezID"
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    } else if(dsIn == "kegg"){
      reg_file <- file.path("hsa_kegg_entrez.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=FALSE, stringsAsFactors = FALSE,
                           col.names = c("targetEntrezID", "regSymbol"))
      reg_dt$targetEntrezID <- gsub("hsa:", "",reg_dt$targetEntrezID )
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    }else {
      reg_file <- file.path(paste0(dsIn, ".v7.0.entrez_processed.txt"))
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, 
                           col.names=c("regSymbol", "targetEntrezID"))
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    }
    
    g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
    g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    stopifnot(!duplicated(g2t_dt$entrezID))
    g2t_vect <- setNames(g2t_dt$region, g2t_dt$entrezID)
    
    reg_dt <- reg_dt[reg_dt$targetEntrezID %in% g2t_dt$entrezID,]
    cat(paste0("with g2t assignment: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt$targetRegion <- g2t_vect[paste0(reg_dt$targetEntrezID)]
    stopifnot(!is.na(reg_dt))
      
    nbrReg_TADs_dt <- aggregate(regSymbol~targetRegion, data=reg_dt, function(x) length(unique(x)))
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      plotTit <- paste0(hicds, "\n", exprds)
      
      result_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds, ]
      result_dt$tad_rank <- rank(result_dt$adjPvalComb, ties="min")
      result_dt$rev_tad_rank <- rank(-result_dt$adjPvalComb, ties="min")
      
      topTADs <- result_dt$region[result_dt$tad_rank <= nTop]
      lastTADs <- result_dt$region[result_dt$rev_tad_rank <= nTop]
      
      
      signifTADs <- result_dt$region[result_dt$adjPvalComb <= tad_signif_thresh]
      
      
      geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata") ))
      stopifnot(geneList %in% g2t_dt$entrezID)
      gByTAD <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      
      reg_dt <- reg_dt[reg_dt$targetEntrezID %in% geneList,]  # update 08.01.20 -> NEED ALSO TO SUBSET THE REGULATED FEATURES !
      
      # 1) # of genes in TAD
      tad_nGenes_dt <- aggregate(entrezID ~ region, data=gByTAD, FUN=function(x) length(x))
      colnames(tad_nGenes_dt)[colnames(tad_nGenes_dt) == "entrezID"] <- "nGenes"
      stopifnot(tad_nGenes_dt$nGenes >= 3)
      
      # 2) # of genes regulated within TAD
      tad_nRegGenes_dt <- aggregate(targetEntrezID~targetRegion, data=reg_dt, FUN=function(x)length(unique(x)) )
      colnames(tad_nRegGenes_dt)[colnames(tad_nRegGenes_dt) == "targetRegion"] <- "region"
      colnames(tad_nRegGenes_dt)[colnames(tad_nRegGenes_dt) == "targetEntrezID"] <- "nRegGenes"
      
      # 3) # of TFs within TAD
      tad_nTFs_dt <- aggregate(regSymbol~targetRegion, data=reg_dt, FUN=function(x)length(unique(x)) )
      colnames(tad_nTFs_dt)[colnames(tad_nTFs_dt) == "targetRegion"] <- "region"
      colnames(tad_nTFs_dt)[colnames(tad_nTFs_dt) == "regSymbol"] <- "nTFs"
      
      plot_dt <- merge(tad_nTFs_dt, merge(tad_nGenes_dt, tad_nRegGenes_dt,by="region"), by="region")
      
      stopifnot(plot_dt$nRegGenes <= plot_dt$nGenes)
      
      plot_dt$nTFs_byGenes <- plot_dt$nTFs/plot_dt$nGenes
      plot_dt$nRegGenes_byGenes <- plot_dt$nRegGenes/plot_dt$nGenes
      
      stopifnot(!duplicated(plot_dt$region))
      
      plot_dt$dotCols <- ifelse(plot_dt$region %in% topTADs, top_col, 
                                ifelse(plot_dt$region %in% lastTADs, last_col, mid_col))
      
      # save(plot_dt, file="plot_dt.Rdata",version=2);stop("ok");
      
      my_x <- plot_dt$nTFs_byGenes
      my_y <- plot_dt$nRegGenes_byGenes
      my_tads <- plot_dt$region
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_ratio_regGenes_nTFs_allTADs_densplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l", family=fontFamily)
      densplot(
        x = my_x,
        y = my_y,
        main=paste0(plotTit),
        xlab = "# TFs in TAD/# genes in TAD",
        ylab = "# reg. genes in TAD/# genes in TAD",
        pch = 16,
        cex = 0.7,
        cex.axis = plotCex,
        cex.lab = plotCex,
        cex.main = plotCex
      )
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n")
      mtext(side=3, text = paste0(dsIn,  " - n =", nrow(plot_dt)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      

      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_ratio_regGenes_nTFs_allTADs_colplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l", family=fontFamily)
      plot(
        x = my_x,
        y = my_y,
        main=paste0(plotTit),
        xlab = "# TFs in TAD/# genes in TAD",
        ylab = "# reg. genes in TAD/# genes in TAD",
        pch = 16,
        col = plot_dt$dotCols,
        cex = 0.7,
        cex.axis = plotCex,
        cex.lab = plotCex,
        cex.main = plotCex
      )
      
      
      points(
        x = plot_dt$nTFs_byGenes[plot_dt$dotCols == top_col],
        y = plot_dt$nRegGenes_byGenes[plot_dt$dotCols == top_col],
        col = top_col,
        pch = 16,
        cex = 1.2
      )
      
      
      points(
        x = plot_dt$nTFs_byGenes[plot_dt$dotCols == last_col],
        y = plot_dt$nRegGenes_byGenes[plot_dt$dotCols == last_col],
        col = last_col,
        pch = 16,
        cex = 1.2
      )
      
      
      
      
      
      
      
      
      
      legend(
        "bottomright",
        pch=16,
        col = c(top_col, last_col),
        legend=c(paste0("top TADs (n=", sum(plot_dt$region %in% topTADs), ")"),
                 paste0("last TADs (n=", sum(plot_dt$region %in% lastTADs), ")")),
        bty="n"
      )
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n")
      mtext(side=3, text = paste0(dsIn,  " - n =", nrow(plot_dt)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      # signif in red
      plot_dt$dotSignifCols <- ifelse(plot_dt$region %in% signifTADs, signif_col, "grey")
                                
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_ratio_regGenes_nTFs_allTADs_signifcolplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l", family=fontFamily)
      plot(
        x = my_x,
        y = my_y,
        main=paste0(plotTit),
        xlab = "# TFs in TAD/# genes in TAD",
        ylab = "# reg. genes in TAD/# genes in TAD",
        pch = 16,
        col = plot_dt$dotSignifCols,
        cex = 0.7,
        cex.axis = plotCex,
        cex.lab = plotCex,
        cex.main = plotCex
      )
      
      points(
        x = plot_dt$nTFs_byGenes[plot_dt$dotSignifCols == signif_col],
        y = plot_dt$nRegGenes_byGenes[plot_dt$dotSignifCols == signif_col],
        col = signif_col,
        pch = 16,
        cex = 1.2
      )
      
      legend(
        "bottomright",
        pch=16,
        col = c(signif_col),
        legend=c(paste0("signif. TADs (p-val <= ", tad_signif_thresh, ")")),
        bty="n"
      )
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n")
      mtext(side=3, text = paste0(dsIn,  " - n =", nrow(plot_dt)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      
      
      # dot color meanCorr
      
      

      
      
      meancolPlot_dt <- merge(plot_dt, result_dt[,c("region", "meanCorr")], by="region")
      stopifnot(nrow(meancolPlot_dt) == nrow(plot_dt))
      
      #Create a function to generate a continuous color palette
      rbPal <- colorRampPalette(c('red','blue'))
      meancolPlot_dt$dotCols <- rev(rbPal(10))[as.numeric(cut(meancolPlot_dt$meanCorr,breaks = 10))]
      

      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_ratio_regGenes_nTFs_allTADs_meanCorrColplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      par(bty="l", family=fontFamily)
      plot(
        x = my_x,
        y = my_y,
        col = meancolPlot_dt$dotCols, 
        main=paste0(plotTit),
        xlab = "# TFs in TAD/# genes in TAD",
        ylab = "# reg. genes in TAD/# genes in TAD",
        pch = 16,
        cex = 0.7,
        cex.axis = plotCex,
        cex.lab = plotCex,
        cex.main = plotCex
      )
      mtext(side=3, text = paste0(dsIn,  " - n =", nrow(meancolPlot_dt)))
      legend("bottomright",
             pch=c(-1,16,16),
             legend=c("meanCorr", round(meancolPlot_dt$meanCorr[which.min(meancolPlot_dt$meanCorr)],2),
                      round(meancolPlot_dt$meanCorr[which.max(meancolPlot_dt$meanCorr)], 2)),
             col=c("black", meancolPlot_dt$dotCols[which.min(meancolPlot_dt$meanCorr)],
                          meancolPlot_dt$dotCols[which.max(meancolPlot_dt$meanCorr)]),
             bty='n')
      abline(lm(my_y~my_x), lty=2, col="grey")
      addCorr(x = my_x, y = my_y, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      my_x_qt <- quantile(my_x, probs=x_qt_val)
      my_y_qt <- quantile(my_y, probs=y_qt_val)
      
      # my_x_qt <- median(my_x)
      # my_y_qt <- median(my_y)
      
      if(hicds == "Barutcu_MCF-10A_40kb") save(plot_dt, file="my_plot_dt.Rdata", version=2)
      
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_ratio_regGenes_nTFs_allTADs_labPlot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      # plotType="svg"
      # outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_ratio_regGenes_nTFs_allTADs_labPlot.", plotType))
      # do.call(plotType, list(outFile, height=7, width=7))
      par(bty="l", family=fontFamily)
      plot(
        x = my_x,
        y = my_y,
        col = "grey", 
        main=paste0(plotTit),
        xlab = "# TFs in TAD/# genes in TAD",
        ylab = "# reg. genes in TAD/# genes in TAD",
        pch = 16,
        cex = 0.2,
        cex.axis = plotCex,
        cex.lab = plotCex,
        cex.main = plotCex
      )
      if(any(my_x <= my_x_qt & my_y >= my_y_qt))
        text(
          x = my_x[my_x <= my_x_qt & my_y >= my_y_qt],
          y = my_y[my_x <= my_x_qt & my_y >= my_y_qt],
          labels = my_tads[my_x <= my_x_qt & my_y >= my_y_qt],
          cex = 0.8
        )
      legend("bottomright", legend=c(paste0("x<=", round(my_x_qt, 2), " (", x_qt_val, "-qt)"), 
                                     paste0("y>=", round(my_y_qt, 2), " (", y_qt_val, "-qt)"),
                                     paste0("(in qt: ", sum(my_x <= my_x_qt & my_y >= my_y_qt),  ")")
                                     ), bty="n") 
      
      
      qtTADs <- my_tads[my_x <= my_x_qt & my_y >= my_y_qt]
      
      
      signifInQt <- sum(qtTADs %in% signifTADs)/ length(signifTADs)
      
      legend("topright",
             legend = c(paste0("signif. in qt: ", sum(qtTADs %in% signifTADs), "/", length(signifTADs), "\n(", round(signifInQt*100, 2), "%)" )),
             bty="n"
      )
      
      
      mtext(side=3, text = paste0(dsIn,  " - n =", nrow(meancolPlot_dt)))
      # legend("bottomright",
      #        pch=c(-1,16,16),
      #        legend=c("meanCorr", round(meancolPlot_dt$meanCorr[which.min(meancolPlot_dt$meanCorr)],2),
      #                 round(meancolPlot_dt$meanCorr[which.max(meancolPlot_dt$meanCorr)], 2)),
      #        col=c("black", meancolPlot_dt$dotCols[which.min(meancolPlot_dt$meanCorr)],
      #              meancolPlot_dt$dotCols[which.max(meancolPlot_dt$meanCorr)]),
      #        bty='n')
      # abline(lm(my_y~my_x), lty=2, col="grey")
      # addCorr(x = my_x, y = my_y, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        qtTADs = paste0(qtTADs, collapse=","),
        qtSignifTADs = paste0(qtTADs[qtTADs %in% signifTADs], collapse=","),
        nInQt = length(qtTADs),
        ratioSignifInQt = round(signifInQt,4),
        stringsAsFactors = FALSE
      )
    }
    exprds_dt
  } # end-for iterating over hicds
  outFile <- file.path(outFolder, paste0("allQt_xQt", x_qt_val, "_yQt", y_qt_val, "_dt.Rdata"))
  save(allQt_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0("allQt_xQt", x_qt_val, "_yQt", y_qt_val, "_dt.Rdata"))
  allQt_dt <- get(load(outFile))
}  

# allQt_dt <- get(load("TFS_BY_TADS_RATIO_C3.TFT/allQt_xQt0.2_yQt0.95_dt.Rdata"))

nDS <- length(unique(file.path(allQt_dt$hicds, allQt_dt$exprds)))

outFile <- file.path(outFolder, paste0(dsIn, "_all_ds_ratioSignifTADsInQt_xQt", x_qt_val, "_yQt", y_qt_val, "_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(allQt_dt$ratioSignifInQt,
        main = paste0("Ratio of signif. TADs (x<=", x_qt_val, "-qt & y>=", y_qt_val, "-qt)"),
        ylab = "Ratio signif. TADs in qt"
        )
mtext(side=3, text = paste0(dsIn,  " - n =", nDS))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0(dsIn, "_all_ds_nTADsInQt_xQt", x_qt_val, "_yQt", y_qt_val, "_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(allQt_dt$nInQt,
        main = paste0("# TADs in qt (x<=", x_qt_val, "-qt & y>=", y_qt_val, "-qt)"),
        ylab = "# TADs in qt"
)
mtext(side=3, text = paste0(dsIn,  " - n =", nDS))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



load("TFS_BY_TADS_RATIO_C3.TFT/allQt_xQt0.2_yQt0.95_dt.Rdata")
dt=allQt_dt
dt = dt[order(dt$ratioSignifInQt, decreasing=T),]          
write.table(dt, file="TFS_BY_TADS_RATIO_C3.TFT/allQt_xQt0.2_yQt0.95_dt_result.txt", sep="\t", col.names=T, row.names=F, quote=F, append=F)


#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))
