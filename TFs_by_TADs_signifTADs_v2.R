startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript TFs_by_TADs_signifTADs_v2.R crisp
# Rscript TFs_by_TADs_signifTADs_v2.R c3.mir
# Rscript TFs_by_TADs_signifTADs_v2.R c3.tft
# Rscript TFs_by_TADs_signifTADs_v2.R c3.all
# Rscript TFs_by_TADs_signifTADs_v2.R trrust
# Rscript TFs_by_TADs_signifTADs_v2.R tftg
# Rscript TFs_by_TADs_signifTADs_v2.R motifmap
# Rscript TFs_by_TADs_signifTADs_v2.R kegg
# Rscript TFs_by_TADs_signifTADs_v2.R chea3_all
# Rscript TFs_by_TADs_signifTADs_v2.R chea3_lung

# 


plotCex <- 1.4

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 9)
plotCex <- 1.4

nTop <- 10

fontFamily <- "Hershey"

require(ggsci)
top_col <- pal_d3()(2)[1]
last_col <- pal_d3()(2)[2]
# yarrr::transparent("grey", trans.val = .6)
mid_col <- "#BEBEBE66"

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
  all_hicds <- all_hicds[!grepl("_RANDOM", all_hicds)]
  all_hicds <- all_hicds[!grepl("_PERMUT", all_hicds)]

  all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))
}

stopifnot(dsIn %in% c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg", "chea3_all", "chea3_lung"))

outFolder <- file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)

buildData <- FALSE

tad_signif_thresh <- 0.01


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
  nRegFeat_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
  
    cat(paste0("> START - ", hicds,"\n"))
  
    if(dsIn == "crisp") {
      reg_file <- file.path("gene_set_library_crisp_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    } else if(dsIn == "chea3_all") {
      reg_file <- file.path("chea3_all_tissues_TFs_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    } else if(dsIn == "chea3_lung") {
      reg_file <- file.path("chea3_lung_TFs_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
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
    
    hicds_reg_dt <- reg_dt  
    rm("reg_dt")

    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      
      if(dsIn == "chea3_lung") {
        if(! (grepl("lusc", exprds) | grepl("luad", exprds))) return(NULL)
      }
      
      plotTit <- paste0(hicds, "\n", exprds)
      
      result_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds, ]
      result_dt$tad_rank <- rank(result_dt$adjPvalComb, ties="min")
      result_dt$rev_tad_rank <- rank(-result_dt$adjPvalComb, ties="min")
      
      geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata") ))
      stopifnot(geneList %in% g2t_dt$entrezID)
      gByTAD <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      
      
      nGbyT <- setNames(as.numeric(table(g2t_dt$region)), names(table(g2t_dt$region)))
      
      reg_dt <- hicds_reg_dt[hicds_reg_dt$targetEntrezID %in% geneList,]  # update 08.01.20 -> NEED ALSO TO SUBSET THE REGULATED FEATURES !
      
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
      
      
      plot_dt$hicds <- hicds
      plot_dt$exprds <- exprds
      
      signif_plot_dt <- merge(plot_dt, result_dt[,c("hicds", "exprds", "region", "adjPvalComb")], by=c("hicds", "exprds", "region"))
      
      stopifnot(signif_plot_dt$region %in% names(nGbyT))
      signif_plot_dt$nGenes <- nGbyT[paste0(signif_plot_dt$region)]
      stopifnot(!is.na(signif_plot_dt$nGenes))
      
      
      signif_plot_dt$signif_lab <- ifelse(signif_plot_dt$adjPvalComb <= tad_signif_thresh, "signif.", "not signif.")
      
     
      nGenes_signif <- signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "signif."]
      nGenes_notSignif <- signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "not signif."]

      
      nTFs_signif <- signif_plot_dt$nTFs[signif_plot_dt$signif_lab == "signif."]
      nTFs_notSignif <- signif_plot_dt$nTFs[signif_plot_dt$signif_lab == "not signif."]
      
      nRegGenes_signif <- signif_plot_dt$nRegGenes[signif_plot_dt$signif_lab == "signif."]
      nRegGenes_notSignif <- signif_plot_dt$nRegGenes[signif_plot_dt$signif_lab == "not signif."]
      
      nTFsOVERnGenes_signif <- signif_plot_dt$nTFs[signif_plot_dt$signif_lab == "signif."]/signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "signif."]
      nTFsOVERnGenes_notSignif <- signif_plot_dt$nTFs[signif_plot_dt$signif_lab == "not signif."]/signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "not signif."]
      
      nRegGenesOVERnGenes_signif <- signif_plot_dt$nRegGenes[signif_plot_dt$signif_lab == "signif."]/signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "signif."]
      nRegGenesOVERnGenes_notSignif <- signif_plot_dt$nRegGenes[signif_plot_dt$signif_lab == "not signif."]/signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "not signif."]
      
      data.frame(
        hicds = hicds,
        exprds = exprds, 
        mean_nTFs_signif = mean(nTFs_signif),
        mean_nTFs_notSignif = mean(nTFs_notSignif),
        mean_nRegGenes_signif = mean(nRegGenes_signif),
        mean_nRegGenes_notSignif = mean(nRegGenes_notSignif),
        mean_nTFsOVERnGenes_signif = mean(nTFsOVERnGenes_signif),
        mean_nTFsOVERnGenes_notSignif = mean(nTFsOVERnGenes_notSignif),
        mean_nRegGenesOVERnGenes_signif = mean(nRegGenesOVERnGenes_signif),
        mean_nRegGenesOVERnGenes_notSignif = mean(nRegGenesOVERnGenes_notSignif),
        mean_nGenes_signif = mean(nGenes_signif),
        mean_nGenes_notSignif = mean(nGenes_notSignif),
        stringsAsFactors = FALSE
      )
    }# end-for iterating over exprds
    exprds_dt
  } # end-for iterating over hicds
  outFile <- file.path(outFolder, "nRegFeat_dt.Rdata")  
  save(nRegFeat_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "nRegFeat_dt.Rdata")  
  nRegFeat_dt <- get(load(outFile))
}  
# load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/nRegFeat_dt.Rdata")
outFile <- file.path(outFolder, paste0("nRegFeat_boxplot_allDS.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(mar=par()$mar+c(9,0,0,0))
boxplot(nRegFeat_dt[,!colnames(nRegFeat_dt) %in% c("hicds", "exprds")], las=2, main=paste0("all obs. ds (n=", length(unique(file.path(nRegFeat_dt$hicds, nRegFeat_dt$exprds))),")"),  cex.axis=0.8)
mtext(side=3, text = paste0(dsIn))
cat(paste0("... written: ", outFile, "\n"))

# load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/nRegFeat_dt.Rdata")

keepCols <- c("mean_nTFs_signif", "mean_nTFs_notSignif", "mean_nGenes_signif", "mean_nGenes_notSignif", "mean_nTFsOVERnGenes_signif", "mean_nTFsOVERnGenes_notSignif")

outFile <- file.path(outFolder, paste0("nRegFeat_boxplot_allDS_keepCols.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(mar=par()$mar+c(9,0,0,0))
boxplot(nRegFeat_dt[, keepCols], las=2, main=paste0("all obs. ds (n=", length(unique(file.path(nRegFeat_dt$hicds, nRegFeat_dt$exprds))),")"),  cex.axis=0.8)
mtext(side=3, text = paste0(dsIn))
cat(paste0("... written: ", outFile, "\n"))



#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))
