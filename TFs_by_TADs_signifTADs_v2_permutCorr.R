startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R crisp
# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R c3.mir
# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R c3.tft
# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R c3.all
# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R trrust
# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R tftg
# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R motifmap
# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R kegg
# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R chea3_all
# Rscript TFs_by_TADs_signifTADs_v2_permutCorr.R chea3_lung

# 


plotCex <- 1.4

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 500, 8)
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
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))
}

stopifnot(dsIn %in% c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg", "chea3_all", "chea3_lung"))

outFolder <- file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_PERMUTCORR_", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)

buildData <- FALSE




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
  permutCorr_nRegFeat_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
  
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
    
    hicds_reg_dt <- reg_dt 
    rm("reg_dt")
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      
    sample_data <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "5sameNbr_runPermutationsCorr", "sample_around_TADs_sameNbr.Rdata") ))
    cat(paste0("... loaded ...\n"))
    
    
    g2t_dt_left <- do.call(rbind, lapply(1:length(sample_data), function(x) {
      if(sample_data[[x]][["nGenes_left"]] == 0) return(NULL)
      data.frame(
        # region = names(sample_data)[x],
        region = paste0(names(sample_data)[x], "_left"),  # otherwise a TAD will have 2x the # of genes
        entrezID = sample_data[[x]][["genes_left"]],
        stringsAsFactors = FALSE)
    }
    ))
    
    g2t_dt_right <- do.call(rbind, lapply(1:length(sample_data), function(x) {
      if(sample_data[[x]][["nGenes_right"]] == 0) return(NULL)
      data.frame(
        # region = names(sample_data)[x],
        region = paste0(names(sample_data)[x], "_right"),  # otherwise a TAD will have 2x the # of genes
        entrezID = sample_data[[x]][["genes_right"]],
        stringsAsFactors = FALSE)
    }
    ))
    
    g2t_dt <- rbind(g2t_dt_left,g2t_dt_right)
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
  
    g2t_vect <- setNames(g2t_dt$region, g2t_dt$entrezID)
    
    reg_dt <- hicds_reg_dt[hicds_reg_dt$targetEntrezID %in% g2t_dt$entrezID,]
    cat(paste0("with g2t assignment: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt$targetRegion <- g2t_vect[paste0(reg_dt$targetEntrezID)]
    stopifnot(!is.na(reg_dt))
      
    nbrReg_TADs_dt <- aggregate(regSymbol~targetRegion, data=reg_dt, function(x) length(unique(x)))
    

      
      if(dsIn == "chea3_lung") {
        if(! (grepl("lusc", exprds) | grepl("luad", exprds))) return(NULL)
      }
      
      plotTit <- paste0(hicds, "\n", exprds)
      
      result_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds, ]
      result_dt$tad_rank <- rank(result_dt$adjPvalComb, ties="min")
      result_dt$rev_tad_rank <- rank(-result_dt$adjPvalComb, ties="min")
      
      geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata") ))
      # stopifnot(geneList %in% g2t_dt$entrezID) # not for permut
      gByTAD <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      
      
      nGbyT <- setNames(as.numeric(table(g2t_dt$region)), names(table(g2t_dt$region)))
      
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
      
      
      plot_dt$hicds <- hicds
      plot_dt$exprds <- exprds
      
      permutCorr_plot_dt <- plot_dt
      
      stopifnot(permutCorr_plot_dt$region %in% names(nGbyT))
      permutCorr_plot_dt$nGenes <- nGbyT[paste0(permutCorr_plot_dt$region)]
      stopifnot(!is.na(permutCorr_plot_dt$nGenes))
      
      
      
      
     
      nGenes_permutCorr <- permutCorr_plot_dt$nGenes

      
      nTFs_permutCorr <- permutCorr_plot_dt$nTFs
      
      nRegGenes_permutCorr <- permutCorr_plot_dt$nRegGenes
      
      nTFsOVERnGenes_permutCorr <- permutCorr_plot_dt$nTFs
      
      nRegGenesOVERnGenes_permutCorr <- permutCorr_plot_dt$nRegGenes
      
      data.frame(
        hicds = hicds,
        exprds = exprds, 
        mean_nTFs_permutCorr = mean(nTFs_permutCorr),
        mean_nRegGenes_permutCorr = mean(nRegGenes_permutCorr),
        mean_nTFsOVERnGenes_permutCorr = mean(nTFsOVERnGenes_permutCorr),
        mean_nRegGenesOVERnGenes_permutCorr = mean(nRegGenesOVERnGenes_permutCorr),
        mean_nGenes_permutCorr = mean(nGenes_permutCorr),
        stringsAsFactors = FALSE
      )
    }# end-for iterating over exprds
    exprds_dt
  } # end-for iterating over hicds
  outFile <- file.path(outFolder, "permutCorr_nRegFeat_dt.Rdata")  
  save(permutCorr_nRegFeat_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "permutCorr_nRegFeat_dt.Rdata")  
  permutCorr_nRegFeat_dt <- get(load(outFile))
}  
# load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/permutCorr_nRegFeat_dt.Rdata")
outFile <- file.path(outFolder, paste0("permutCorr_nRegFeat_boxplot_allDS.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(mar=par()$mar+c(9,0,0,0))
boxplot(permutCorr_nRegFeat_dt[,!colnames(permutCorr_nRegFeat_dt) %in% c("hicds", "exprds")],
        las=2, 
        main=paste0("all ds (n=", length(unique(file.path(permutCorr_nRegFeat_dt$hicds, permutCorr_nRegFeat_dt$exprds))),")"),  
        cex.main = plotCex, cex.lab = plotCex,
        cex.axis=0.8)
mtext(side=3, text = paste0("permutCorr - ", dsIn))
cat(paste0("... written: ", outFile, "\n"))

# load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/permutCorr_nRegFeat_dt.Rdata")

keepCols <- c("mean_nTFs_permutCorr", "mean_nGenes_permutCorr", "mean_nTFsOVERnGenes_permutCorr")

outFile <- file.path(outFolder, paste0("permutCorr_nRegFeat_boxplot_allDS_keepCols.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(mar=par()$mar+c(9,0,0,0))
boxplot(permutCorr_nRegFeat_dt[, keepCols], las=2, 
        main=paste0("all ds (n=", length(unique(file.path(permutCorr_nRegFeat_dt$hicds, permutCorr_nRegFeat_dt$exprds))),")"),  
        cex.main = plotCex, cex.lab = plotCex,
        cex.axis=0.8)
mtext(side=3, text = paste0("permutCorr - ", dsIn))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))
