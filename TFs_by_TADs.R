startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript TFs_by_TADs.R crisp 
# Rscript TFs_by_TADs.R c3.mir
# Rscript TFs_by_TADs.R c3.mir LG1_40kb TCGAluad_norm_luad
# Rscript TFs_by_TADs.R c3.tft
# Rscript TFs_by_TADs.R c3.all
# Rscript TFs_by_TADs.R trrust
# Rscript TFs_by_TADs.R tftg
# Rscript TFs_by_TADs.R motifmap LG1_40kb TCG
# Rscript TFs_by_TADs.R motifmap 
# Rscript TFs_by_TADs.R kegg LG1_40kb

plotType <- "png"
myHeight <- 400
myWidth <- 600
plotCex <- 1.4


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

outFolder <- file.path(paste0("TFS_BY_TADS_", toupper(dsIn)))
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



if(buildData){
  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
  
    cat(paste0("> START - ", hicds,"\n"))
  
    # tadFile <- file.path(hicds, "genes2tad", "all_assigned_regions.txt")
    # tad_DT <- read.delim(tadFile, header=FALSE, stringsAsFactors = FALSE, col.names = c("chromo", "region", "start", "end"))
    # tad_DT <- tad_DT[grepl("_TAD", tad_DT$region),]
    
    if(dsIn == "crisp") {
      reg_file <- file.path("gene_set_library_crisp_processed.txt")
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
    
    ds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      pval_comb_file <- file.path("PIPELINE/OUTPUT_FOLDER/LG1_40kb/TCGAluad_norm_luad/11sameNbr_runEmpPvalCombined/emp_pval_combined.Rdata") 
      pval_comb <- get(load(pval_comb_file))
      adj_pval_comb <- p.adjust(pval_comb, method="BH")
      pval_dt <- data.frame(targetRegion = names(adj_pval_comb), adjPvalComb=adj_pval_comb, stringsAsFactors = FALSE)
      
      
      tf_pval_dt <- merge(nbrReg_TADs_dt, pval_dt, by="targetRegion")
      stopifnot(!duplicated(tf_pval_dt$targetRegion))
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_pval_vs_nReg.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      densplot(
        main=paste0(hicds, " - ", exprds),
        xlab = "adj. comb. pval [-log10]",
        ylab = "# reg. elements",
        x = -log10(tf_pval_dt$adjPvalComb),
        y = tf_pval_dt$regSymbol,
        cex.lab = plotCex,
        cex.axis = plotCex
      )
      mtext(side=3, text = paste0(dsIn, " - # TADs = ", nrow(tf_pval_dt)))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      # tf_pval_dt$hicds <- hicds
      # tf_pval_dt$exprds <- exprds
      tf_pval_dt
    }
    ds_dt
  } # end-for iterating over hicds
  outFile <- file.path(outFolder, "all_dt.Rdata")
  save(all_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))

} else {
  inFile <- file.path(outFolder, "all_dt.Rdata")
  all_dt <- get(load(inFile))
}  


outFile <- file.path(outFolder, paste0("allDS_pval_vs_nReg.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  main=paste0("all DS"),
  xlab = "adj. comb. pval [-log10]",
  ylab = "# reg. elements",
  x = -log10(all_dt$adjPvalComb),
  y = all_dt$regSymbol,
  cex.lab = plotCex,
  cex.axis = plotCex
)
mtext(side=3, text = paste0(dsIn, " - # TADs = ", nrow(all_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))
