startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript TFs_stats.R

outFolder <- "TFS_STATS"
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- 400
myWidth <- 600
plotCex <- 1.4


dsIn <- "crisp"


setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)



stat_dt <- foreach(dsIn = c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg"), .combine='rbind') %dopar% {
  
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
  
  nTFs <- length(unique(reg_dt$regSymbol))
  
  nReg_by_TF <- setNames(as.numeric(table(reg_dt$regSymbol)), names(table(reg_dt$regSymbol)))
  
  data.frame(
    TF_set = dsIn,
    nTFs = nTFs,
    meanReg = round(mean(nReg_by_TF), 4),
    medReg = round(median(nReg_by_TF), 4),
    stringsAsFactors = FALSE
    )
}

outFile <- file.path(outFolder, "stat_dt.Rdata")
save(stat_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
     
outFile <- file.path(outFolder, "stat_dt.txt")
write.table(stat_dt, col.names=TRUE, row.names = FALSE, sep="\t", quote=F, append=F, file = outFile)
cat(paste0("... written: ", outFile, "\n"))



  
