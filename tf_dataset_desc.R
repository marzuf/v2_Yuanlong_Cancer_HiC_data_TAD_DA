startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)
require(dplyr)

# Rscript tf_dataset_desc.R

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- ifelse(plotType=="png", 500, 8)
plotCex <- 1.4

hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAluad_norm_luad"
dsIn <- "chea3_lung"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 0 | length(args) == 2)

if(length(args) == 2) {
  all_hicds <- args[1]
  all_exprds <- args[2]
  
} else {
  all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))
}

stopifnot(dsIn %in% c("crisp",  "trrust", "motifmap", "chea3_all", "chea3_lung"))

outFolder <- file.path(paste0("TF_DATASET_DESC"))#, toupper(dsIn)))
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
  
  
  out_dt <- foreach(dsIn = c("crisp",  "trrust", "motifmap", "chea3_all", "chea3_lung"), .combine='rbind') %do% {  

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
      reg_dt$targetSymbol <- reg_dt$targetEntrezID
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
    
    h_e_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
      

    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata") ))
      
      h_e_reg_dt <- reg_dt[reg_dt$targetEntrezID %in% geneList,]
      
        data.frame(
          hicds = hicds,
          exprds = exprds,
          tf_ds = dsIn,
          nUniqueTarget = length(unique(h_e_reg_dt$targetSymbol)),
          nUniqueReg = length(unique(h_e_reg_dt$regSymbol)),
          stringsAsFactors=FALSE
        )
        
      } # end-foreach iterating exprds
    exprds_dt
      
    } # end-for iterating over hicds
    h_e_dt
  } # end-for each iterating dsIN
  
  
  outFile <- file.path(outFolder, paste0("out_dt.Rdata"))
  save(out_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  # outFile="TF_RATIO_OBS_CHEA3_LUNG/all_ratio_data.Rdata"
  outFile <- file.path(outFolder, paste0("out_dt.Rdata"))
  out_dt <- get(load(outFile))
}  

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path(outFolder, paste0("nUniqueTarget_tfDS.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot_multiDens(
   by(out_dt, out_dt$tf_ds, function(x) log10(x$"nUniqueTarget")),
plotTit = "# unique targets"
)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("nUniqueReg_tfDS.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot_multiDens(
  by(out_dt, out_dt$tf_ds, function(x) log10(x$"nUniqueReg")), plotTit = "# unique reg. elements"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))
