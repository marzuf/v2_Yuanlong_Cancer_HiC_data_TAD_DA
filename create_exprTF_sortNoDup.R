startTime <- Sys.time()

cat(paste0("> Rscript create_exprTF_sortNoDup.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)



# Rscript create_exprTF_sortNoDup.R c3.tft
# Rscript create_exprTF_sortNoDup.R chea3_lung

# 1 file for all datasets !

### UPDATE sortNoDup 30.06.2018
# -> sort rows of the gene coord. table to ensure alphabetical order of the genes !!
#    so that after melt the 1st gene will always be the 1st in alphabetical order


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
setDir <- ifelse(SSHFS, "~/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

dsIn <- "chea3_lung"
stopifnot(length(args) == 1)
args <- commandArgs(trailingOnly = TRUE)
dsIn <- args[1]

outFold <- file.path("CREATE_EXPRTF_SORTNODUP", dsIn)
dir.create(outFold, recursive=TRUE)

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)

stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)

stopifnot(dsIn == "chea3_lung")


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
} else if(dsIn == "chea3_all") {
  reg_file <- file.path("chea3_all_tissues_TFs_processed.txt")
  reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
  cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
  reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
  # regSymbol targetSymbol targetEntrezID 
  # regSymbol targetSymbol targetEntrezID
  # 1    TFAP2A        RAB26          25837
  # 2    TFAP2A        LRFN4          78999
} else if(dsIn == "chea3_lung") {
  reg_file <- file.path("chea3_lung_TFs_processed.txt")
  reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
  cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
  reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
  
  
  
  reg_dt <- reg_dt[reg_dt$regSymbol %in% names(symb2entrez),]
  cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  reg_dt$regEntrezID <- symb2entrez[reg_dt$regSymbol]
  reg_dt$regEntrezID <- as.character(reg_dt$regEntrezID)
  
  
  
} else if(dsIn == "trrust"){
  reg_file <- file.path("trrust_rawdata.human.tsv")
  reg_dt <- read.delim(reg_file, sep="\t", header=FALSE, stringsAsFactors = FALSE,
                       col.names = c("regSymbol", "targetSymbol", "direction", "ID"))
  cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
  cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
  reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
  # regSymbol targetSymbol  direction       ID targetEntrezID
  # 1      AATF          BAX Repression 22909821            581
  # 2      AATF       CDKN1A    Unknown 17157788           1026
} else if(dsIn == "tftg") {
  reg_file <- file.path("tftg_db_all_processed.txt")
  reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, 
                       col.names=c("regSymbol", "targetEntrezID"))
  cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  # regSymbol targetEntrezID
  # 1 Alx1_homeodomain_DBD_dimeric_13_1             34
  # 2 Alx1_homeodomain_DBD_dimeric_13_1         163882
} else if(dsIn == "motifmap"){
  reg_file <- file.path("MOTIFMAP_ALLGENES/overlapDT_bp.Rdata")
  reg_dt <- get(load(reg_file))
  colnames(reg_dt)[colnames(reg_dt)=="entrezID"] <- "targetEntrezID"
  cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  # refID                                 queryID overlapBp targetEntrezID     regSymbol
  # 1 340348_chr7_128784712_128809534  LM1_RFX1=RFX1_chr7_128800229_128800247        19         340348 LM1_RFX1=RFX1
  # 2       8514_chr1_6052358_6161253      LM1_RFX1=RFX1_chr1_6133012_6133030        19           8514 LM1_RFX1=RFX1
} else if(dsIn == "kegg"){
  reg_file <- file.path("hsa_kegg_entrez.txt")
  reg_dt <- read.delim(reg_file, sep="\t", header=FALSE, stringsAsFactors = FALSE,
                       col.names = c("targetEntrezID", "regSymbol"))
  reg_dt$targetEntrezID <- gsub("hsa:", "",reg_dt$targetEntrezID )
  cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  # targetEntrezID     regSymbol
  # 1          10327 path:hsa00010
  # 2            124 path:hsa00010
}else {
  reg_file <- file.path(paste0(dsIn, ".v7.0.entrez_processed.txt"))
  reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, 
                       col.names=c("regSymbol", "targetEntrezID"))
  cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  # regSymbol targetEntrezID
  # 1 AAANWWTGC_UNKNOWN           4208
  # 2 AAANWWTGC_UNKNOWN            481
}


stopifnot("regSymbol" %in% colnames(reg_dt))
stopifnot("targetEntrezID" %in% colnames(reg_dt))

permutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(permutFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(permutFolder, x)))




foo <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
  cat(paste0("... start: ", hicds, "\n"))
  foo <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
    
    if(dsIn == "chea3_lung") {
      if(! (grepl("lusc", exprds) | grepl("luad", exprds))) return(NULL)
    }
    
    ## => change here, take only the expressed TFs
    de_dt <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "1_runGeneDE", "DE_topTable.Rdata") ))
    geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata") ))
    
    
    ds_reg_dt <- reg_dt[reg_dt$regEntrezID %in% de_dt$genes | reg_dt$regEntrezID %in% geneList,]
    
    cat(paste0("expr TF: nrow(reg_dt)", "\t=\t", nrow(ds_reg_dt), "\n"))
    



  tf_dt <- do.call(rbind, 
  
    by(ds_reg_dt, ds_reg_dt$regSymbol, function(x) {
      reg_entrez <- unique(x$targetEntrezID)
      tf <- unique(x$regSymbol)
      stopifnot(length(tf) == 1)
      all_cmbs <- combn(reg_entrez, m = 2)
      stopifnot(nrow(all_cmbs) == 2)
      stopifnot(ncol(all_cmbs) == 0.5*(length(reg_entrez)-1) * length(reg_entrez))
      
      data.frame(
        gene1_tmp = all_cmbs[1,], 
        gene2_tmp = all_cmbs[2,],
        tf = tf,
        stringsAsFactors = FALSE
      )
      
    })
    )
  
  head(tf_dt)
  
  tf_dt$gene1_tmp <- as.character(tf_dt$gene1_tmp)
  tf_dt$gene2_tmp <- as.character(tf_dt$gene2_tmp)
  tf_dt$gene1 <-as.character(pmin(tf_dt$gene1_tmp, tf_dt$gene2_tmp))
  tf_dt$gene2 <-as.character(pmax(tf_dt$gene1_tmp, tf_dt$gene2_tmp))
  
  all_tf_pairs <- tf_dt[,c("gene1", "gene2", "tf")]
  
  stopifnot(all_tf_pairs$gene1 < all_tf_pairs$gene2)
  
  
  outFile <- file.path(outFold, hicds, exprds, paste0("all_tf_pairs.Rdata"))
  dir.create(dirname(outFile), recursive = TRUE)
  
  save(all_tf_pairs, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))

  } # end-foreach hicds
} # end-foreach exprds
  
######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
