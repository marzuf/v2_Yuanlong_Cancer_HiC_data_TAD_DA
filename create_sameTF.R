startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
require(doMC)
registerDoMC(40)

# Rscript create_sameTF.R 

plotType <- "png"
myHeight <- 400
myWidth <- 600
plotCex <- 1.4


all_db <- c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg")
all_db <- c("crisp", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg")

# all_db="c3.mir"

foo <- foreach(dsIn = all_db) %dopar% {
  
  outFolder <- file.path(paste0("CREATE_SAMETF_", toupper(dsIn)))
  dir.create(outFolder, recursive = TRUE)
  
  
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
  
  # build the sameTF dt
  sameTF_dt <- do.call(rbind, by(reg_dt, reg_dt$regSymbol, function(x) {
    all_regGenes <- unique(x$targetEntrezID)
    if(length(all_regGenes) == 1) return(NULL)
    all_cmbs <- combn(all_regGenes, m=2)
    genes1 <- as.character(all_cmbs[1,])
    genes2 <- as.character(all_cmbs[2,])
    data.frame(
      gene1 = pmin(genes1,genes2),
      gene2 = pmax(genes1,genes2),
      sameTF = 1,
      stringsAsFactors = FALSE
    )
  }))
  rownames(sameTF_dt) <- NULL
  
  outFile <- file.path(outFolder, "sameTF_dt.Rdata")
  save(sameTF_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
}

