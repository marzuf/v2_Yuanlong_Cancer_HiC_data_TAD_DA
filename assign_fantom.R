startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript assign_fantom.R
# 


plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- ifelse(plotType=="png", 400, 7)
plotCex <- 1.4

nTop <- 10

fontFamily <- "Hershey"




outFolder <- file.path(paste0("ASSIGN_FANTOM_", toupper(dsIn)))
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


fantom_dt <- read.delim("41586_2014_BFnature13182_MOESM86_ESM_A promoter-level mammalian expression atlas_processed.csv", header=TRUE)