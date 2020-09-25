### CREATE TABLE
# -----------------------
#  TAD_ID | genes | meanFC | meanCorr | adj. comb. p-value | signif. at FDR 0.1 ?  | signif. at FDR 0.2  

# Rscript create_final_table_purityFilter.R

# 12.08.2019 -> added start and end positions

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "create_final_table_purityFilter.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script11same_name <- "11sameNbrPF_runEmpPvalCombined/"


source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(reshape2)

pval001_col <- "dodgerblue3"
pval005_col <- "darkorange3"
fdr01_col <- "goldenrod"
fdr02_col <- "darkolivegreen"
fdr_pval_col <- "indianred4"


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 8


pipFolder<- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CREATE_FINAL_TABLE_PURITYFILTER" # "CREATE_FINAL_TABLE_10000" for 10000 permut
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
hicds <- args[1]
exprds <- args[2]

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb_dt$symbol <- as.character(entrez2symb_dt$symbol)
stopifnot(!duplicated(entrez2symb_dt$entrezID))


fdr_thresh1 <- 0.1
fdr_thresh2 <- 0.2

buildTable <- TRUE


all_hicds <- list.files(pipOutFolder)

all_hicds <- all_hicds[!( grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]

all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


hicds = all_hicds[1]

if(buildTable) {
  
  
  all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start building DT for :", hicds, " - ", exprds, "\n")

      comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      if(!file.exists(comb_empPval_file)) return(NULL)
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))

      # ADJUST THE PVAL
      tad_adjCombPval <- p.adjust(comb_empPval, method="BH")
      stopifnot(names(tad_adjCombPval) == names(comb_empPval))

      exprds_hicds_dt <- data.frame(
        hicds=hicds,
        exprds=exprds,
        region = as.character(names(tad_adjCombPval)),
        adjPvalComb = as.numeric(tad_adjCombPval),
        stringsAsFactors = FALSE
      )
      exprds_hicds_dt <- exprds_hicds_dt[order(exprds_hicds_dt$adjPvalComb),]
      stopifnot(!is.na(exprds_hicds_dt))
      exprds_hicds_dt
    } # end-foreach iterating over exprds
    hicds_dt
  } # end-foreach iterating over hicds
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  save(all_result_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile,  "\n"))
  

    
} else {
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  all_result_dt <- eval(parse(text = load(outFile)))
}




#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))









