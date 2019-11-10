options(scipen=100)

SSHFS=F

setDir <- "/media/electron"
setDir <- ""

# Rscript presentation_figures4.R

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "presentation_figures4.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(RColorBrewer)
require(reshape2)
require(ggplot2)
require(ggpubr)
# require(gplots)
registerDoMC(4)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb <- setNames(entrez2symb_dt$symbol, entrez2symb_dt$entrezID)


buildTable <- FALSE

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightHeat <- myHeight * 1.8
myWidthHeat <- myWidth * 1.8

myWidthGG <- 12
myHeightGG <- 7

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script11same_name <- "11sameNbr_runEmpPvalCombined"

tieMeth <- "min"


mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

outFolder <- file.path("PRESENTATION_FIGURES4")
dir.create(outFolder, recursive = TRUE)


final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

final_table_DT$cmpType <- all_cmps[paste0(final_table_DT$exprds)]
stopifnot(!is.na(final_table_DT$cmpType))

xx <- final_table_DT[,c("hicds", "exprds")]
xx <- unique(xx)
xx <- xx[order(xx$hicds, xx$exprds),]
xx



tad_signif_col <- "tad_adjCombPval"
gene_signif_col <- "adj.P.Val"

signif_column <- "adjPvalComb"
signifThresh <- 0.01


tad_pval_thresh <- 0.01
gene_pval_thresh <- 0.05
tad_pval_thresh_log10 <- -log10(tad_pval_thresh)
gene_pval_thresh_log10 <- -log10(gene_pval_thresh)
file_suffix <- paste0("tad_pval_thresh", tad_pval_thresh, "_gene_pval_thresh", gene_pval_thresh)

minOverlapBpRatio <- 0.8
minIntersectGenes <- 3

nRegionLolli <- 10

cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> signifThresh\t=\t", signifThresh, "\n"))
cat(paste0("> minOverlapBpRatio\t=\t", minOverlapBpRatio, "\n"))
cat(paste0("> minIntersectGenes\t=\t", minIntersectGenes, "\n"))
cat(paste0("> nRegionLolli\t=\t", nRegionLolli, "\n"))


tad_plot_list <- list()

data_cmpType = "norm_vs_tumor"

plotTit <- setNames(c("normal vs. tumor", "subtypes", "wild-type vs. mutant", "all"),
                    c("norm_vs_tumor", "subtypes", "wt_vs_mut", "")
                    )

for(data_cmpType in c("norm_vs_tumor", "subtypes", "wt_vs_mut", "")) {
  
  if(data_cmpType == "") {
    nDS <- length(unique(paste0(final_table_DT$hicds,final_table_DT$exprds)))
    tmpDS<- length(unique(paste0(final_table_DT$exprds)))
  } else {
    nDS <- length(unique(paste0(final_table_DT$hicds[final_table_DT$cmpType==data_cmpType],final_table_DT$exprds[final_table_DT$cmpType==data_cmpType])))
    tmpDS <- length(unique(paste0(final_table_DT$exprds[final_table_DT$cmpType==data_cmpType])))
  }
  stopifnot(nDS > 0)
  
  
  inFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", data_cmpType)
  stopifnot(dir.exists(inFolder))
  
  inFile <- file.path(inFolder, paste0(ifelse(data_cmpType=="", data_cmpType, paste0(data_cmpType, "_")), 
                      "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio,
                      "_minInterGenes", minIntersectGenes, ".Rdata"))
  stopifnot(file.exists(inFile))
  conserved_signif_tads <- get(load(inFile))
  
  tmp <- lapply(conserved_signif_tads, function(x) unique(basename(dirname(x))))
  exprds_countConserv <- abs(sort(-lengths(tmp)))
  stopifnot(exprds_countConserv >= 1)
  tmp_countConserv_dt <- data.frame(exprds_countConserv)
  tmp_countConserv_dt$region <- factor(rownames(tmp_countConserv_dt), levels=rownames(tmp_countConserv_dt))
  tmp_countConserv_dt$conservRatio <- tmp_countConserv_dt$exprds_countConserv/tmpDS
  stopifnot(tmp_countConserv_dt$conservRatio <= 1)
  
  
  
  countConserv <- abs(sort(-lengths(conserved_signif_tads)))
  stopifnot(countConserv > 1)
  countConserv_dt <- data.frame(countConserv)
  countConserv_dt$region <- factor(rownames(countConserv_dt), levels=rownames(countConserv_dt))
  countConserv_dt$conservRatio <- countConserv_dt$countConserv/nDS
  stopifnot(countConserv_dt$conservRatio <= 1)
  
  
  # > head(countConserv_dt)
  # countConserv              region conservRatio
  # conserved_region_22           10 conserved_region_22       0.6250
  # conserved_region_9             7  conserved_region_9       0.4375
  # conserved_region_10            6 conserved_region_10       0.3750
  # conserved_region_23            6 conserved_region_23       0.3750
  # conserved_region_17            5 conserved_region_17       0.3125
  # conserved_region_30            5 conserved_region_30       0.3125
  # > head(tmp_countConserv_dt)
  # exprds_countConserv              region conservRatio
  # conserved_region_22                   4 conserved_region_22          0.8
  # conserved_region_36                   4 conserved_region_36          0.8
  # conserved_region_9                    3  conserved_region_9          0.6
  # conserved_region_16                   3 conserved_region_16          0.6
  # conserved_region_17                   3 conserved_region_17          0.6
  # conserved_region_30                   3 conserved_region_30          0.6
  # > 
    
  outFile <- file.path(outFolder, paste0(data_cmpType, "_count_conserved_regions.", plotType))
  
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  
  par(bty="l")
  plot(countConserv_dt$countConserv, type="h",
       main = paste0(plotTit[data_cmpType]),
       xlab="Conserved regions",
       ylab="# of datasets",
       cex.main=axisCex,
       cex.axis=axisCex,
       cex.lab=axisCex,
       lwd=5
       )
  legend("topright",
         bty="n",
         legend=paste0("# tot. datasets = ", nDS), cex=1.2)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}







  
######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


