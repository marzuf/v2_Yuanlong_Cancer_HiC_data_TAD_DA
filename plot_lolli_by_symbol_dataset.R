options(scipen=100)

SSHFS=F

# Rscript plot_lolli_by_symbol_dataset.R <gene_symbol> <hicds> <exprds>
# Rscript plot_lolli_by_symbol_dataset.R MMP2 LG1_40kb TCGAlusc_norm_lusc

#norm vs. tumor limma missed: JMJD8 MT1L MAPK7 FCGR2C NUP93, SORD, SLAMF7, LILRA2, ACTL10, RRP12, BRICD5, GSTM4, PGAM1, RHBDL1, FBXL16, MT1A, MT1E, MT1F, MT1G, MT1H
gene_symbol <- "IFIT3"
gene_symbol <- "JMJD8"
gene_entrez <- "3437"
gene_symbol <- "MMP2"
hicds <- "LG1_40kb"
exprds <- "TCGAlusc_norm_lusc"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
gene_symbol <- args[1]
hicds <- args[2]
exprds <- args[3]

script_name <- "plot_lolli_by_symbol_dataset.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

require(gridExtra)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("plot_lolliTAD_funct.R")
source("my_heatmap.2.R")

buildTable <- TRUE

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightHeat <- myHeight * 1.8
myWidthHeat <- myWidth * 1.8

myWidthGG <- 12
myHeightGG <- 12


source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

script0_name <- "0_prepGeneData"

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

stopifnot(gene_symbol %in% entrez2symb)
gene_entrez <- names(entrez2symb[entrez2symb == gene_symbol])

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- file.path("PLOT_LOLLI_BY_SYMBOL_DATASET")
dir.create(outFolder, recursive = TRUE)

#all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))
all_datasets <- file.path(hicds, exprds)

final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))

signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
final_dt[, paste0(signifcol)] <- final_dt[, paste0(signif_column)] <= signifThresh

### retrieve the TAD in which the gene is located
ds=all_datasets[1]
ds=all_datasets[2]
# return a value if present in the dataset 
gene_tads <- foreach(ds = all_datasets, .combine='c') %dopar% {
  hicds_file <- file.path(dirname(ds), "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(hicds_file))
  g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
  g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
  stopifnot(nrow(g2t_dt) > 0 )
  gene_tad <- g2t_dt$region[g2t_dt$entrezID == gene_entrez]
  if(length(gene_tad) == 0) return(NULL)
  
  match_idx <- which(final_dt$hicds == dirname(ds) & final_dt$exprds == basename(ds) & final_dt$region == gene_tad)
  # might be zero if not in pipeline_region !
  if(length(match_idx) == 0) return(NULL)

  return(setNames(gene_tad, ds))

#  if(final_dt[final_dt$hicds == dirname(ds) & final_dt$exprds == basename(ds) & final_dt$region == gene_tad, paste0(signifcol)]) {
#    return(setNames(gene_tad, ds))
#  } else {
#    return(NULL)
#  }
}

if(length(gene_tads) == 0) {
  stop("--- gene not used in pipeline\n")
}

stopifnot(length(gene_tads) > 0)


#######################################################################################################################################
######################################################################################################################### PLOT - lolli plot nTop regions all DS
#######################################################################################################################################                   
nTotDS <- length(all_datasets)
nTotExprds <- length(unique(basename(all_datasets)))
nDS <- length(gene_tads)
nExprds <- length(unique(basename(names(gene_tads))))

nPlotted <- length(gene_tads)

outHeightGG <- min(c(7 * nPlotted/2, 49))
outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)

plotList <- list()
i_tad=1

# foo <- foreach(i_tad = 1:nPlotted) %dopar% {
  for(i_tad in 1:nPlotted) {
  
  hicds <- dirname(names(gene_tads[i_tad]))
  exprds <- basename(names(gene_tads[i_tad]))
  tad <- as.character(gene_tads[i_tad])
  
  mytit <- paste0( hicds, " - ", exprds, " - ", tad)
  
  
  plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                        hicds = hicds,
                                        all_TADs = tad,
                                        orderByLolli = "startPos", mytitle=mytit)
} # end-for iterating over TADs to plot


outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", gene_symbol, "_", gene_entrez, "_all_ds_signif_TADs_lolli.", plotType))

mytit <- paste0("All signif. TADs containing ", gene_symbol, " (", gene_entrez, ")", " (", nDS, "/", nTotDS, "; ", nExprds, "/", nTotExprds, ")")
all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
outHeightGG <- min(c(7 * nPlotted/2, 49))
outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)

ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
cat("... written: ", outFile, "\n")


## plot by cmpType
#cmp = unique(all_cmps)[1]
#for(cmp in unique(all_cmps)) {
#  
#  
#  toKeepDS <- names(all_cmps[all_cmps==cmp])
#  
#  cmp_datasets <- all_datasets[basename(all_datasets) %in% toKeepDS]
#  cmp_tads <- gene_tads[basename(names(gene_tads)) %in% toKeepDS]
#  
#  if(length(cmp_tads) == 0) next
#  
#  # cmp_datasets <- names(cmp_tads)
#  
#  cmp_nTotDS <- length(cmp_datasets)
#  cmp_nTotExprds <- length(unique(basename(cmp_datasets)))
#  cmp_nDS <- length(cmp_tads)
#  cmp_nExprds <- length(unique(basename(names(cmp_tads))))
#  
#  cmp_nPlotted <- length(cmp_tads)
#  
#  outHeightGG <- min(c(7 * cmp_nPlotted/2, 49))
#  outHeightGG <- ifelse(cmp_nPlotted < 3, outHeightGG*1.5,outHeightGG)
#  outWidthGG <- ifelse(cmp_nPlotted == 1, 20/2, 20)
#  
#  plotList <- list()
#  i_tad=1
#  
#  for(i_tad in 1:cmp_nPlotted) {
#  # foo <- foreach(i_tad = 1:cmp_nPlotted) %dopar% {
#    
#    hicds <- dirname(names(cmp_tads[i_tad]))
#    exprds <- basename(names(cmp_tads[i_tad]))
#    tad <- as.character(cmp_tads[i_tad])
#    
#    mytit <- paste0( hicds, " - ", exprds, " - ", tad)
#    
#    
#    plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
#                                          hicds = hicds,
#                                          all_TADs = tad,
#                                          orderByLolli = "startPos", mytitle=mytit)
#  } # end-for iterating over TADs to plot
#  
#  
#  outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", gene_symbol, "_", gene_entrez, "_", cmp, "_signif_TADs_lolli.", plotType))
#  
#  mytit <- paste0(cmp, " signif. TADs containing ", gene_symbol, " (", gene_entrez, ")", " (", cmp_nDS, "/", cmp_nTotDS, "; ", cmp_nExprds, "/", cmp_nTotExprds, ")")
#  all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(cmp_nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
#  outHeightGG <- min(c(7 * cmp_nPlotted/2, 49))
#  outHeightGG <- ifelse(cmp_nPlotted < 3, outHeightGG*1.5,outHeightGG)
#  outWidthGG <- ifelse(cmp_nPlotted == 1, 20/2, 20)
#  
#  ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
#  cat("... written: ", outFile, "\n")
#  
#  
#  
#}



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

