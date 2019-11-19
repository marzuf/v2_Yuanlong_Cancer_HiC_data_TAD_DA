options(scipen=100)

setDir=""

# Rscript plot_signif_TADs_overallTop.R
# Rscript plot_signif_TADs_overallTop.R wt_vs_mut
# Rscript plot_signif_TADs_overallTop.R norm_vs_tumor
# Rscript plot_signif_TADs_overallTop.R subtypes

script_name <- "plot_signif_TADs_overallTop.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

buildData <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(ggplot2)
require(reshape2)


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("plot_lolliTAD_funct.R")
source("../Cancer_HiC_data_TAD_DA/colors_utils.R")
source("subtype_cols.R")



plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.4

nToPlot <- 10
nPlotted <- nToPlot

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)

pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))


outFolder <- "PLOT_SIGNIF_TADS_OVERALLTOP"
dir.create(outFolder, recursive=TRUE)


args <- commandArgs(trailingOnly = TRUE)

toplot <- args[1]

final_DT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

final_DT <- final_DT[order(final_DT$adjPvalComb),]

fileprefix <- ""

if(length(args) == 1) {
  final_DT$cmpType <- all_cmps[final_DT$exprds]
  stopifnot(!is.na(final_DT$cmpType))
  stopifnot(toplot %in% final_DT$cmpType)
  final_DT <- final_DT[final_DT$cmpType %in% toplot,]
  fileprefix <- paste0(toplot, "_")  
}


### BUILD SIGNIF ALONG FDR THRESH
cat("... start retrieving FDR signif. TADs\n")


plotList <- list()

for(i_tad in 1:nToPlot) {
  plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = final_DT$exprds[i_tad],
                                    hicds = final_DT$hicds[i_tad],
                                    all_TADs = final_DT$region[i_tad],
                                    mytitle = paste0(final_DT$exprds[i_tad], "-",final_DT$hicds[i_tad],"-",final_DT$region[i_tad]," (", formatC(final_DT$adjPvalComb[i_tad], format = "e", digits = 2), ")"),
                                    orderByLolli = "startPos")
} # end-for iterating over TADs to plot

outFile <- file.path(outFolder, paste0(fileprefix, "_overallTop_adjCombPvalSignifTADs",  "_nToPlot", nToPlot, ".", plotType ))
mytit <- paste0("overall Top Adj Pval Comb")
all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nToPlot == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
outHeightGG <- min(c(7 * nPlotted/2, 49))
outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
cat("... written: ", outFile, "\n")

      
      
      

##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

