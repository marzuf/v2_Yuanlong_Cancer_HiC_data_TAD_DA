startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R crisp
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R c3.mir
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R c3.tft
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R c3.all
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R trrust
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R tftg
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R motifmap
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R kegg
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R chea3_all
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_median.R chea3_lung

# 


plotCex <- 1.4

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 9)
plotCex <- 1.4

nTop <- 10

fontFamily <- "Hershey"

require(ggsci)
top_col <- pal_d3()(2)[1]
last_col <- pal_d3()(2)[2]
# yarrr::transparent("grey", trans.val = .6)
mid_col <- "#BEBEBE66"

x_qt_val <- 0.2
y_qt_val <- 0.95


dsIn <- " c3.tft"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1 | length(args) == 3)
dsIn <- args[1]
if(length(args) == 3) {
  all_hicds <- args[2]
  all_exprds <- args[3]
  
} else {
  all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds <- all_hicds[!grepl("_RANDOM", all_hicds)]
all_hicds <- all_hicds[!grepl("_PERMUT", all_hicds)]

  all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))
}

stopifnot(dsIn %in% c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg", "chea3_all", "chea3_lung"))

outFolder <- file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_CMP_OBS_PERMUT_MEDIAN_", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

obs_data <- get(load(file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_MEDIAN_", toupper(dsIn)), "nRegFeat_dt.Rdata")))
permutCorr_data <- get(load(file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_PERMUTCORR_MEDIAN_", toupper(dsIn)), "permutCorr_nRegFeat_dt.Rdata")))
permutG2t_data <- get(load(file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_PERMUTG2T_MEDIAN_", toupper(dsIn)), "permutG2t_nRegFeat_dt.Rdata")))

plot_dt <- merge(obs_data, merge(permutG2t_data, permutCorr_data, by =c("hicds", "exprds")), by =c("hicds", "exprds"))

plot_cols <- c(
  "median_nGenes_signif"  ,"median_nGenes_notSignif", "median_nGenes_permutCorr","median_nGenes_permutG2t",
  "median_nRegGenes_signif","median_nRegGenes_notSignif","median_nRegGenes_permutCorr","median_nRegGenes_permutG2t",
  "median_nRegGenesOVERnGenes_signif","median_nRegGenesOVERnGenes_notSignif","median_nRegGenesOVERnGenes_permutCorr","median_nRegGenesOVERnGenes_permutG2t",
  "median_nTFsOVERnGenes_signif","median_nTFsOVERnGenes_notSignif","median_nTFsOVERnGenes_permutCorr","median_nTFsOVERnGenes_permutG2t",
  "median_nTFs_signif","median_nTFs_notSignif","median_nTFs_permutCorr","median_nTFs_permutG2t"
)

stopifnot(plot_cols %in% colnames(plot_dt))

# load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/permutCorr_nRegFeat_dt.Rdata")
outFile <- file.path(outFolder, paste0("cmp_obs_permut_nRegFeat_boxplot_allDS.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(mar=par()$mar+c(9,0,0,0))
boxplot(plot_dt[,plot_cols],
        las=2, 
        main=paste0("all ds (n=", length(unique(file.path(plot_dt$hicds, plot_dt$exprds))),") - ", dsIn),  
        cex.main = plotCex , cex.lab = plotCex,
        cex.axis=0.8)
mtext(side=3, text = paste0(dsIn))
cat(paste0("... written: ", outFile, "\n"))

# load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/permutCorr_nRegFeat_dt.Rdata")

keepCols <- c("median_nTFs_signif", "median_nTFs_notSignif", "median_nTFs_permutCorr","median_nTFs_permutG2t",  
              "median_nGenes_signif", "median_nGenes_notSignif", "median_nGenes_permutCorr","median_nGenes_permutG2t",  
              "median_nTFsOVERnGenes_signif", "median_nTFsOVERnGenes_notSignif", "median_nTFsOVERnGenes_permutCorr", "median_nTFsOVERnGenes_permutG2t")

stopifnot(keepCols %in% colnames(plot_dt))

outFile <- file.path(outFolder, paste0("cmp_obs_permut_nRegFeat_boxplot_allDS_keepCols.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(mar=par()$mar+c(9,0,0,0))
boxplot(plot_dt[, keepCols], las=2, 
        main=paste0("all ds (n=", length(unique(file.path(permutCorr_nRegFeat_dt$hicds, permutCorr_nRegFeat_dt$exprds))),") - ", dsIn),  
        cex.main = plotCex , cex.lab = plotCex,
        cex.axis=0.8)
mtext(side=3, text = paste0(dsIn))
cat(paste0("... written: ", outFile, "\n"))



#####################################################################
cat("*** DONE\n")
