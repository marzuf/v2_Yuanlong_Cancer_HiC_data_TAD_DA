startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_allRandom.R c3.tft
# Rscript TFs_by_TADs_signifTADs_v2_cmp_obs_permut_allRandom.R chea3_lung

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
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))
}

stopifnot(dsIn %in% c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg", "chea3_all", "chea3_lung"))

outFolder <- file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_CMP_OBS_PERMUT_ALLRANDOM_", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

obs_data <- get(load(file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_", toupper(dsIn)), "nRegFeat_dt.Rdata")))
permutCorr_data <- get(load(file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_PERMUTCORR_", toupper(dsIn)), "permutCorr_nRegFeat_dt.Rdata")))
permutG2t_data <- get(load(file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_PERMUTG2T1000_", toupper(dsIn)), "permutG2t_nRegFeat_dt.Rdata")))

all_random_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T" )
rd_patt = all_random_patterns[1]

plot_dt <- merge(obs_data, merge(permutG2t_data, permutCorr_data, by =c("hicds", "exprds")), by =c("hicds", "exprds"))

for(rd_patt in all_random_patterns) {
  rd_dt <- permutCorr_data <- get(load(file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_", rd_patt, "_V2_", toupper(dsIn)), "nRegFeat_dt.Rdata")))
  colnames(rd_dt)[!colnames(rd_dt) %in% c("hicds", "exprds")] <- paste0(colnames(rd_dt)[!colnames(rd_dt) %in% c("hicds", "exprds")], "_", rd_patt)
  rd_dt$hicds <- gsub(paste0("_", rd_patt, "_40kb"), "_40kb", rd_dt$hicds)
  
  plot_dt <- merge(plot_dt, rd_dt, by =c("hicds", "exprds"))
}
stopifnot(nrow(plot_dt) == nrow(obs_data))


all_plot_cols <- list(
  c("mean_nGenes_signif"  ,"mean_nGenes_notSignif", "mean_nGenes_permutCorr","mean_nGenes_permutG2t",paste0("mean_nGenes_", all_random_patterns)),
    c("mean_nRegGenes_signif","mean_nRegGenes_notSignif","mean_nRegGenes_permutCorr","mean_nRegGenes_permutG2t",paste0("mean_nRegGenes_", all_random_patterns)),
      c("mean_nRegGenesOVERnGenes_signif","mean_nRegGenesOVERnGenes_notSignif","mean_nRegGenesOVERnGenes_permutCorr","mean_nRegGenesOVERnGenes_permutG2t", paste0("mean_nRegGenesOVERnGenes_", all_random_patterns)),
        c("mean_nTFsOVERnGenes_signif","mean_nTFsOVERnGenes_notSignif","mean_nTFsOVERnGenes_permutCorr","mean_nTFsOVERnGenes_permutG2t",paste0("mean_nTFsOVERnGenes_", all_random_patterns)),
          c("mean_nTFs_signif","mean_nTFs_notSignif","mean_nTFs_permutCorr","mean_nTFs_permutG2t", paste0("mean_nTFs_", all_random_patterns)))


stopifnot(unlist(all_plot_cols) %in% colnames(plot_dt))

for(plot_cols in all_plot_cols) {

  plotTit <- unique(gsub("^(.+?_.+?)_.+", "\\1", unlist(plot_cols)))
  stopifnot(length(plotTit) == 1)
  
  # load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/permutCorr_nRegFeat_dt.Rdata")
  outFile <- file.path(outFolder, paste0("cmp_obs_permut_", plotTit, "_boxplot_allDS_allRandom.", plotType))  
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(mar=par()$mar+c(9,0,0,0))
  boxplot(plot_dt[,plot_cols],
          las=2, 
          main=paste0("all ds (n=", length(unique(file.path(plot_dt$hicds, plot_dt$exprds))),") - ", dsIn),  
          cex.main = plotCex , cex.lab = plotCex,
          cex.axis=0.8)
  mtext(side=3, text = paste0(dsIn))
  cat(paste0("... written: ", outFile, "\n"))
  
    
}
# 
# # load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/permutCorr_nRegFeat_dt.Rdata")
# 
# keepCols <- c("mean_nTFs_signif", "mean_nTFs_notSignif", "mean_nTFs_permutCorr","mean_nTFs_permutG2t",  
#               "mean_nGenes_signif", "mean_nGenes_notSignif", "mean_nGenes_permutCorr","mean_nGenes_permutG2t",  
#               "mean_nTFsOVERnGenes_signif", "mean_nTFsOVERnGenes_notSignif", "mean_nTFsOVERnGenes_permutCorr", "mean_nTFsOVERnGenes_permutG2t")
# 
# stopifnot(keepCols %in% colnames(plot_dt))
# 
# outFile <- file.path(outFolder, paste0("cmp_obs_permut_nRegFeat_boxplot_allDS_keepCols.", plotType))  
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# par(mar=par()$mar+c(9,0,0,0))
# boxplot(plot_dt[, keepCols], las=2, 
#         main=paste0("all ds (n=", length(unique(file.path(permutCorr_nRegFeat_dt$hicds, permutCorr_nRegFeat_dt$exprds))),") - ", dsIn),  
#         cex.main = plotCex , cex.lab = plotCex,
#         cex.axis=0.8)
# mtext(side=3, text = paste0(dsIn))
# cat(paste0("... written: ", outFile, "\n"))
# 
# 

#####################################################################
cat("*** DONE\n")
