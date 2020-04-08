
# Rscript investigate_fcc_auc.R

require(foreach)
require(doMC)
registerDoMC(40)

outFolder <- "INVESTIGATE_FCC_AUC"
dir.create(outFolder, recursive = TRUE)

plotType <- "svg"
myHeight <- myWidth <- 7

buildTable <- TRUE

hicds <- "LG1_40kb"
exprds <- "TCGAluad_norm_luad"

rd_type <- "RANDOMMIDPOS"

pipFolder <- "PIPELINE/OUTPUT_FOLDER/"

all_hicds <- list.files(file.path(pipFolder))
all_hicds <- all_hicds[!grepl("_RANDOM", all_hicds)]
all_hicds <- all_hicds[!grepl("_PERMUT", all_hicds)]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

if(buildTable) {
  
  hicds <- all_hicds[18]
  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    
    exprds = all_exprds[[paste0(hicds)]][1]
    expr_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      
      g2t_dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), stringsAsFactors = FALSE, header=FALSE,
                           col.names=c("entrezID", "chromo", "start", "end", "region"))
      g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
      
      
      rd_g2t_dt <- read.delim(file.path(gsub("_40kb", paste0("_", rd_type, "_40kb"), hicds), "genes2tad", "all_genes_positions.txt"), stringsAsFactors = FALSE,
                              col.names=c("entrezID", "chromo", "start", "end", "region"), header=FALSE)
      rd_g2t_dt$entrezID <- as.character(rd_g2t_dt$entrezID)
      
      
      
      geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
      
      rd_geneList <- get(load(file.path(pipFolder, gsub("_40kb", paste0("_", rd_type, "_40kb"), hicds), exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
      
      stopifnot(geneList %in% g2t_dt$entrezID)
      stopifnot(rd_geneList %in% g2t_dt$entrezID)
      
      g2t_obs_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]      
      g2t_obs <- setNames(as.numeric(table(g2t_obs_dt$region)), names(table(g2t_obs_dt$region)))
      stopifnot(grepl("_TAD", names(g2t_obs)))
      
      g2t_rd_dt <- rd_g2t_dt[rd_g2t_dt$entrezID %in% rd_geneList,]      
      g2t_rd <- setNames(as.numeric(table(g2t_rd_dt$region)), names(table(g2t_rd_dt$region)))
      stopifnot(grepl("_TAD", names(g2t_rd)))
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        nbrTADs_obs = length(g2t_obs),
        meanG2t_obs = mean(g2t_obs),
        nbrTADs_rd = length(g2t_rd),
        meanG2t_rd = mean(g2t_rd),
        stringsAsFactors = FALSE
      )
      
      
    }
    expr_dt
  }
  
  outFile <- file.path(outFolder, "all_dt.Rdata")
  save(all_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {
  outFile <- file.path(outFolder, "all_dt.Rdata")
  all_dt <- get(load(outFile))
  
}

feature_dt <- all_dt

rd_fcc_auc_dt <- get(load(paste0("../FIGURES_V5_YUANLONG/RANDOM_FCC_AUC_RATIO_RANDOMMIDPOS_V2/all_auc_ratio_dt_", rd_type, ".Rdata")))
permg2t_fcc_auc_dt <- get(load(paste0("../FIGURES_V2_YUANLONG/BARPLOT_FCC_AUC_RATIO/all_dt.Rdata")))


plot_dt <- merge(permg2t_fcc_auc_dt[,c("hicds", "exprds", "fcc_auc")], by=c("hicds", "exprds"), all=TRUE,
                 merge(feature_dt, rd_fcc_auc_dt, by=c("hicds", "exprds"), all=TRUE))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")




plot_dt$diff_auc <- plot_dt$fcc_auc - plot_dt$rd_fcc_auc
plot_dt$diff_nTADs <- plot_dt$nbrTADs_obs - plot_dt$nbrTADs_rd
plot_dt$diff_g2t <- plot_dt$meanG2t_obs - plot_dt$meanG2t_rd

plot_dt$dotcols <- all_cols[all_cmps[plot_dt$exprds]]

#########################################
yvar <- "diff_auc"
my_y <- plot_dt[,paste0(yvar)]
ylab <- "Delta FCC AUC ratio"


for(xvar in c("diff_nTADs", "diff_g2t")) {
  
  if(xvar=="diff_nTADs") xlab <- "Delta # of TADs"
  if(xvar=="diff_g2t") xlab <- "Delta mean genes/TADs"
  
  my_x <- plot_dt[,paste0(xvar)]
  
  outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "_", rd_type, "_dotplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x  = my_x,
    y  = my_y,
    main = paste0(ylab, " ~ ", xlab),
    xlab=xlab,
    ylab=ylab,
    cex.main=plotCex,
    cex.lab=plotCex,
    cex.axis=plotCex,
    col=plot_dt$dotcols,
    pch=16,
    cex=0.7
  )
  addCorr(x=my_x, y=my_y, bty="n")
  mtext(side=3, text = paste0("PERMG2T vs. ", rd_type))
  #  TCGAcoad_msi_mss
  # xdt[xdt$PERMG2T_FCC_AUC >= 1.5 & xdt$RANDOMMIDPOSSTRICT_FCC_AUC <= 1.2,]
  # ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1
  
  
  hicds1 <- "ENCSR504OTV_transverse_colon_40kb"
  exprds1 <- "TCGAcoad_msi_mss"
  
  plot_dt$dataset <- paste0(plot_dt$hicds, "\n", plot_dt$exprds)
  
  text(x = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,paste0(xvar)],
         y = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,paste0(yvar)],
         labels = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,"dataset"])
  
  
  hicds2 <- "ENCSR312KHQ_SK-MEL-5_40kb"
  exprds2 <- "TCGAskcm_wt_mutCTNNB1"
  
  text(x = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,paste0(xvar)],
       y = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,paste0(yvar)],
       labels = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,"dataset"])
  
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}

#########################################

yvar <- "nbrTADs_obs"
my_y <- plot_dt[,paste0(yvar)]
ylab <- "# TADs (PERMG2T=OBS.)"

xvar <- "nbrTADs_rd"
my_x <- plot_dt[,paste0(xvar)]
xlab <- paste0("# TADs (", rd_type, ")")

outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "_", rd_type, "_dotplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x  = my_x,
  y  = my_y,
  main = paste0(ylab, " ~ ", xlab),
  xlab=xlab,
  ylab=ylab,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex,
  col=plot_dt$dotcols,
  pch=16,
  cex=0.7
)
curve(x*1, lty=2, add=T, col="grey")
addCorr(x=my_x, y=my_y, bty="n")
mtext(side=3, text = paste0("PERMG2T vs. ", rd_type))
#  TCGAcoad_msi_mss
# xdt[xdt$PERMG2T_FCC_AUC >= 1.5 & xdt$RANDOMMIDPOSSTRICT_FCC_AUC <= 1.2,]
# ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1


hicds1 <- "ENCSR504OTV_transverse_colon_40kb"
exprds1 <- "TCGAcoad_msi_mss"

plot_dt$dataset <- paste0(plot_dt$hicds, "\n", plot_dt$exprds)

text(x = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,paste0(xvar)],
     y = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,paste0(yvar)],
     labels = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,"dataset"])


hicds2 <- "ENCSR312KHQ_SK-MEL-5_40kb"
exprds2 <- "TCGAskcm_wt_mutCTNNB1"

text(x = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,paste0(xvar)],
     y = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,paste0(yvar)],
     labels = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,"dataset"])


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#########################################


yvar <- "meanG2t_obs"
my_y <- plot_dt[,paste0(yvar)]
ylab <- "avg. genes by TAD (PERMG2T=OBS.)"

xvar <- "meanG2t_rd"
my_x <- plot_dt[,paste0(xvar)]
xlab <- paste0("avg. genes by TAD (", rd_type, ")")

outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "_", rd_type, "_dotplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x  = my_x,
  y  = my_y,
  main = paste0(ylab, " ~ ", xlab),
  xlab=xlab,
  ylab=ylab,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex,
  col=plot_dt$dotcols,
  pch=16,
  cex=0.7
)
curve(x*1, lty=2, add=T, col="grey")
addCorr(x=my_x, y=my_y, bty="n")
mtext(side=3, text = paste0("PERMG2T vs. ", rd_type))
#  TCGAcoad_msi_mss
# xdt[xdt$PERMG2T_FCC_AUC >= 1.5 & xdt$RANDOMMIDPOSSTRICT_FCC_AUC <= 1.2,]
# ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1


hicds1 <- "ENCSR504OTV_transverse_colon_40kb"
exprds1 <- "TCGAcoad_msi_mss"

plot_dt$dataset <- paste0(plot_dt$hicds, "\n", plot_dt$exprds)

text(x = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,paste0(xvar)],
     y = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,paste0(yvar)],
     labels = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,"dataset"])


hicds2 <- "ENCSR312KHQ_SK-MEL-5_40kb"
exprds2 <- "TCGAskcm_wt_mutCTNNB1"

text(x = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,paste0(xvar)],
     y = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,paste0(yvar)],
     labels = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,"dataset"])


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#########################################

all_cols <- colnames(plot_dt)
all_cols <- all_cols[!all_cols %in% c("hicds", "exprds", "dataset", "dotcols")]

for(i in 1:(length(all_cols)-1)) {
  for(j in (i+1):length(all_cols)){
    
    xvar <- all_cols[i]
    yvar <- all_cols[j]
    
    my_y <- plot_dt[,paste0(yvar)]
    ylab <- yvar
    
    my_x <- plot_dt[,paste0(xvar)]
    xlab <- xvar
    
    outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "_", rd_type, "_raw_dotplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(
      x  = my_x,
      y  = my_y,
      main = paste0(ylab, " ~ ", xlab),
      xlab=xlab,
      ylab=ylab,
      cex.main=plotCex,
      cex.lab=plotCex,
      cex.axis=plotCex,
      col=plot_dt$dotcols,
      pch=16,
      cex=0.7
    )
    addCorr(x=my_x, y=my_y, bty="n")
    # mtext(side=3, text = paste0("PERMG2T vs. ", rd_type))
    #  TCGAcoad_msi_mss
    # xdt[xdt$PERMG2T_FCC_AUC >= 1.5 & xdt$RANDOMMIDPOSSTRICT_FCC_AUC <= 1.2,]
    # ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1
    
    
    hicds1 <- "ENCSR504OTV_transverse_colon_40kb"
    exprds1 <- "TCGAcoad_msi_mss"
    
    plot_dt$dataset <- paste0(plot_dt$hicds, "\n", plot_dt$exprds)
    
    text(x = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,paste0(xvar)],
         y = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,paste0(yvar)],
         labels = plot_dt[plot_dt$hicds == hicds1 & plot_dt$exprds == exprds1,"dataset"])
    
    
    hicds2 <- "ENCSR312KHQ_SK-MEL-5_40kb"
    exprds2 <- "TCGAskcm_wt_mutCTNNB1"
    
    text(x = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,paste0(xvar)],
         y = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,paste0(yvar)],
         labels = plot_dt[plot_dt$hicds == hicds2 & plot_dt$exprds == exprds2,"dataset"])
    
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
    
  }
}




