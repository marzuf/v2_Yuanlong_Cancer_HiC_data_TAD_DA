
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
# Rscript meanCorr_pval_deltaRank.R

plotType <- "svg"
myHeight <- myWidth <- 7
plotType <- "png"
myHeight <- myWidth <- 400

outFolder <- "MEANCORR_PVAL_DELTARANK"
dir.create(outFolder, recursive = TRUE)

all_hicds_init <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds <- all_hicds_init[grepl("RANDOMMIDPOS", all_hicds_init)]
# => will need to add RANDOMMIDPOSDISC ONCE I HAVE PVALS
all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

hicds = all_hicds[1]
all_ds_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    meanCorr <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata")))
    
    empPval <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
    adj_empPval <- p.adjust(empPval, method="BH")
    
    tad_withRank_dt <- get(load(file.path(gsub("RANDOMMIDPOSDISC", "RANDOMMIDPOS", hicds), "ASSIGNED_REGIONS_DELTARANK", "all_assigned_regions_withDeltaRank.Rdata")))
    stopifnot(names(meanCorr) %in% tad_withRank_dt$region)
    stopifnot(setequal(names(meanCorr), names(adj_empPval)))
    
    tadRank <- setNames(tad_withRank_dt$abs_delta_rank, tad_withRank_dt$region)
    
    out_dt <- data.frame(
      hicds = hicds,
      exprds = exprds,
      region = names(meanCorr),
      meanCorr = as.numeric(meanCorr[names(meanCorr)]),
      adjPval = as.numeric(adj_empPval[names(meanCorr)]),
      absDeltaRank = as.numeric(tadRank[names(meanCorr)]),
      stringsAsFactors = FALSE
      )
    stopifnot(!is.na(out_dt))
    out_dt
  }
  exprds_dt
}

outFile <- file.path(outFolder, "all_ds_dt.Rdata")
save(all_ds_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

all_ds_dt$hicds_type <-  gsub(".+_(RANDOM.+)_40kb", "\\1", all_ds_dt$hicds)

ex_hicds <- "ENCSR489OCU_NCI-H460"
ex_exprds <- "TCGAluad_norm_luad"

h_t = unique(all_ds_dt$hicds_type)[1]
for(h_t in unique(all_ds_dt$hicds_type)) {
  
  plot_dt <- all_ds_dt[all_ds_dt$hicds_type == h_t,]
  plot_dt_sub <- plot_dt[grepl(ex_hicds, plot_dt$hicds) & plot_dt$exprds == ex_exprds,]
  
  nDS <- length(unique(file.path(plot_dt$hicds, plot_dt$exprds)))
  
  source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
  
  plotCex <- 1.4
  
  my_x <- plot_dt$absDeltaRank
  my_x_sub <- plot_dt_sub$absDeltaRank
  
  yvar = c("meanCorr")[1]
  for(yvar in c("meanCorr", "adjPval")){
    my_y <-  plot_dt[,paste0(yvar)]
    my_y_sub <-  plot_dt_sub[,paste0(yvar)]
    
    outFile <- file.path(outFolder, paste0(yvar, "_vs_absDeltaRank_", h_t, "_allDS.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    par(bty="L")
    densplot(
      x = my_x,
      y= my_y,
      cex=0.7,
      main = paste0(h_t, " - ", yvar, " vs. absDeltaRank"),
      xlab = "absDeltaRank",
      ylab = yvar,
      cex.main=plotCex,
      cex.lab = plotCex,
      cex.axis = plotCex
    )
    mtext(side=3, text = paste0("all DS - n=", nDS, " - #TADs=", nrow(plot_dt)))
    addCorr(x = my_x, y = my_y, bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    outFile <- file.path(outFolder, paste0(yvar, "_vs_absDeltaRank_", h_t, "_", ex_hicds, "_", ex_exprds, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    par(bty="L")
    densplot(
      x = my_x_sub,
      y= my_y_sub,
      cex=0.7,
      main = paste0(h_t, " - ", yvar, " vs. absDeltaRank"),
      xlab = "absDeltaRank",
      ylab = yvar,
      cex.main=plotCex,
      cex.lab = plotCex,
      cex.axis = plotCex
    )
    mtext(side=3, text = paste0(ex_hicds, " - ", ex_exprds, " - #TADs=", nrow(plot_dt)))
    addCorr(x = my_x_sub, y = my_y_sub, bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
  }
  
  
  
  
  
  
}





