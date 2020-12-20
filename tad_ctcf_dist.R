
# Rscript tad_ctcf_dist.R

library("readxl")
library(doMC)
library(foreach)
library(stringr)

require(ggpubr)
require(ggsci)


pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_obs_exprds <- sapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))

registerDoMC(40)
# runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878" #PIPELINE/OUTPUT_FOLDER/GM12878_40kb/TCGAluad_norm_luad/11sameNbr_runEmpPvalCombined/"
# hicds <- "GM12878_40kb"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("ctcf_da_utils.R")

plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.2

plotTypeGG <- "svg"
ggHeight <- 6
ggWidth <- 5

pthresh <- 0.01

fontFamily <- "Hershey"

buildTable <- F

runFolder <- "." 
outFolder <- file.path("TAD_CTCF_DIST")
dir.create(outFolder, recursive = TRUE)


inFile <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
final_dt <- get(load(inFile))
final_dt$chr <- gsub("(chr.+)_TAD.+","\\1", final_dt$region)
stopifnot(final_dt$chr %in% paste0("chr", 1:22))

init_ctcf_dt <- read_excel("13059_2020_2108_MOESM2_ESM.xlsx", sheet="CTCFs")
init_ctcf_dt <- as.data.frame(init_ctcf_dt)
init_ctcf_dt <- init_ctcf_dt[, 1:7]
init_ctcf_dt$chr <- as.character(init_ctcf_dt$chr)
init_ctcf_dt$midpos <- (init_ctcf_dt$start+init_ctcf_dt$end)/2
forward_dt <- init_ctcf_dt[init_ctcf_dt$orientation == ">",]
reverse_dt <- init_ctcf_dt[init_ctcf_dt$orientation == "<",]
stopifnot(nrow(forward_dt) + nrow(reverse_dt) == nrow(init_ctcf_dt))

if(buildTable) {
  ctcf_bd_dist_dt <- foreach(i = 1:nrow(final_dt), .combine='rbind') %dopar% {
    
    hicds <- final_dt$hicds[i]
    exprds <- final_dt$exprds[i]
    region <- final_dt$region[i]
    adjPvalComb <- final_dt$adjPvalComb[i]
    tad_start <- final_dt$start[i]
    tad_end <- final_dt$end[i]
    chr <- final_dt$chr[i]
    stopifnot(chr %in% init_ctcf_dt$chr)
    
    # find closest 
    i_tadStart_closestAll <- which.min(abs(tad_start-init_ctcf_dt$midpos[init_ctcf_dt$chr == chr]))
    i_tadEnd_closestAll <- which.min(abs(tad_end-init_ctcf_dt$midpos[init_ctcf_dt$chr == chr]))
    i_tadStart_closestForward <- which.min(abs(tad_start-forward_dt$midpos[forward_dt$chr == chr]))
    i_tadEnd_closestReverse <- which.min(abs(tad_end-reverse_dt$midpos[reverse_dt$chr == chr]))

    tadStart_closestAll_dist <- min(abs(tad_start-init_ctcf_dt$midpos[init_ctcf_dt$chr == chr]))
    tadEnd_closestAll_dist <- min(abs(tad_end-init_ctcf_dt$midpos[init_ctcf_dt$chr == chr]))
    tadStart_closestForward_dist <- min(abs(tad_start-forward_dt$midpos[forward_dt$chr == chr]))
    tadEnd_closestReverse_dist <- min(abs(tad_end-reverse_dt$midpos[reverse_dt$chr == chr]))
    
    out_dt <- data.frame(
      hicds=hicds,
      exprds=exprds,
      region=region,
      adjPvalComb=adjPvalComb,
      
      tadStart_closestAll_dist=tadStart_closestAll_dist,
      tadStart_closestAll_MotifScore=init_ctcf_dt$MotifScore[i_tadStart_closestAll],
      tadStart_closestAll_ChipSeqScore=init_ctcf_dt$ChipSeqScore[i_tadStart_closestAll],
      
      tadEnd_closestAll_dist=tadEnd_closestAll_dist,
      tadEnd_closestAll_MotifScore=init_ctcf_dt$MotifScore[i_tadEnd_closestAll],
      tadEnd_closestAll_ChipSeqScore=init_ctcf_dt$ChipSeqScore[i_tadEnd_closestAll],
      
      tadStart_closestForward_dist=tadStart_closestForward_dist,
      tadStart_closestForward_MotifScore=forward_dt$MotifScore[i_tadStart_closestForward],
      tadStart_closestForward_ChipSeqScore=forward_dt$ChipSeqScore[i_tadStart_closestForward],
      
      tadEnd_closestReverse_dist=tadEnd_closestReverse_dist,
      tadEnd_closestReverse_MotifScore=reverse_dt$MotifScore[i_tadEnd_closestReverse],
      tadEnd_closestReverse_ChipSeqScore=reverse_dt$ChipSeqScore[i_tadEnd_closestReverse],
      
      stringsAsFactors = FALSE
    )
    stopifnot(!is.na(out_dt))
    out_dt
  }
  
  outFile <- file.path(outFolder, "ctcf_bd_dist_dt.Rdata")
  save(ctcf_bd_dist_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,"\n"))

}else {
  outFile <- file.path(outFolder, "ctcf_bd_dist_dt.Rdata")
  ctcf_bd_dist_dt <- get(load(outFile))
}

ctcf_bd_dist_dt$signif_lab <- ifelse(ctcf_bd_dist_dt$adjPvalComb <= pthresh, "signif.", "not signif.")
ctcf_bd_dist_dt$adjPvalComb_log10 <- -log10(ctcf_bd_dist_dt$adjPvalComb)

all_cols <- c("tadStart_closestAll", "tadEnd_closestAll", "tadStart_closestForward", "tadEnd_closestReverse")
all_metrics <- c("_dist", "_MotifScore", "_ChipSeqScore")
curr_col="tadStart_closestAll"

for(metrics in all_metrics)  {
  
  for(init_curr_col in all_cols) {
    
    curr_col <- paste0(init_curr_col, metrics)
    
    if(metrics == "_dist") {
      ctcf_bd_dist_dt[paste0(curr_col, "_log10")] <- log10(ctcf_bd_dist_dt[paste0(curr_col, "")])
      p <- ggdensity(ctcf_bd_dist_dt,
                     x = paste0(curr_col, "_log10"),
                     y = "..density..",
                     # combine = TRUE,                  # Combine the 3 plots
                     xlab = paste0(curr_col, "_log10"),
                     # add = "median",                  # Add median line.
                     rug = FALSE,                      # Add marginal rug
                     color = "signif_lab",
                     fill = "signif_lab",
                     palette = "d3"
      ) 
      outFile <- file.path(outFolder, paste0(curr_col, "_bySignif_densityplot_log10.", plotTypeGG))
      ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
      cat(paste0("... written: ", outFile,  "\n"))
      
      yvar <- paste0(curr_col, "_log10")
      xvar <- paste0("adjPvalComb_log10")
      
      outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "_", "log10_densplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      do_densplot_withCorr(xvar, yvar, ctcf_bd_dist_dt)
      mtext(side=3, text = paste0("nSignif=", sum(ctcf_bd_dist_dt$signif_lab=="signif."), 
                                  "; nNotSignif=", sum(ctcf_bd_dist_dt$signif_lab=="not signif."), 
                                  "; nDS=", length(unique(ctcf_bd_dist_dt$hicds, ctcf_bd_dist_dt$exprds))))
      title(main = paste0(yvar, " vs. ", xvar))
      foo <- dev.off()
      cat(paste0("... written: ", outFile,"\n"))
      
    }
    
    p <- ggdensity(ctcf_bd_dist_dt,
                   x = paste0(curr_col),
                   y = "..density..",
                   # combine = TRUE,                  # Combine the 3 plots
                   xlab = paste0(curr_col),
                   # add = "median",                  # Add median line.
                   rug = FALSE,                      # Add marginal rug
                   color = "signif_lab",
                   fill = "signif_lab",
                   palette = "d3"
    ) 
    outFile <- file.path(outFolder, paste0(curr_col, "_bySignif_densityplot.", plotTypeGG))
    ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
    cat(paste0("... written: ", outFile,  "\n"))
    
    yvar <- paste0(curr_col)
    xvar <- paste0("adjPvalComb_log10")
    
    outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "_densplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    do_densplot_withCorr(xvar, yvar, ctcf_bd_dist_dt)
    mtext(side=3, text = paste0("nSignif=", sum(ctcf_bd_dist_dt$signif_lab=="signif."), 
                                "; nNotSignif=", sum(ctcf_bd_dist_dt$signif_lab=="not signif."), 
                                "; nDS=", length(unique(ctcf_bd_dist_dt$hicds, ctcf_bd_dist_dt$exprds))))
    title(main = paste0(yvar, " vs. ", xvar))
    foo <- dev.off()
    cat(paste0("... written: ", outFile,"\n"))
    
    
  }
  
}




# between adjacent midpos -> highest forward left and highest reverse right

# cluster ctcf -> nbr by TADs
               

