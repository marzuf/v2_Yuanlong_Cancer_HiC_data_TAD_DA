
# Rscript tad_ctcf_boundary.R

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

buildTable <- TRUE

runFolder <- "." 
outFolder <- file.path("TAD_CTCF_BOUNDARY")
dir.create(outFolder, recursive = TRUE)


inFile <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
final_dt <- get(load(inFile))
final_dt$chr <- gsub("(chr.+)_TAD.+", "\\1", final_dt$region)
stopifnot(final_dt$chr %in% paste0("chr", 1:22))
final_dt$tad_midpos <- (final_dt$start + final_dt$end)/2 
final_dt$dataset <- file.path(final_dt$hicds, final_dt$exprds)

init_ctcf_dt <- read_excel("13059_2020_2108_MOESM2_ESM.xlsx", sheet="CTCFs")
init_ctcf_dt <- as.data.frame(init_ctcf_dt)
init_ctcf_dt <- init_ctcf_dt[, 1:7]
init_ctcf_dt$chr <- as.character(init_ctcf_dt$chr)
init_ctcf_dt$midpos <- (init_ctcf_dt$start+init_ctcf_dt$end)/2
forward_dt <- init_ctcf_dt[init_ctcf_dt$orientation == ">",]
reverse_dt <- init_ctcf_dt[init_ctcf_dt$orientation == "<",]
stopifnot(nrow(forward_dt) + nrow(reverse_dt) == nrow(init_ctcf_dt))

all_chrs <- unique(final_dt$chr)
chr="chr1"

if(buildTable) {
  
  ### need to split by dataset and sort
  ctcf_bd_peaks_dt <- foreach(dataset = unique(final_dt$dataset), .combine='rbind') %do% {
    
    hicds <- dirname(dataset)
    tad_dt <- read.delim(file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt"), stringsAsFactors = FALSE, 
                         header=F, col.names = c("chromo", "region", "start", "end"))
    # tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
    stopifnot(is.numeric(tad_dt$start))
    stopifnot(is.numeric(tad_dt$end))
    tad_dt$tad_midpos <- (tad_dt$start + tad_dt$end)/2
    
    dataset_ctcf_bd_peaks_dt <- foreach(chr = all_chrs, .combine='rbind') %do% {
      
      chr_tad_dt <- tad_dt[tad_dt$chr == chr,]
      chr_tad_dt <- chr_tad_dt[order(chr_tad_dt$start, chr_tad_dt$end),]
      stopifnot(nrow(chr_tad_dt) > 0)
      stopifnot(nrow(chr_tad_dt) > 3)
      
      all_chr_dt <- foreach(i = 2:(nrow(chr_tad_dt)-1), .combine='rbind') %dopar% {
        
        
        stopifnot(i-1 > 0)
        stopifnot(i+1 <= nrow(chr_tad_dt))
        
        hicds <- dirname(dataset)
        exprds <- basename(dataset)
        region <- chr_tad_dt$region[i]
      
        tad_start <- chr_tad_dt$start[i]
        tad_end <- chr_tad_dt$end[i]
      
        adjPvalComb <- final_dt$adjPvalComb[final_dt$dataset == dataset &
                                              final_dt$region == region]
        adjPvalComb <- ifelse(length(adjPvalComb) == 0, NA, adjPvalComb)
          
        stopifnot(grepl(chr, chr_tad_dt$region))
        
        curr_midpos <- chr_tad_dt$tad_midpos[i]
        prev_midpos <- chr_tad_dt$tad_midpos[i-1]
        next_midpos <- chr_tad_dt$tad_midpos[i+1]
        stopifnot(!is.na(prev_midpos))
        stopifnot(!is.na(next_midpos))
        # retrieve highest peaks 
        
        subLeft_ctcf_dt <- init_ctcf_dt[init_ctcf_dt$chr == chr &
                                          init_ctcf_dt$midpos >= prev_midpos &
                                          init_ctcf_dt$midpos <= curr_midpos,]
        
        subRight_ctcf_dt <- init_ctcf_dt[init_ctcf_dt$chr == chr &
                                           init_ctcf_dt$midpos >= curr_midpos &
                                           init_ctcf_dt$midpos <= next_midpos,]
        
        
        # find highest peaks - LEFT
        leftHighestMotifScore <- max(subLeft_ctcf_dt$MotifScore)
        leftHighestChipSeqScore <- max(subLeft_ctcf_dt$ChipSeqScore)
        leftHighestMotifScoreForward <- max(subLeft_ctcf_dt$MotifScore[subLeft_ctcf_dt$orientation == ">"])
        leftHighestChipSeqScoreForward <- max(subLeft_ctcf_dt$ChipSeqScore[subLeft_ctcf_dt$orientation == ">"])
        
        # find highest peaks - RIGHT
        rightHighestMotifScore <- max(subRight_ctcf_dt$MotifScore)
        rightHighestChipSeqScore <- max(subRight_ctcf_dt$ChipSeqScore)
        rightHighestMotifScoreForward <- max(subRight_ctcf_dt$MotifScore[subRight_ctcf_dt$orientation == "<"])
        rightHighestChipSeqScoreForward <- max(subRight_ctcf_dt$ChipSeqScore[subRight_ctcf_dt$orientation == "<"])
        
        
        
        data.frame(
          hicds=hicds,
          exprds=exprds,
          region=region,
          adjPvalComb=adjPvalComb,
          
          leftHighestMotifScore = ifelse(is.infinite(leftHighestMotifScore), NA,leftHighestMotifScore ),
          leftHighestChipSeqScore =ifelse(is.infinite(leftHighestChipSeqScore), NA,leftHighestChipSeqScore ),
          leftHighestMotifScoreForward = ifelse(is.infinite(leftHighestMotifScoreForward), NA, leftHighestMotifScoreForward),
          leftHighestChipSeqScoreForward = ifelse(is.infinite(leftHighestChipSeqScoreForward), NA, leftHighestChipSeqScoreForward),
          
          # find highest peaks - RIGHT
          rightHighestMotifScore = ifelse(is.infinite(rightHighestMotifScore), NA,rightHighestMotifScore ),
          rightHighestChipSeqScore = ifelse(is.infinite(rightHighestChipSeqScore), NA,rightHighestChipSeqScore ),
          rightHighestMotifScoreForward= ifelse(is.infinite(rightHighestMotifScoreForward), NA,rightHighestMotifScoreForward ),
          rightHighestChipSeqScoreForward = ifelse(is.infinite(rightHighestChipSeqScoreForward), NA, rightHighestChipSeqScoreForward),
          
          
          stringsAsFactors = FALSE
        )
      }
      all_chr_dt
    }
    dataset_ctcf_bd_peaks_dt
  }
  
  outFile <- file.path(outFolder, "ctcf_bd_peaks_dt.Rdata")
  save(ctcf_bd_peaks_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,"\n"))
  
}else {
  outFile <- file.path(outFolder, "ctcf_bd_peaks_dt.Rdata")
  ctcf_bd_peaks_dt <- get(load(outFile))
}

bound_dt <- ctcf_bd_peaks_dt[grepl("BOUND", ctcf_bd_peaks_dt$region),]
stopifnot(is.na(bound_dt$adjPvalComb))

pip_dt <- ctcf_bd_peaks_dt
pip_dt$region_id <- file.path(pip_dt$hicds, pip_dt$exprds, pip_dt$region)
pip_dt <- pip_dt[pip_dt$region_id %in% file.path(final_dt$hicds, final_dt$exprds, final_dt$region),]
stopifnot(!is.na(pip_dt$adjPvalComb))

pip_dt$signif_lab <- ifelse(pip_dt$adjPvalComb <= pthresh, "signif.", "not signif.")

pip_dt$adjPvalComb_log10 <- -log10(pip_dt$adjPvalComb)

plot_cols <- colnames(pip_dt)[grepl("Score", colnames(pip_dt))]

plot_col="leftHighestMotifScore"

for(plot_col in plot_cols) {
  # ctcf_bd_dist_dt[paste0(curr_col, "_log10")] <- log10(ctcf_bd_dist_dt[paste0(curr_col, "")])
  # 
  p <- ggdensity(pip_dt,
                 x = paste0(plot_col),
                 y = "..density..",
                 # combine = TRUE,                  # Combine the 3 plots
                 xlab = paste0(plot_col),
                 # add = "median",                  # Add median line.
                 rug = FALSE,                      # Add marginal rug
                 color = "signif_lab",
                 fill = "signif_lab",
                 palette = "d3"
  )
  outFile <- file.path(outFolder, paste0(plot_col, "_bySignif_densityplot.", plotTypeGG))
  ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  yvar <- paste0(plot_col)
  xvar <- paste0("adjPvalComb_log10")
  
  outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_densplot_withCorr(xvar, yvar, pip_dt)
  mtext(side=3, text = paste0("nSignif=", sum(pip_dt$signif_lab=="signif."), 
                              "; nNotSignif=", sum(pip_dt$signif_lab=="not signif."), 
                              "; nDS=", length(unique(pip_dt$hicds, pip_dt$exprds))))
  title(main = paste0(yvar, " vs. ", xvar))
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
  
  
  # p <- ggdensity(pip_dt,
  #                x = paste0(curr_col, "_log10"),
  #                y = "..density..",
  #                # combine = TRUE,                  # Combine the 3 plots
  #                xlab = paste0(curr_col, "_log10"),
  #                # add = "median",                  # Add median line.
  #                rug = FALSE,                      # Add marginal rug
  #                color = "signif_lab",
  #                fill = "signif_lab",
  #                palette = "d3"
  # ) 
  # outFile <- file.path(outFolder, paste0(curr_col, "_bySignif_densityplot_log10.", plotTypeGG))
  # ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
  # cat(paste0("... written: ", outFile,  "\n"))
  
  
}




# between adjacent midpos -> highest forward left and highest reverse right

# cluster ctcf -> nbr by TADs


