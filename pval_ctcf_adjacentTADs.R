# Rscript pval_ctcf_adjacentTADs.R

require(foreach)
require(doMC)
registerDoMC(40)

plotType <- "png"
myWidth <- myHeight <- 400
plotCex <- 1.2

nPermut <- 1000

outFolder <- file.path("PVAL_CTCF_ADJACENTTADS")
dir.create(outFolder, recursive = TRUE)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

ctcf_bd_dt <- get(load("TAD_CTCF_BOUNDARY/ctcf_bd_peaks_dt.Rdata")) # here i have all of TADs
# > nrow(ctcf_bd_dt)
# [1] 603127
ctcf_bd_dt$regionID <- file.path(ctcf_bd_dt$hicds, ctcf_bd_dt$exprds, ctcf_bd_dt$region)

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata")) # here only TADs kept in poipeline
# [1] 100884
final_dt$regionID <- file.path(final_dt$hicds, final_dt$exprds, final_dt$region)

nDS <- length(unique(file.path(final_dt$hicds, final_dt$exprds)))

# stopifnot(final_dt$regionID %in% ctcf_bd_dt$regionID)  # not true because remove first and last TADs for the boundary analysis
missing_regs <- final_dt$region[!final_dt$regionID %in% ctcf_bd_dt$regionID]
missing_chrs <- gsub("(chr.+)_TAD.+", "\\1", missing_regs)
# but still should not have more than 2 missing per chromo per datasets
stopifnot(table(missing_chrs) <= 2*nDS)

# foreach dataset
# foreach chromo
# order the TADs
# foreach domain i: cmp score_right with diff in pval with score du domain i+1
# or equivalent cmp score_left with diff in pvla with domain i-1

# might be na in the scores if nothing in between ### control during iteration because when switch chromo not true
# tmp_dt <- ctcf_bd_peaks_dt[,c("rightHighestChipSeqScore", "leftHighestChipSeqScore")]
# tmp_dt <- na.omit(tmp_dt)
# tmp1 <- tmp_dt$rightHighestChipSeqScore[1:(nrow(tmp_dt)-1)]
# tmp2 <- tmp_dt$leftHighestChipSeqScore[2:(nrow(tmp_dt))]

merged_dt <- merge(ctcf_bd_dt, final_dt, by=c("hicds", "exprds", "region", "regionID", "adjPvalComb"), all=FALSE)
merged_dt$chr <- gsub("(chr.+)_TAD.+", "\\1", merged_dt$region)
stopifnot(merged_dt$chr %in% paste0("chr", 1:22))
stopifnot(!duplicated(merged_dt$regionID))

merged_dt$tmp_ds <- file.path(merged_dt$hicds, merged_dt$exprds, merged_dt$chr)


hicds="Barutcu_MCF-10A_40kb"
exprds="TCGAbrca_lum_bas"
chromo="chr21"

diff_var = c("meanCorr")

score_var = c("rightHighestChipSeqScore")

for(diff_var in c("meanLogFC", "meanCorr", "adjPvalComb")) {
  
  for(score_var in c("rightHighestMotifScore", "rightHighestChipSeqScore")){
    
    ds= "Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/chr1"
    all_ds_dt <- foreach(ds = unique(merged_dt$tmp_ds), .combine='rbind') %dopar% {
      
      hicds <- dirname(dirname(ds))
      exprds <- basename(dirname(ds))
      chromo <- basename(ds)
      sub_dt <- merged_dt[merged_dt$hicds == hicds & merged_dt$exprds == exprds & merged_dt$chr == chromo,]
      stopifnot(nrow(sub_dt) > 0)
      stopifnot(is.numeric(c(sub_dt$start, sub_dt$end)))
      sub_dt <- sub_dt[order(sub_dt$start, sub_dt$end),]
      
      
      absDiffValues <- abs(diff(sub_dt[,paste0(diff_var)]))
      scoreValues <- sub_dt[1:(nrow(sub_dt)-1), paste0(score_var)]
      
      stopifnot(length(scoreValues) == length(absDiffValues))
      
      data.frame(
        ds = dirname(ds),
        ds_full = ds,
        absDiffValues = absDiffValues,
        scoreValues = scoreValues,
        stringsAsFactors = FALSE
      )

    }

    plotTit <- paste0("diff. ", diff_var, " vs. ", score_var)
    subTit <- paste0("# ds = ", length(unique(all_ds_dt$ds)))
    
    outFile <- file.path(outFolder, paste0("diff", "_", diff_var, "_vs_", score_var, "_densplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x=all_ds_dt$scoreValues,
      y= all_ds_dt$absDiffValues,
      xlab=paste0(score_var),
      ylab=paste0("abs. diff. ", diff_var),
      main=plotTit,
      cex.main=plotCex,
      cex.axis = plotCex,
      cex.lab = plotCex
    )
    addCorr(x=all_ds_dt$scoreValues, y = all_ds_dt$absDiffValues, bty="n")
    mtext(side=3, text=subTit)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
  

    
}


hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAluad_mutKRAS_mutEGFR"
chromo="chr10"

cum_var = c("adjPvalComb")

all_chromos <- c("chr10", "chr17")

for(chromo in all_chromos){
  
  sub_dt <- merged_dt[merged_dt$hicds == hicds & merged_dt$exprds == exprds & merged_dt$chr == chromo,]
  sub_dt <- sub_dt[order(sub_dt$start, sub_dt$end),]
  
  
  for(cum_var in c("meanLogFC", "meanCorr", "adjPvalComb")) {
    
    
    plotTit <- paste0(hicds, " - ", exprds, " - ", chromo)
    
    subTit <- c(cum_var, " obs. vs. permut. (n=", nPermut, ")")
    
    myylab <- paste0("cumsum diff. ", cum_var)
    myxlab <- paste0("position along ", chromo)
    
    cumsumDiffVar <- cumsum(abs(diff(sub_dt[,cum_var])))
    
    
    cumsumDiffVar_rd <- foreach(i_perm=1:nPermut, .combine='rbind') %dopar% {
      shuff_dt <- sub_dt[sample(1:nrow(sub_dt)),]
      stopifnot(setequal(sub_dt$regionID, shuff_dt$regionID))
      cumsum(abs(diff(shuff_dt[,cum_var])))
    }
    xvar <- 1:length(cumsumDiffVar)
    
    subTit <- paste0(hicds, " - ", exprds, " - ", chromo)
    plotTit <- paste0("cumsum abs. diff. ", cum_var)
    
    outFile <- file.path(outFolder, paste0("cumSumAbsDiff_", cum_var, "_alongGenome_", hicds, "_", exprds, "_", chromo, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
    plot(
      x=xvar,
      y= cumsumDiffVar,
      xlab=myxlab,
      ylab=myylab,
      type="l",
      main=plotTit
    )
    foo <- apply(cumsumDiffVar_rd, 1,lines, x=xvar, col="grey")
    lines(
      x=xvar,
      y= cumsumDiffVar,
      col="red"
    )
    mtext(side=3, text=subTit)
    legend("topleft", lty=1, col=c("red", "grey"), legend=c("obs.", paste0("permut.\nn=", nPermut)), bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
}


