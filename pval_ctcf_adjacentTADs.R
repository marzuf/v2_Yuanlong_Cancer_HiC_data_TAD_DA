# Rscript pval_ctcf_adjacentTADs.R

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



hicds="Barutcu_MCF-10A_40kb"
exprds="TCGAbrca_lum_bas"
chromo="chr21"

sub_dt <- merged_dt[merged_dt$hicds == hicds & merged_dt$exprds == exprds & merged_dt$chr == chromo,]
sub_dt <- sub_dt[order(sub_dt$start, sub_dt$end),]


diff_var = c("meanCorr")

score_var = c("rightHighestChipSeqScore")

for(diff_var in c("meanLogFC", "meanCorr", "adjPvalComb")) {
  
  for(score_var in c("rightHighestMotifScore", "rightHighestChipSeqScore"))
  
  absDiffValues <- abs(diff(sub_dt[,paste0(diff_var)]))
  scoreValues <- sub_dt[1:(nrow(sub_dt)-1), paste0(score_var)]

  stopifnot(length(scoreValues) == length(absDiffValues))
  
  plot(
    x=scoreValues,
    y= absDiffValues
  )
    
}













hicds="Barutcu_MCF-10A_40kb"
exprds="TCGAbrca_lum_bas"
chromo="chr21"

sub_dt <- merged_dt[merged_dt$hicds == hicds & merged_dt$exprds == exprds & merged_dt$chr == chromo,]
sub_dt <- sub_dt[order(sub_dt$start, sub_dt$end),]



cum_var = c("adjPvalComb")

for(cum_var in c("meanLogFC", "meanCorr", "adjPvalComb")) {
  
  
  plotTit <- paste0(hicds, " - ", exprds, " - ", chromo)
  
  subTit <- c(cum_var, " obs. vs. permut. (n=", nPermut, ")")
  
  myylab <- paste0("cumsum diff. ", cum_var)
  myxlab <- paste0("position along ", chromo)
  
  cumsumDiffVar <- cumsum(abs(diff(sub_dt[,cum_var])))
  
  nPermut <- 100
  
  cumsumDiffVar_rd <- foreach(i_perm=1:nPermut, .combine='rbind') %dopar% {
    shuff_dt <- sub_dt[sample(1:nrow(sub_dt)),]
    stopifnot(setequal(sub_dt$regionID, shuff_dt$regionID))
    cumsum(abs(diff(shuff_dt[,cum_var])))
  }
  xvar <- 1:length(cumsumDiffVar)
  
  plot(
    x=xvar,
    y= cumsumDiffVar,
    xlab=myxlab,
    ylab=myylab,
    type="l"
  )
  foo <- apply(cumsumDiffVar_rd, 1,lines, x=xvar, col="grey")
  lines(
    x=xvar,
    y= cumsumDiffVar,
    col="red"
  )
  mtext(side=3, text=subTit)
}


