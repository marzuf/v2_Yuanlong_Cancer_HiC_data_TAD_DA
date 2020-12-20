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
