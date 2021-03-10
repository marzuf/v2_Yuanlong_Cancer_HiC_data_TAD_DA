require(foreach)
require(doMC)
registerDoMC(40)

# Rscript revision_aggChipPeaks_data.R

outFolder <- "REVISION_AGGCHIPPEAKS_DATA"
dir.create(outFolder, recursive = TRUE)
pvalthresh <- 0.01

runFolder <- "."
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
resultData$region_id <- file.path(resultData$dataset, resultData$region)
resultData$signif_lab <- ifelse(resultData$adjPvalComb <= pvalthresh, "signif", "notsignif")
resultData$direction_lab <- ifelse(resultData$meanLogFC <0, "down", 
                                   ifelse(resultData$meanLogFC >0, "up",NA))
stopifnot(!is.na(resultData$direction_lab))
resultData$signif_lab <- paste0(resultData$signif_lab, ".", resultData$direction_lab)
signif_labs <- setNames(resultData$signif_lab, resultData$region_id)
resultData$chromo <- gsub("(chr.+)_TAD.+", "\\1", resultData$region)
stopifnot(resultData$chromo %in% paste0("chr", 1:22))

exprds <- "TCGAprad_norm_prad"

ref_hicds <- "GSE118514_RWPE1_40kb"
match_hicds <- "GSE118514_22Rv1_40kb"

hicds_chips <- c( "GSE118514_RWPE1_40kb" = "RWPE1_ENCFF149SDU",
                "GSE118514_22Rv1_40kb" = "22Rv1_ENCFF655WXZ")

peakFolder <- "prostate_chip_seq"
ref_peaks_file <- file.path(peakFolder, paste0(hicds_chips[paste0(ref_hicds)], "_repPeaks.bed"))
ref_peaks_dt <- read.delim(ref_peaks_file, header=FALSE, col.names = 
                             c("chrom" ,"chromStart", "chromEnd" ,"name","score","strand" ,"signalValue","pValue","qValue","peak"))
stopifnot(is.numeric(ref_peaks_dt$chromStart))
stopifnot(is.numeric(ref_peaks_dt$chromEnd))
stopifnot(resultData$chromo %in% ref_peaks_dt$chrom )
# take only signif. peaks
sum(ref_peaks_dt$qValue < -log10(0.05)) # 18
ref_peaks_dt <- ref_peaks_dt[ref_peaks_dt$qValue >= -log10(0.05),]
# [1] 95495

match_peaks_file <- file.path(peakFolder, paste0(hicds_chips[paste0(match_hicds)], "_repPeaks.bed"))
match_peaks_dt <- read.delim(match_peaks_file, header=FALSE, col.names = 
                             c("chrom" ,"chromStart", "chromEnd" ,"name","score","strand" ,"signalValue","pValue","qValue","peak"))
stopifnot(is.numeric(match_peaks_dt$chromStart))
stopifnot(is.numeric(match_peaks_dt$chromEnd))
stopifnot(resultData$chromo %in% match_peaks_dt$chrom )
# take only signif. peaks
sum(match_peaks_dt$qValue < -log10(0.05)) # 21
match_peaks_dt <- match_peaks_dt[match_peaks_dt$qValue >= -log10(0.05),]
# [1] 45955
# 

ref_final_dt <- resultData[resultData$hicds == ref_hicds & 
                             resultData$exprds == exprds, ]
stopifnot(nrow(ref_final_dt)> 0)

ref_final_dt$match_meanQvalue <-ref_final_dt$match_nbrPeaks <- ref_final_dt$ref_meanQvalue <- ref_final_dt$ref_nbrPeaks <- NA



i_tad=1
tad_withPeaks_dt <- foreach(i_tad = 1:nrow(ref_final_dt), .combine='rbind') %dopar% {
  
  if(i_tad==4) save(ref_final_dt, file="ref_final_dt.Rdata", version=2)
  if(i_tad==4) save(match_peaks_dt, file="match_peaks_dt.Rdata", version=2)
  if(i_tad==4) save(ref_peaks_dt, file="ref_peaks_dt.Rdata", version=2)
  
  
  curr_chromo <- ref_final_dt$chromo[i_tad]
  curr_start <- ref_final_dt$start[i_tad]
  curr_end <- ref_final_dt$end[i_tad]

  sub_match_dt <- match_peaks_dt[match_peaks_dt$chrom == curr_chromo &
                               match_peaks_dt$start >= curr_start &
                               match_peaks_dt$end <= curr_end ,
                               ]
  # if(nrow(sub_match_dt) > 0 ) {
    ref_final_dt$match_meanQvalue[i_tad] <- mean(sub_match_dt$qValue)
    ref_final_dt$match_nbrPeaks[i_tad] <- nrow(sub_match_dt)
    
  # } else {
  #   ref_final_dt$match_meanQvalue[i_tad] <- NA
  #   ref_final_dt$match_nbrPeaks[i_tad] <- NA
  #   
  # }
  
  
  sub_ref_dt <- ref_peaks_dt[ref_peaks_dt$chrom == curr_chromo &
                               ref_peaks_dt$start >= curr_start &
                               ref_peaks_dt$end <= curr_end ,
  ]
  
  # if(nrow(sub_match_dt) > 0 ) {
    ref_final_dt$ref_meanQvalue[i_tad] <- mean(sub_ref_dt$qValue)
    ref_final_dt$ref_nbrPeaks[i_tad] <- nrow(sub_ref_dt)
  #   
  #   
  # } else {
  #   ref_final_dt$ref_meanQvalue[i_tad] <- NA
  #   ref_final_dt$ref_nbrPeaks[i_tad] <- NA
  #   
  # }
  # 
  
  
  ref_final_dt[i_tad,]
}

# tad_withPeaks_dt <- ref_final_dt
outFile <- file.path(outFolder, paste0(ref_hicds, "_TADs_match_", match_hicds, "_", "tad_withPeaks_dt.Rdata"))
save(tad_withPeaks_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))