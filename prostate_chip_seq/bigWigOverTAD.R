pvalthresh <- 0.05
runFolder <- ".."
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

#   name - name field from bed, which should be unique
#   size - size of bed (sum of exon sizes
#   covered - # bases within exons covered by bigWig
#   sum - sum of values over all bases covered
#   mean0 - average over bases with non-covered bases counting as zeroes
#   mean - average over just covered bases

hicds <- "GSE118514_22Rv1_40kb"
exprds <- "TCGAprad_norm_prad"

sub_result_dt <- resultData[resultData$hicds == hicds & resultData$exprds == exprds,]


bw_over_tad_dt <- read.delim("chip_22Rv1_cover_22Rv1_TADs.bed", 
                             header=F,
                             col.names=c("region", "size", "covered", "sum", "mean0", "mean"),
                             stringsAsFactors = FALSE)
bw_over_tad_dt$region_id <- file.path(hicds, exprds, bw_over_tad_dt$region)

merge_dt <- merge(sub_result_dt, bw_over_tad_dt, by="region_id", all.x=T, all.y=F)
stopifno(!is.na(merge_dt))
merge_dt$signif_lab2 <- gsub("notsignif.+", "notsignif", merge_dt$signif_lab)
ggboxplot(merge_dt,
          y = "mean",
          x="signif_lab2")


