# take the 22rv1 tad cover by 22rv1

bw_over_tad_dt <- read.delim("chip_22Rv1_cover_22Rv1_TADs.bed", 
                             header=F,
                             col.names=c("region", "size", "covered", "sum", "mean0", "mean"),
                             stringsAsFactors = FALSE)
bw_over_tad_dt$region_id <- file.path(hicds, exprds, bw_over_tad_dt$region)


# take the 22rv1 tad cover by rwpe1

merge on region id
take the ratio of the coverage 

# -> compare if the methylation is different between normal and cancer for the signif ones