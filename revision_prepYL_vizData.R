runFolder <- "."
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

hicds_of_interest <- c( "GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "GSE118514_22Rv1_40kb")
# exprds_of_interest <- c("TCGAprad_norm_prad")


# dt <- resultData[resultData$hicds %in% hicds_of_interest & resultData$exprds %in% exprds_of_interest ,]
# dt <- dt[order(dt$adjPvalComb),]
# 
# out_dt = do.call(rbind, by(dt, dt$dataset, function(x) {tmp=x[1:3,];tmp$tad_rank <- rank(tmp$adjPvalComb); tmp}))
# out_dt$chromo <- gsub("(chr.+)_.+", "\\1", out_dt$region)
# 
# write_dt1 <- out_dt[,c("chromo", "start", "end")]
# 
# write_dt <- cbind(write_dt1, write_dt1)
# write_dt$extra <- "0,0,255"
# 
# write.table(write_dt, file="prostate_tads_top3_viz.txt", col.names=F, row.names=F, sep="\t", quote=F)

outFolder <- "TAD_LIST_FOR_HICDS"
dir.create(outFolder, recursive = TRUE)
# Rscript revision_prepYL_vizData.R

for(hicds in hicds_of_interest) {
  
  region_file <- file.path(hicds, "genes2tad", "all_assigned_regions.txt")
  region_dt <- read.delim(region_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo",  "region", "start", "end"))
  # take only TAD regions
  region_dt <- region_dt[grepl("_TAD", region_dt$region),]
  
  out_dt <- region_dt
  
  # out_dt = do.call(rbind, by(dt, dt$dataset, function(x) {tmp=x[1:3,];tmp$tad_rank <- rank(tmp$adjPvalComb); tmp}))
  # out_dt$chromo <- gsub("(chr.+)_.+", "\\1", out_dt$region)
  
  write_dt1 <- out_dt[,c("chromo", "start", "end")]
  
  write_dt <- cbind(write_dt1, write_dt1)
  write_dt$extra <- "0,0,255"
  
  write.table(write_dt, file=file.path(outFolder, paste0(hicds, "_TADs_for_hicdc.txt")), col.names=F, row.names=F, sep="\t", quote=F)
  
}

