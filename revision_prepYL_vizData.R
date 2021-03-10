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

hicds_of_interest <- c( "GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "GSE118514_22Rv1_40kb")
# hicds_of_interest <- c( "GSE118514_RWPE1_40kb")
exprds_of_interest <- c("TCGAprad_norm_prad")


dt <- resultData[resultData$hicds %in% hicds_of_interest & resultData$exprds %in% exprds_of_interest ,]
dt <- dt[order(dt$adjPvalComb),]

out_dt = do.call(rbind, by(dt, dt$dataset, function(x) {tmp=x[1:3,];tmp$tad_rank <- rank(tmp$adjPvalComb); tmp}))
out_dt$chromo <- gsub("(chr.+)_.+", "\\1", out_dt$region)
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

# for(hicds in hicds_of_interest) {
#   
#   region_file <- file.path(hicds, "genes2tad", "all_assigned_regions.txt")
#   region_dt <- read.delim(region_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo",  "region", "start", "end"))
#   # take only TAD regions
#   region_dt <- region_dt[grepl("_TAD", region_dt$region),]
#   
#   out_dt <- region_dt
#   
#   # out_dt = do.call(rbind, by(dt, dt$dataset, function(x) {tmp=x[1:3,];tmp$tad_rank <- rank(tmp$adjPvalComb); tmp}))
#   # out_dt$chromo <- gsub("(chr.+)_.+", "\\1", out_dt$region)
#   
#   write_dt1 <- out_dt[,c("chromo", "start", "end")]
#   
#   write_dt <- cbind(write_dt1, write_dt1)
#   write_dt$extra <- "0,0,255"
#   
#   write.table(write_dt, file=file.path(outFolder, paste0(hicds, "_TADs_for_hicdc.txt")), col.names=F, row.names=F, sep="\t", quote=F)
#   
# }

for(hicds in hicds_of_interest) {
  
  curr_signif_labs <- signif_labs[grepl(exprds_of_interest, names(signif_labs)) & grepl(hicds, names(signif_labs))]
  
  region_file <- file.path(hicds, "genes2tad", "all_assigned_regions.txt")
  region_dt <- read.delim(region_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo",  "region", "start", "end"))
  # take only TAD regions
  region_dt <- region_dt[grepl("_TAD", region_dt$region),]
  region_dt$region_id <- file.path(hicds, exprds_of_interest, region_dt$region)
  stopifnot(names(curr_signif_labs) %in% region_dt$region_id)
  region_dt$signif_lab <- curr_signif_labs[paste0(region_dt$region_id)]
  
  stopifnot(length(na.omit(region_dt$signif_lab)) == 
              sum(resultData$hicds == hicds & resultData$exprds == exprds_of_interest))
  stopifnot(sum(na.omit(region_dt$signif_lab) == "signif.") == 
              sum(resultData$hicds == hicds & resultData$exprds == exprds_of_interest & resultData$signif_lab == "signif."))
  
  stopifnot(setequal(region_dt$region_id[!is.na(region_dt$signif_lab)],
                     resultData$region_id[resultData$hicds == hicds & resultData$exprds == exprds_of_interest]
                     ))
  
  out_dt <- region_dt
  
  # out_dt = do.call(rbind, by(dt, dt$dataset, function(x) {tmp=x[1:3,];tmp$tad_rank <- rank(tmp$adjPvalComb); tmp}))
  # out_dt$chromo <- gsub("(chr.+)_.+", "\\1", out_dt$region)
  
  write_dt1 <- out_dt[,c("chromo", "start", "end")]
  
  write_dt <- cbind(write_dt1, write_dt1)
  write_dt$extra <- "0,0,255"
  stopifnot(file.path(out_dt$chromo,out_dt$start, out_dt$end) ==  file.path(write_dt$chromo,write_dt$start, write_dt$end))
  write_dt$signif_lab <- out_dt$signif_lab
  
  write.table(write_dt, file=file.path(outFolder, paste0(hicds, "_TADs_for_hicdc_withSignif.txt")), col.names=F, row.names=F, sep="\t", quote=F)
  
}
