runFolder <- "."
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

hicds_of_interest <- c( "GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "GSE118514_22Rv1_40kb")
exprds_of_interest <- c("TCGAprad_norm_prad")


dt <- resultData[resultData$hicds %in% hicds_of_interest & resultData$exprds %in% exprds_of_interest ,]
dt <- dt[order(dt$adjPvalComb),]

out_dt = do.call(rbind, by(dt, dt$dataset, function(x) {tmp=x[1:3,];tmp$tad_rank <- rank(tmp$adjPvalComb); tmp}))
out_dt$chromo <- gsub("(chr.+)_.+", "\\1", out_dt$region)

write_dt1 <- out_dt[,c("chromo", "start", "end")]

write_dt <- cbind(write_dt1, write_dt1)
write_dt$extra <- "0,0,255"

write.table(write_dt, file="prostate_tads_top3_viz.txt", col.names=F, row.names=F, sep="\t", quote=F)