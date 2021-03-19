outFolder <- file.path("PREP_PRAD_DATA")
dir.create(outFolder,recursive = TRUE)

# Rscript prep_prad_data.R

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
dt <- dt[order(dt$hicds, dt$exprds, dt$adjPvalComb),]


cptmt_dt <- get(load("REVISION_EXPRESSIONLEVEL_CPTMTS/tad2cptmt_dt.Rdata"))
stopifnot(cptmt_dt$startCptmtLabel==cptmt_dt$tadCptmtLabel)
head(cptmt_dt)

out_dt <- merge(dt, cptmt_dt[,c("hicds", "exprds", "region", "tad_eightCptmtLab", "tadCptmtNormRank","tadCptmtScore")],
                by=c("hicds", "exprds", "region"),all.x=T,all.y=F)
out_dt <- out_dt[order(out_dt$hicds, out_dt$exprds, out_dt$adjPvalComb),]
rownames(out_dt) <- NULL
stopifnot(!is.na(out_dt))

change_cptmt_dt <- get(load("REVISION_CHANGES_CPTMTLABELS_V2/plot_dt.Rdata"))
change_cptmt_dt <- change_cptmt_dt[change_cptmt_dt$matching_hicds %in% hicds_of_interest & 
                                   change_cptmt_dt$matching_exprds %in% exprds_of_interest &
                                     change_cptmt_dt$ref_hicds %in% hicds_of_interest & 
                                     change_cptmt_dt$ref_exprds %in% exprds_of_interest,]

change_cptmt_dt$hicds <- change_cptmt_dt$ref_hicds
change_cptmt_dt$exprds <- change_cptmt_dt$ref_exprds
change_cptmt_dt$region <- change_cptmt_dt$refID

out_dt2 <- merge(dt, change_cptmt_dt[,c("hicds", "exprds", "region",
                                 "ref_region_cptmt_eight", "ref_region_pval",
                                 "matching_hicds", "matching_exprds",
                                 "matching_region_cptmt_eight", "matching_region_pval", "rankDiff")],
                by=c("hicds", "exprds", "region"),all.x=T,all.y=F)
# stopifnot(!is.na(out_dt2))
which(apply(out_dt2, 1, function(x)any(is.na(x))))
out_dt2 <- out_dt2[order(out_dt2$hicds, out_dt2$exprds, out_dt2$ref_region_pval),]
rownames(out_dt2) <- NULL

out_dt$meanLogFC <- round(out_dt$meanLogFC,4)
out_dt$meanCorr <- round(out_dt$meanCorr,4)
# out_dt$adjPvalComb <- round(out_dt$adjPvalComb,4)
# out_dt$ratioDown <- round(out_dt$ratioDown,4)
out_dt$tadCptmtNormRank <- round(out_dt$tadCptmtNormRank,4)
out_dt$tadCptmtScore <- as.numeric(as.character(out_dt$tadCptmtScore))
# stopifnot(!is.na(out_dt$tadCptmtScore))
out_dt$tadCptmtScore <- round(out_dt$tadCptmtScore,4)




out_dt2$ref_region_pval <- round(out_dt2$ref_region_pval, 4)
out_dt2$matching_region_pval <- round(out_dt2$matching_region_pval, 4)
out_dt2$meanLogFC <- round(out_dt2$meanLogFC, 4)
out_dt2$meanCorr <- round(out_dt2$meanCorr, 4)
# out_dt2$ratioDown <- round(out_dt2$ratioDown,4)
# out_dt2$adjPvalComb <- round(out_dt2$adjPvalComb, 4)
out_dt2$rankDiff <- round(out_dt2$rankDiff, 4)


out_dt$dataset <- file.path(out_dt$hicds, out_dt$exprds)
head_out_dt = do.call(rbind, by(out_dt, out_dt$dataset, function(x) {tmp=x[1:10,];tmp$tad_rank <- rank(tmp$adjPvalComb); tmp}))

out_dt2$dataset <- file.path(out_dt2$hicds, out_dt2$exprds)
# head_out_dt2 = do.call(rbind, by(out_dt2, out_dt2$dataset, function(x) {tmp=x[1:10,];tmp$tad_rank <- rank(tmp$ref_region_pval); tmp}))

head_out_dt2 <- out_dt2[out_dt2$region_id %in% head_out_dt$region_id,]


write.table(out_dt, file=file.path(outFolder, "prad_results.txt"), sep="\t", quote=F, col.names=T, row.names=F)
write.table(out_dt2, file=file.path(outFolder, "prad_results_withMatch.txt"), sep="\t", quote=F, col.names=T, row.names=F)

write.table(head_out_dt, file=file.path(outFolder, "prad_results_top10.txt"), sep="\t", quote=F, col.names=T, row.names=F)
write.table(head_out_dt2, file=file.path(outFolder, "prad_results_withMatch_top10.txt"), sep="\t", quote=F, col.names=T, row.names=F)

