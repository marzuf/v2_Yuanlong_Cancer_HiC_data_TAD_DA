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

out_dt <- merge(dt, cptmt_dt[,c("hicds", "exprds", "region", "tadCptmtLabel", "tadCptmtNormRank","tadCptmtScore")],
                by=c("hicds", "exprds", "region"),all.x=T,all.y=F)
rownames(out_dt) <- NULL
stopifnot(!is.na(out_dt))