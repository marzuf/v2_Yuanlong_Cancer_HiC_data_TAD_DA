# Rscript signif_cmp.R

outFolder <- "SIGNIF_CMP"
dir.create(outFolder, recursive = TRUE)

signifThresh <- 0.01

result_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
result_dt <- result_dt[order(result_dt$adjPvalComb),]
result_dt$dataset <- file.path(result_dt$hicds, result_dt$exprds)
result_dt <- result_dt[result_dt$adjPvalComb <= signifThresh,]
result_dt$id <- file.path(result_dt$hicds, result_dt$exprds, result_dt$region)

partial_result_dt <- get(load("CREATE_FINAL_TABLE_PARTIAL/all_result_dt.Rdata"))
partial_result_dt <- partial_result_dt[order(partial_result_dt$adjPvalComb),]
partial_result_dt$dataset <- file.path(partial_result_dt$hicds, partial_result_dt$exprds)
partial_result_dt <- partial_result_dt[partial_result_dt$adjPvalComb <= signifThresh,]
partial_result_dt$id <- file.path(partial_result_dt$hicds, partial_result_dt$exprds, partial_result_dt$region)

full_only_dt <- result_dt[!result_dt$id %in% partial_result_dt$id,]
full_only_dt <- full_only_dt[order(full_only_dt$hicds, full_only_dt$exprds, full_only_dt$region),]
partial_only_dt <- partial_result_dt[!partial_result_dt$id %in% result_dt$id,]
partial_only_dt <- partial_only_dt[order(partial_only_dt$hicds, partial_only_dt$exprds, partial_only_dt$region),]

outFile <- file.path(outFolder, paste0("signif_byDS_full_only_dt.txt"))
write.table(full_only_dt, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("signif_byDS_partial_only_dt.txt"))
write.table(partial_only_dt, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
cat(paste0("... written: ", outFile, "\n"))
