# Rscript top10_cmp.R

outFolder <- "TOP10_CMP"
dir.create(outFolder, recursive = TRUE)

result_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
result_dt <- result_dt[order(result_dt$adjPvalComb),]
result_dt$dataset <- file.path(result_dt$hicds, result_dt$exprds)
result_dt_top10 <- do.call(rbind, by(result_dt, result_dt$dataset, function(x) x[1:10,]))
rownames(result_dt_top10) <- NULL
result_dt_top10$id <- file.path(result_dt_top10$hicds, result_dt_top10$exprds, result_dt_top10$region)

partial_result_dt <- get(load("CREATE_FINAL_TABLE_PARTIAL/all_result_dt.Rdata"))
partial_result_dt <- partial_result_dt[order(partial_result_dt$adjPvalComb),]
partial_result_dt$dataset <- file.path(partial_result_dt$hicds, partial_result_dt$exprds)
partial_result_dt_top10 <- do.call(rbind, by(partial_result_dt, partial_result_dt$dataset, function(x) x[1:10,]))
rownames(partial_result_dt_top10) <- NULL
partial_result_dt_top10$id <- file.path(partial_result_dt_top10$hicds, partial_result_dt_top10$exprds, partial_result_dt_top10$region)

full_only_dt <- result_dt_top10[!result_dt_top10$id %in% partial_result_dt_top10$id,]
full_only_dt <- full_only_dt[order(full_only_dt$hicds, full_only_dt$exprds, full_only_dt$region),]
partial_only_dt <- partial_result_dt_top10[!partial_result_dt_top10$id %in% result_dt_top10$id,]
partial_only_dt <- partial_only_dt[order(partial_only_dt$hicds, partial_only_dt$exprds, partial_only_dt$region),]

outFile <- file.path(outFolder, paste0("top10_byDS_full_only_dt.txt"))
write.table(full_only_dt, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("top10_byDS_partial_only_dt.txt"))
write.table(partial_only_dt, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
cat(paste0("... written: ", outFile, "\n"))
