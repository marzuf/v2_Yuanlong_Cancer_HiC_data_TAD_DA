outFolder <- file.path("TOPTABLE_PARTIALDIFF")
dir.create(outFolder, recursive = TRUE)

signifTADthresh <- 0.01

full_final_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
full_final_dt <- get(load(full_final_file))

partial_final_file <- file.path("CREATE_FINAL_TABLE_PARTIAL/all_result_dt.Rdata")
partial_final_dt <- get(load(partial_final_file))


full_final_dt_signif <- full_final_dt[full_final_dt$adjPvalComb <= signifTADthresh,]
partial_final_dt_signif <- partial_final_dt[partial_final_dt$adjPvalComb <= signifTADthresh,]

full_final_dt_signif$id <- file.path(full_final_dt_signif$hicds, full_final_dt_signif$exprds, full_final_dt_signif$region)
partial_final_dt_signif$id <- file.path(partial_final_dt_signif$hicds, partial_final_dt_signif$exprds, partial_final_dt_signif$region)

full_final_dt_signif <- full_final_dt_signif[order(full_final_dt_signif$adjPvalComb),]
partial_final_dt_signif <- partial_final_dt_signif[order(partial_final_dt_signif$adjPvalComb),]

top10_full_final_dt_signif <- full_final_dt_signif[1:10,]
top10_partial_final_dt_signif <- partial_final_dt_signif[1:10,]

signifFullOnly_dt <- full_final_dt_signif[!full_final_dt_signif$id %in% partial_final_dt_signif$id,]
signifPartialOnly_dt <- partial_final_dt_signif[! partial_final_dt_signif$id %in% full_final_dt_signif$id,]

top10_signifFullOnly_dt <- top10_full_final_dt_signif[! top10_full_final_dt_signif$id %in% top10_partial_final_dt_signif$id,]
top10_signifPartialOnly_dt <- top10_partial_final_dt_signif[! top10_partial_final_dt_signif$id %in% top10_full_final_dt_signif$id,]

outRegion_signifFullOnly_dt <- aggregate(region~hicds+exprds, data =signifFullOnly_dt, function(x) paste0(x, collapse="/"))
outGenes_signifFullOnly_dt <- aggregate(region_genes~hicds+exprds, data =signifFullOnly_dt, function(x) paste0(x, collapse="/"))
out_signifFullOnly_dt <- merge(outRegion_signifFullOnly_dt, outGenes_signifFullOnly_dt, by=c("hicds", "exprds"))
out_signifFullOnly_dt <- out_signifFullOnly_dt[order(out_signifFullOnly_dt$hicds, out_signifFullOnly_dt$exprds),]

outRegion_signifPartialOnly_dt <- aggregate(region~hicds+exprds, data =signifPartialOnly_dt, function(x) paste0(x, collapse="/"))
outGenes_signifPartialOnly_dt <- aggregate(region_genes~hicds+exprds, data =signifPartialOnly_dt, function(x) paste0(x, collapse="/"))
out_signifPartialOnly_dt <- merge(outRegion_signifPartialOnly_dt, outGenes_signifPartialOnly_dt, by=c("hicds", "exprds"))
out_signifPartialOnly_dt <- out_signifPartialOnly_dt[order(out_signifPartialOnly_dt$hicds, out_signifPartialOnly_dt$exprds),]

outRegion_top10_signifFullOnly_dt <- aggregate(region~hicds+exprds, data =top10_signifFullOnly_dt, function(x) paste0(x, collapse="/"))
outGenes_top10_signifFullOnly_dt <- aggregate(region_genes~hicds+exprds, data =top10_signifFullOnly_dt, function(x) paste0(x, collapse="/"))
out_top10_signifFullOnly_dt <- merge(outRegion_top10_signifFullOnly_dt, outGenes_top10_signifFullOnly_dt, by=c("hicds", "exprds"))
out_top10_signifFullOnly_dt <- out_top10_signifFullOnly_dt[order(out_top10_signifFullOnly_dt$hicds, out_top10_signifFullOnly_dt$exprds),]

outRegion_top10_signifPartialOnly_dt <- aggregate(region~hicds+exprds, data =top10_signifPartialOnly_dt, function(x) paste0(x, collapse="/"))
outGenes_top10_signifPartialOnly_dt <- aggregate(region_genes~hicds+exprds, data =top10_signifPartialOnly_dt, function(x) paste0(x, collapse="/"))
out_top10_signifPartialOnly_dt <- merge(outRegion_top10_signifPartialOnly_dt, outGenes_top10_signifPartialOnly_dt, by=c("hicds", "exprds"))
out_top10_signifPartialOnly_dt <- out_top10_signifPartialOnly_dt[order(out_top10_signifPartialOnly_dt$hicds, out_top10_signifPartialOnly_dt$exprds),]


outFile <- file.path(outFolder, "signifFullOnly_dt.txt")
write.table(out_signifFullOnly_dt, file=outFile, col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE, sep="\t")
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "top10_signifFullOnly_dt.txt")
write.table(out_top10_signifFullOnly_dt, file=outFile, col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE, sep="\t")
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "signifPartialOnly_dt.txt")
write.table(out_signifPartialOnly_dt, file=outFile, col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE, sep="\t")
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "top10_signifPartialOnly_dt.txt")
write.table(out_top10_signifPartialOnly_dt, file=outFile, col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE, sep="\t")
cat(paste0("... written: ", outFile, "\n"))


















