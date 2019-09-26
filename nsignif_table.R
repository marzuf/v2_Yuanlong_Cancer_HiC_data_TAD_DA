inFile <- "CREATE_FINAL_TABLE/all_result_dt.Rdata"
inDT <- get(load(inFile))

pvalthresh1 <- 0.01
pvalthresh2 <- 0.05

inDT$signifPval001 <- inDT$adjPvalComb <= 0.01
inDT$signifPval005 <- inDT$adjPvalComb <= 0.05

inDT$signifPval001_signifFDR02 <- inDT$signifPval001 & inDT$signifFDR_0.2

dt1 <- aggregate(signifFDR_0.1 ~hicds+exprds, data=inDT, sum)
dt2 <- aggregate(signifFDR_0.2 ~hicds+exprds, data=inDT, sum)
dt3 <- aggregate(signifPval001 ~hicds+exprds, data=inDT, sum)
dt4 <- aggregate(signifPval005 ~hicds+exprds, data=inDT, sum)
dt5 <- aggregate(signifPval001_signifFDR02 ~hicds+exprds, data=inDT, sum)


idcols <- c("hicds", "exprds")

out_dt <- merge(dt5,merge(dt4, merge(dt3, merge(dt1, dt2, by=idcols), by=idcols), by=idcols),by=idcols)

out_dt <- out_dt[order(out_dt$signifPval001, out_dt$signifFDR_0.2, decreasing=TRUE),]

outFolder <- "NSIGNIF_TABLE"
dir.create(outFolder, recursive = TRUE)
outFile <- file.path(outFolder, paste0("nSignifTADs_table.txt"))
write.table(out_dt, file=outFile, col.names = TRUE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... written: ", outFile, "\n"))