dt1 <- read.delim("hgTables_g38_gaps.txt", header=TRUE, stringsAsFactors = FALSE)
dt2 <- read.delim("hgTables_g38_centromeres.txt", header=TRUE, stringsAsFactors = FALSE)

dt2$size <- dt2$chromEnd-dt2$chromStart
dt2$ix <- dt2$n <- dt2$bridge <- NA
dt2$name <- NULL
dt2$type <- "centromere"

stopifnot(setequal(colnames(dt1), colnames(dt2)))

out_dt <- rbind(dt1, dt2[,colnames(dt2)])
write.table(out_dt, file="hgTables_g38_merged.txt", append=F, sep="\t", quote=F, row.names=F, col.names=T)