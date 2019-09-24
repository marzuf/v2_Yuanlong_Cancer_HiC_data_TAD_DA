# Rscript prep_homer_file.R

outFolder <- file.path("HOMER")

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$chromo <- gsub("(chr.+)_.+","\\1", final_table_DT$region)


final_table_DT$dataset <- file.path(final_table_DT$hicds, final_table_DT$exprds)
all_ds <- unique(final_table_DT$dataset)

signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
stopifnot(signif_column %in% colnames(final_table_DT))
final_table_DT[,paste0(signifcol)] <- final_table_DT[,paste0(signif_column)] <= signifThresh



signifcols <- c("signifFDR_0.2", signifcol)

#    Column1: Unique Peak ID
#    Column2: chromosome
#    Column3: starting position
#    Column4: ending position
#    Column5: Strand (+/- or 0/1, where 0="+", 1="-")



ds = all_ds[1]
for(ds in all_ds) {
  
  ds_dt <- final_table_DT[final_table_DT$dataset == ds,]
  
  signcol=signifcols[1]
  for(signcol in signifcols) {
    
    
    signif_dt <- ds_dt[ds_dt[,paste0(signcol)],]
    
    if(nrow(signif_dt) == 0) next
    
    out_dt <- signif_dt[,c("region", "chromo", "start", "end")]
    out_dt$strand <- "+"
    
    out_dt <- out_dt[order(out_dt$chromo, out_dt$start, out_dt$end),]
    
    outFile <- file.path(outFolder, dirname(ds), basename(ds), paste0(dirname(ds), "_", basename(ds), "_", signcol, "_plus.txt"))
    dir.create(dirname(outFile), recursive = TRUE)
    write.table(out_dt, file=outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    cat(paste0("... written: ", outFile, "\n"))
    
    out_dt$strand <- "-"
    outFile <- file.path(outFolder, dirname(ds), basename(ds), paste0(dirname(ds), "_", basename(ds), "_", signcol, "_minus.txt"))
    dir.create(dirname(outFile), recursive = TRUE)
    write.table(out_dt, file=outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
  
  
  
}

#warnings()