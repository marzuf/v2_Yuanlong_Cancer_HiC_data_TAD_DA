# Rscript revision_lookGOgenes.R

all_infiles <- c("REVISION_GO_SIGNIF/keepnotPF_entrez_signif_enrich.Rdata",
                 "REVISION_GO_SIGNIF/discardPF_entrez_signif_enrich.Rdata",
                 "REVISION_GO_BYCPTMT/B22_discardPF_entrez_signif_enrich.Rdata",
                 "REVISION_GO_BYCPTMT/B22_keepnotPF_entrez_signif_enrich.Rdata",
                 "REVISION_GO_BYCPTMT_WITHOUTPURITYFILTER/B22_withoutPF_entrez_signif_enrich.Rdata")
infile <- "REVISION_GO_SIGNIF/keepnotPF_entrez_signif_enrich.Rdata"

for(infile in all_infiles) {
  
  outfile <- gsub(".Rdata", "_topTable.txt", infile)
  
  dt=get(load(infile))
  result_dt <- dt@result
  
  nTop <- 10
  sort_col <- "p.adjust"
  
  result_dt <- result_dt[order(result_dt[,paste0(sort_col)]),]
  crop_dt <- result_dt[1:nTop,c("Description", sort_col, "geneID")]
  
  setDir <- "/media/electron"
  setDir <- ""
  entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
  gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
  gff_dt$entrezID <- as.character(gff_dt$entrezID)
  stopifnot(!duplicated(gff_dt$entrezID))
  entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
  
  # x=crop_dt$geneID[1]
  crop_dt$geneSymbol <- sapply(crop_dt$geneID, function(x) {
    all_entrez <- unlist(strsplit(x, split="/"))
    stopifnot(all_entrez %in% names(entrez2symb))
    all_symbols <- entrez2symb[paste0(all_entrez)]
    paste0(all_symbols, collapse="/")
  })
  
  out_dt <- crop_dt
  out_dt$geneID <- NULL
  out_dt[,paste0(sort_col)] <- round(out_dt[,paste0(sort_col)] ,4) 
  write.table(out_dt, file = outfile, col.names=F, row.names=F, sep="\t", quote=F)
  cat(paste0("... written ", outfile, "\n"))
  
}
