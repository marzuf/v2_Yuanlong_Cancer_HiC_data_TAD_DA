options(scipen=100)

SSHFS=F


# Rscript genes_behind_go.R

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  data_cmpType <- ""
  file_prefix <- ""
} else if(length(args) == 1) {
  data_cmpType <- args[1]  
  stopifnot(data_cmpType %in% c("norm_vs_tumor", "subtypes", "wt_vs_mut"))
  file_prefix <- paste0(data_cmpType, "_")
} else {
  stop("error\n") 
}

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)



outFolder <- file.path("GENES_BEHIND_GO", data_cmpType)
dir.create(outFolder, recursive = TRUE)

inFolder <- file.path("GO_SIGNIF_ACROSS_HICDS_v2", data_cmpType)
stopifnot(dir.exists(inFolder))

inFile <- file.path(inFolder, paste0(file_prefix, "conserved_signif_enrich_resultDT.Rdata"))
stopifnot(file.exists(inFile))
conserved_signif_enrich_resultDT <- get(load(inFile))

inFile <- file.path(inFolder, paste0(file_prefix, "not_conserved_signif_enrich_resultDT.Rdata" ))
stopifnot(file.exists(inFile))
not_conserved_signif_enrich_resultDT <- get(load(inFile))


conservedDT <- conserved_signif_enrich_resultDT[, c("Description", "p.adjust", "geneID")]
conservedDT$geneSymbols <- sapply(1:nrow(conservedDT), function(i) {
  genes_entrez <- unlist(strsplit(x=conservedDT$geneID[i], split="/"))
  stopifnot(genes_entrez %in% names(entrez2symb))
  genes_symbols <- as.character(entrez2symb[paste0(genes_entrez)])
  paste0(sort(genes_symbols), collapse="/")
})
conservedDT <- conservedDT[, c("Description", "p.adjust", "geneSymbols")]

outFile <- file.path(outFolder, "conserved_signif_enrich_resultDT_withSymbols.txt")
write.table(conservedDT, file=outFile, col.names = TRUE, row.names=FALSE, quote=F, append=F, sep="\t")
cat(paste0("... written: ", outFile, "\n"))

not_conservedDT <- not_conserved_signif_enrich_resultDT[, c("Description", "p.adjust", "geneID")]
not_conservedDT$geneSymbols <- sapply(1:nrow(not_conservedDT), function(i) {
  genes_entrez <- unlist(strsplit(x=not_conservedDT$geneID[i], split="/"))
  stopifnot(genes_entrez %in% names(entrez2symb))
  genes_symbols <- as.character(entrez2symb[paste0(genes_entrez)])
  paste0(sort(genes_symbols), collapse="/")
})
not_conservedDT <- not_conservedDT[, c("Description", "p.adjust", "geneSymbols")]
outFile <- file.path(outFolder, "not_conserved_signif_enrich_resultDT_withSymbols.txt")
write.table(not_conservedDT, file=outFile, col.names = TRUE, row.names=FALSE, quote=F, append=F, sep="\t")
cat(paste0("... written: ", outFile, "\n"))




