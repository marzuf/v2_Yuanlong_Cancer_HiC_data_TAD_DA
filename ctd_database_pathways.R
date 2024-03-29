plotType <- "png"
myWidth <- 600
myHeight <- 400


# Rscript ctd_database_pathways.R

outFolder <- "CTD_DATABASE_PATHWAYS"
dir.create(outFolder, recursive = TRUE)

# CTD_genes.csv.gz  CTD_genes_diseases.csv.gz  CTD_genes_pathways.csv.gz
setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


gene_tad_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
stopifnot(gene_tad_dt$entrezID %in% names(entrez2symb))
gene_tad_dt$geneSymbol <- entrez2symb[paste0(gene_tad_dt$entrezID)]


result_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))




ctd_genes_pathways_dt <- read.delim("CTD_genes_pathways.csv", header=F, skip=29, sep=",", stringsAsFactors=F,
                                    col.names=c("GeneSymbol","GeneID","PathwayName","PathwayID"))
ctd_genes_pathways_dt$GeneID <- as.character(ctd_genes_pathways_dt$GeneID)


result_dt <- result_dt[order(result_dt$adjPvalComb),]
result_dt$dataset <- file.path(result_dt$hicds, result_dt$exprds)
result_dt_top10 <- do.call(rbind, by(result_dt, result_dt$dataset, function(x) x[1:10,]))
rownames(result_dt_top10) <- NULL
result_dt_top10$id <- file.path(result_dt_top10$hicds, result_dt_top10$exprds, result_dt_top10$region)


gene_tad_dt$id <- file.path(gene_tad_dt$hicds, gene_tad_dt$exprds, gene_tad_dt$region)

gene_tad_dt_top10 <- gene_tad_dt[gene_tad_dt$id %in% result_dt_top10$id,]

#### DO FOR PATHWAYS
require(dplyr)

result_and_pathways_v1 <- inner_join(gene_tad_dt_top10[, c("hicds", "exprds", "region", "id", "entrezID", "geneSymbol", "region", "tad_adjCombPval")],
                                  ctd_genes_pathways_dt[,c("GeneSymbol", "GeneID", "PathwayName")], by=c("entrezID" = "GeneID"))

# 26738
# result_and_pathways_v2 <- inner_join(gene_tad_dt_top10[, c("hicds", "exprds", "region", "id", "entrezID", "geneSymbol", "region", "tad_adjCombPval")],
#                                      ctd_genes_pathways_dt[,c("GeneSymbol", "GeneID", "PathwayName")], by=c("geneSymbol" = "GeneSymbol"))
# nrow(result_and_pathways_v2)
# 24243

nrow(result_and_pathways_v1)
result_and_pathways_dt <- result_and_pathways_v1

countPathways_dt <- data.frame(
  nCount = as.numeric(table(result_and_pathways_dt$PathwayName)),
  pathway = names(table(result_and_pathways_dt$PathwayName)),
  stringsAsFactors = FALSE
)
countPathways_dt <- countPathways_dt[order(countPathways_dt$nCount, decreasing=TRUE),]

outFile <- file.path(outFolder, paste0("barplot_top10pathways_top10tads.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(oma=c(10,1,1,1))
barplot(countPathways_dt$nCount[1:10], names.arg = countPathways_dt$pathway[1:10], las=2)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



