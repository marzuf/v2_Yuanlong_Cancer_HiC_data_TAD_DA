plotType <- "png"
myWidth <- 600
myHeight <- 400

startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

# Rscript ctd_database_dataset.R TCGAluad_norm_luad
# Rscript ctd_database_dataset.R TCGAlusc_norm_lusc


args <- commandArgs(trailingOnly = TRUE)
exprds <- args[1]

outFolder <- file.path("CTD_DATABASE", exprds)
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

stopifnot(exprds %in% result_dt$exprds)
result_dt <- result_dt[result_dt$exprds == exprds,]


# ctd_genes_dt <- read.delim("CTD_genes.csv", header=F)

ctd_genes_diseases_dt <- read.delim("CTD_genes_diseases.csv", 
                                    skip=29, header=F, sep=",", stringsAsFactors=F,
                                    col.names=c("GeneSymbol","GeneID","DiseaseName","DiseaseID","DirectEvidence","InferenceChemicalName","InferenceScore","OmimIDs","PubMedIDs"))
ctd_genes_diseases_dt$GeneID <- as.character(ctd_genes_diseases_dt$GeneID)

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

#### DO FOR DISEASES

result_and_diseases_dt <- inner_join(gene_tad_dt_top10[, c("hicds", "exprds", "region", "id", "entrezID", "geneSymbol", "region", "tad_adjCombPval")],
                                     ctd_genes_diseases_dt[,c("GeneSymbol", "GeneID", "DiseaseName")], by=c("entrezID" = "GeneID"))

countDiseases_dt <- data.frame(
  nCount = as.numeric(table(result_and_diseases_dt$DiseaseName)),
  disease = names(table(result_and_diseases_dt$DiseaseName)),
  stringsAsFactors = FALSE
)
countDiseases_dt <- countDiseases_dt[order(countDiseases_dt$nCount, decreasing=TRUE),]

outFile <- file.path(outFolder, paste0("barplot_top10diseases_top10tads.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(oma=c(10,1,1,1))
barplot(countDiseases_dt$nCount[1:10], names.arg = countDiseases_dt$disease[1:10], las=2)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


cat(paste0("... end - ", Sys.time(), "\n"))

