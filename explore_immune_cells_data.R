
setDir <- "/media/electron"
setDir <- ""

# Rscript explore_immune_cells_data.R

# IMMUNE_CELLS_DATA/GSE119234_LN_BcellSubset_raw_counts.txt symbols
# IMMUNE_CELLS_DATA/GSE69239_GEOExpression.txt


exprFile <- file.path("IMMUNE_CELLS_DATA", "GSE119234_LN_BcellSubset_raw_counts.txt")
exprFile <- file.path("IMMUNE_CELLS_DATA", "GSE69239_GEOExpression.txt")
exprFile <- file.path("IMMUNE_CELLS_DATA", "GSE115655_BCells.TMM.cpm.txt")

geneExprAggFun <- "median"

source("../MANUSCRIPT_FIGURES/code/TAD_DE_utils.R")
source("../MANUSCRIPT_FIGURES/COCODATA/R/scores_func.R")

outFolder <- file.path("EXPLORE_IMMUNE_CELLS_DATA")
dir.create(outFolder, recursive = TRUE)

tadpos_file <- file.path(setDir, "/mnt/etemp/marie/OLDER_v2_Yuanlong_Cancer_HiC_data_TAD_DA/GM12878_40kb/genes2tad/all_assigned_regions.txt")
tadpos_dt <- read.delim(tadpos_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]

genepos_file <- file.path(setDir, "/mnt/etemp/marie/OLDER_v2_Yuanlong_Cancer_HiC_data_TAD_DA/GM12878_40kb/genes2tad/all_genes_positions.txt")
g2t_dt <- read.delim(genepos_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)

ensembl_file <- file.path(setDir, "/mnt/etemp/marie/MANUSCRIPT_FIGURES/code/data/final_entrez2ensembl.txt")
ensembl_dt <- read.delim(ensembl_file, header=TRUE, stringsAsFactors = FALSE)
ensembl_dt$entrezID <- as.character(ensembl_dt$entrezID)
ensembl_dt$ensemblID <- as.character(ensembl_dt$ensemblID)
entrez2ensembl <- setNames(ensembl_dt$ensemblID, ensembl_dt$entrezID)

dup_ensembl <- ensembl_dt$ensemblID[duplicated(ensembl_dt$ensemblID)]
noEnsemblDup_ensembl_dt <- ensembl_dt[!ensembl_dt$ensemblID %in% dup_ensembl,]
stopifnot(!duplicated(noEnsemblDup_ensembl_dt$ensemblID))
ensembl2entrez <- setNames(noEnsemblDup_ensembl_dt$entrezID, noEnsemblDup_ensembl_dt$ensemblID)

gimap_symbols <- gff_dt$symbol[grepl("^GIMAP", gff_dt$symbol)]
gimap_entrez <- gff_dt$entrezID[grepl("^GIMAP", gff_dt$symbol)]
gimap_tads <- g2t_dt$region[g2t_dt$entrezID %in% gimap_entrez]

if(grepl("GSE119234", exprFile)) {
  outPrefix <- "GSE119234"
  expr_dt <- read.delim(exprFile, header=TRUE, stringsAsFactors = FALSE)  
  colnames(expr_dt)[1] <- "symbol"
  stopifnot(any(expr_dt$symbol %in% gff_dt$symbol))
  cat(paste0("recover entrezID:\t", sum(expr_dt$symbol %in% names(symb2entrez)), "/", nrow(expr_dt), "\n"))
  expr_dt <- expr_dt[expr_dt$symbol %in% names(symb2entrez),]
  expr_dt$entrezID <- symb2entrez[paste0(expr_dt$symbol)]
  stopifnot(expr_dt$entrezID %in% gff_dt$entrezID)
  
  # if duplicated ID, take the mean
  expr_dt$symbol <- NULL
  agg_expr_dt <- aggregate(.~entrezID, data=expr_dt, FUN=geneExprAggFun)
  rownames(agg_expr_dt) <- as.character(agg_expr_dt$entrezID)
  agg_expr_dt$entrezID <- NULL
  qqnorm_expr_dt <- t(apply(agg_expr_dt, 1, quantNorm))
  stopifnot(all(dim(agg_expr_dt) == dim(qqnorm_expr_dt)))
  rownames(qqnorm_expr_dt) <- rownames(agg_expr_dt)
  colnames(qqnorm_expr_dt) <- colnames(agg_expr_dt)
  stopifnot(rownames(qqnorm_expr_dt) %in% gff_dt$entrezID)
  
  g2t_dt <- g2t_dt[g2t_dt$entrezID %in% rownames(qqnorm_expr_dt),]
  stopifnot(g2t_dt$entrezID %in% rownames(qqnorm_expr_dt))
  
  
  
} else if(grepl("GSE69239", exprFile)) {
  outPrefix <- "GSE69239"
  expr_dt <- read.delim(exprFile, header=TRUE, stringsAsFactors = FALSE)  
  expr_dt$ensemblID <- gsub("(.+)\\..+", "\\1", expr_dt$Gene.id)
  stopifnot(any(expr_dt$ensemblID %in% names(ensembl2entrez)))
  expr_dt$entrezID <- ensembl2entrez[paste0(expr_dt$ensemblID)]
  cat(paste0("recover entrezID:\t", sum(is.na(expr_dt$entrezID)) , "/", nrow(expr_dt), "\n"))
  # if duplicated ID, take the mean
  expr_dt$Gene.id <- NULL
  expr_dt$ensemblID <- NULL
  agg_expr_dt <- aggregate(.~entrezID, data=expr_dt, FUN=geneExprAggFun)
  rownames(agg_expr_dt) <- as.character(agg_expr_dt$entrezID)
  agg_expr_dt$entrezID <- NULL
  qqnorm_expr_dt <- t(apply(agg_expr_dt, 1, quantNorm))
  stopifnot(all(dim(agg_expr_dt) == dim(qqnorm_expr_dt)))
  rownames(qqnorm_expr_dt) <- rownames(agg_expr_dt)
  colnames(qqnorm_expr_dt) <- colnames(agg_expr_dt)
  stopifnot(rownames(qqnorm_expr_dt) %in% ensembl_dt$entrezID)
  stopifnot(rownames(qqnorm_expr_dt) %in% gff_dt$entrezID)
  
  g2t_dt <- g2t_dt[g2t_dt$entrezID %in% rownames(qqnorm_expr_dt),]
  stopifnot(g2t_dt$entrezID %in% rownames(qqnorm_expr_dt))
  
  
    
} else if(grepl("GSE115655", exprFile)) {
  
  outPrefix <- "GSE115655"
  
  expr_dt <- read.delim(exprFile, header=TRUE, stringsAsFactors = FALSE)  
  expr_dt$ensemblID <- gsub("(.+)\\..+", "\\1", expr_dt$Ensemble.ID)
  stopifnot(any(expr_dt$ensemblID %in% names(ensembl2entrez)))
  expr_dt$entrezID <- ensembl2entrez[paste0(expr_dt$ensemblID)]
  cat(paste0("recover entrezID:\t", sum(is.na(expr_dt$entrezID)) , "/", nrow(expr_dt), "\n"))
  # if duplicated ID, take the mean
  expr_dt$Ensemble.ID <- NULL
  expr_dt$ensemblID <- NULL
  agg_expr_dt <- aggregate(.~entrezID, data=expr_dt, FUN=geneExprAggFun)
  rownames(agg_expr_dt) <- as.character(agg_expr_dt$entrezID)
  agg_expr_dt$entrezID <- NULL
  qqnorm_expr_dt <- t(apply(agg_expr_dt, 1, quantNorm))
  stopifnot(all(dim(agg_expr_dt) == dim(qqnorm_expr_dt)))
  rownames(qqnorm_expr_dt) <- rownames(agg_expr_dt)
  colnames(qqnorm_expr_dt) <- colnames(agg_expr_dt)
  stopifnot(rownames(qqnorm_expr_dt) %in% ensembl_dt$entrezID)
  stopifnot(rownames(qqnorm_expr_dt) %in% gff_dt$entrezID)
  
  g2t_dt <- g2t_dt[g2t_dt$entrezID %in% rownames(qqnorm_expr_dt),]
  stopifnot(g2t_dt$entrezID %in% rownames(qqnorm_expr_dt))
  

}

qqnorm_expr_dt <- qqnorm_expr_dt[paste0(g2t_dt$entrezID),]
cat(paste0("final number of assigned entrezID:\t", nrow(qqnorm_expr_dt), "\n"))
tad_meanCorr_dt <- get_meanCorr(gene2tad_dt=g2t_dt, exprd_dt=qqnorm_expr_dt)

tad_annot_dt <- tad_meanCorr_dt[tad_meanCorr_dt$region %in% gimap_tads,]
if(nrow(tad_annot_dt) > 0) {
  annot_dt <- merge(tad_annot_dt, g2t_dt[,c("entrezID", "region")], by="region")
  annot_dt$symbol <- entrez2symb[paste0(annot_dt$entrezID)]
  stopifnot(!is.na(annot_dt$symbol))
  agg_annot_dt <- aggregate(symbol~region, FUN=function(x)paste0(x, collapse=","), data=annot_dt)
  legend_dt <- merge(tad_annot_dt, agg_annot_dt, by="region")
  
  ratioGimap <- aggregate(symbol~region, function(x) mean(grepl("GIMAP", x)), data=annot_dt)
  # plot only TADs with >= 1/3 of TADs are GIMAPs
  keep_tads <- ratioGimap$region[ratioGimap$symbol >= 1/3]
  legend_dt <- legend_dt[legend_dt$region %in% keep_tads,]
  
}

outFile <- file.path(outFolder, paste0(outPrefix, "_meanCorr_GIMAP_TAD.png"))
png(outFile, height=400, width=600)
plot(density(tad_meanCorr_dt$meanCorr), main="Intra-TAD corr. and GIMAP TAD(s)", xlab="mean intra-TAD corr.")
if(nrow(tad_annot_dt) > 0) {
  abline(v=legend_dt$meanCorr, col="red")
  text(x=legend_dt$meanCorr, y=max(density(tad_meanCorr_dt$meanCorr)$y), 
       adj=c(1,1),
       col="red", labels=paste0(legend_dt$region, "\n", legend_dt$symbol))
  mtext(side=1, at=legend_dt$meanCorr, text = round(legend_dt$meanCorr, 2), col = "red")

}
mtext(side=3, text=paste0(outPrefix))
foo <- dev.off()


