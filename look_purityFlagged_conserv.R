# Rscript look_purityFlagged_conserv.R

x = get(load("TAD_MATCHING_PURITYTAGGED_ACROSS_HICDS_ALLMATCH_v2/aran/log10/conserved_regions_with_genes_purityTagged_tads_0.05_minBpRatio0.8_minInterGenes3.Rdata"))
x$nDS <- as.numeric(sapply(x$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
x <- x[order(x$nDS, decreasing = T),]
conservedPurityTagged <- x[,c("conserved_region", "nDS", "intersect_genes_symbol")]
rownames(conservedPurityTagged) <- NULL
head(conservedPurityTagged, 20)
nrow(conservedPurityTagged)
conservedPurityTagged_aran <- conservedPurityTagged

x = get(load("TAD_MATCHING_PURITYTAGGED_ACROSS_HICDS_ALLMATCH_v2/EPIC/log10/conserved_regions_with_genes_purityTagged_tads_0.05_minBpRatio0.8_minInterGenes3.Rdata"))
x$nDS <- as.numeric(sapply(x$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
x <- x[order(x$nDS, decreasing = T),]
conservedPurityTagged <- x[,c("conserved_region", "nDS", "intersect_genes_symbol")]
rownames(conservedPurityTagged) <- NULL
head(conservedPurityTagged, 20)
nrow(conservedPurityTagged)
conservedPurityTagged_epic <- conservedPurityTagged

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

png("cmp_nDS_conserved_purityFlagged.png", height=400, width=600)
plot_multiDens(
  list(
       aranData = conservedPurityTagged_aran$nDS,
       epicData = conservedPurityTagged_epic$nDS
  ), my_xlab = "# datasets conserved", plotTit = "Recurrent purityTagged TADs - # datasets conservation"
)
foo <- dev.off()
