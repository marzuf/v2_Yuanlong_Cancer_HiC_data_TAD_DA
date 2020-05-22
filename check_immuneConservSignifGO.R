dt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
colnames(dt1)
# [1] "ID"             "Description"    "GeneRatio"      "BgRatio"        "pvalue"        
# [6] "p.adjust"       "qvalue"         "geneID"         "Count"          "log10_pval"    
# [11] "geneRatio"      "bgRatio"        "foldEnrichment"

dt2 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signiflog10_pval_GO_BP_log10_pval_barplot_dt.Rdata"))
colnames(dt2)
# "log10_pval"  "Description" "labs"        "plot_labs"

dt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
colnames(dt1b)
all.equal(dt1, dt1b)