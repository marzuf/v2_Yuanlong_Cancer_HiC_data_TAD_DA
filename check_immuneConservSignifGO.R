dt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
colnames(dt1)
# [1] "ID"             "Description"    "GeneRatio"      "BgRatio"        "pvalue"        
# [6] "p.adjust"       "qvalue"         "geneID"         "Count"          "log10_pval"    
# [11] "geneRatio"      "bgRatio"        "foldEnrichment"

dt2 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signiflog10_pval_GO_BP_log10_pval_barplot_dt.Rdata"))
colnames(dt2)
# "log10_pval"  "Description" "labs"        "plot_labs"

dt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/conserved_signif_enrich_resultDT.Rdata"))
colnames(dt1b)
all.equal(dt1, dt1b)

dt1c <- get(load("GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
colnames(dt1c)
all.equal(dt1, dt1c)
all.equal(dt1b, dt1c)


lt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/all_gene_list.Rdata"))
str(lt1)

lt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2/all_gene_list.Rdata"))
str(lt1b)

all.equal(lt1,lt1b)

tad_l1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/go_all_conserved_signif_tads_genes.Rdata"))
str(tad_l1)

tad_l1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/go_all_conserved_signif_tads_genes.Rdata"))
str(tad_l1b)

all.equal(tad_l1,tad_l1b)


enricher_l1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich.Rdata"))
str(enricher_l1)

enricher_l1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/conserved_signif_enrich.Rdata"))
str(enricher_l1b)

all.equal(enricher_l1,enricher_l1b)


dt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT_0.Rdata"))
# colnames(dt1)
dt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/conserved_signif_enrich_resultDT_0.Rdata"))
# colnames(dt1b)
all.equal(dt1, dt1b)

dt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT_1.Rdata"))
# colnames(dt1)
dt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/conserved_signif_enrich_resultDT_1.Rdata"))
# colnames(dt1b)
all.equal(dt1, dt1b)

dt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT_2.Rdata"))
# colnames(dt1)
dt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/conserved_signif_enrich_resultDT_2.Rdata"))
# colnames(dt1b)
all.equal(dt1, dt1b)


dt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT_0.Rdata"))
# colnames(dt1)
dt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/conserved_signif_enrich_resultDT_0.Rdata"))
# colnames(dt1b)
all.equal(dt1, dt1b)
rm(list=ls())

dt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
# colnames(dt1)
dt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/conserved_signif_enrich_resultDT.Rdata"))
# colnames(dt1b)
all.equal(dt1, dt1b) # TRUE
rm(list=ls())

dt1 <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
# colnames(dt1)
dt1c <- get(load("GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
# colnames(dt1b)
all.equal(dt1, dt1c) # FALSE
rm(list=ls())


dt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/conserved_signif_enrich_resultDT.Rdata"))
# colnames(dt1b)
# all.equal(dt1, dt1b)
dt1c <- get(load("GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
# colnames(dt1b)
all.equal(dt1b, dt1c) # TRUE

dt1a <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
dt1a <- dt1a[order(dt1a$ID,dt1a$p.adjust, dt1a$qvalue, dt1a$foldEnrichment),]
dt1b <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_V2/conserved_signif_enrich_resultDT.Rdata"))
dt1b <- dt1b[order(dt1b$ID,dt1b$p.adjust, dt1b$qvalue, dt1b$foldEnrichment),]
dt1c <- get(load("GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
dt1c <- dt1c[order(dt1c$ID,dt1c$p.adjust, dt1c$qvalue, dt1c$foldEnrichment),]

dt1a_id_pvals <- setNames(dt1a$log10_pval, dt1a$ID)
dt1b_id_pvals <- setNames(dt1b$log10_pval, dt1b$ID)
dt1b_id_pvals <- dt1b_id_pvals[paste0(names(dt1a_id_pvals))]
dt1c_id_pvals <- setNames(dt1c$log10_pval, dt1c$ID)
dt1c_id_pvals <- dt1c_id_pvals[paste0(names(dt1a_id_pvals))]

all(round(dt1a_id_pvals,4) == round(dt1b_id_pvals, 4))
all(round(dt1a_id_pvals,4) == round(dt1c_id_pvals, 4))
all(round(dt1b_id_pvals,4) == round(dt1c_id_pvals, 4))

# => OK !

all.equal(dt1a,dt1b)
all.equal(dt1b,dt1c)





