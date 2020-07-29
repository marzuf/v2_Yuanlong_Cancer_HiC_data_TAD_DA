purity_ds <- "aran_pkg"
pm <- "CPE"

dt0 <- get(load(file.path(runFolder,
                          "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER",
                           pm, "log10", 
                          "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))

dt1 <- get(load(file.path(runFolder,
                                                   "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER_PKG",
                                                   purity_ds, pm, "log10", 
                                                   "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))


head(dt0$intersect_genes_symbol)
head(dt1$intersect_genes_symbol)
# > head(dt0$intersect_genes_symbol)
# [1] "AKR1C1,AKR1C2,AKR1C3"                                             "MT2A,MT1E,MT1F,MT1G,MT1X"                                        
# [3] "METRN,FAM173A,CCDC78,HAGHL,NARFL,RPUSD1,CHTF18"                   "RHBDL1,STUB1,JMJD8,WDR24,FBXL16,METRN,FAM173A,CCDC78,HAGHL,NARFL"
# [5] "FBXL6,SLC52A2,ADCK5,CPSF1"                                        "MRPS17,GBAS,PSPH,CCT6A,SUMF2"                                    
# > head(dt1$intersect_genes_symbol)
# [1] "AKR1C1,AKR1C2,AKR1C3"                                             "MT2A,MT1E,MT1F,MT1G,MT1X"                                        
# [3] "METRN,FAM173A,CCDC78,HAGHL,NARFL,RPUSD1,CHTF18"                   "RHBDL1,STUB1,JMJD8,WDR24,FBXL16,METRN,FAM173A,CCDC78,HAGHL,NARFL"
# [5] "FBXL6,SLC52A2,ADCK5,CPSF1"                                        "MRPS17,GBAS,PSPH,CCT6A,SUMF2"                                    

################################################
purity_ds <- "aran_pkg"
pm <- "ESTIMATE"

dt1 <- get(load(file.path(runFolder,
                          "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER_PKG",
                         purity_ds, pm, "log10", 
                          "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))

dt0 <- get(load(file.path(runFolder,
                          "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER",
                          "log10", 
                          "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))


head(dt0$intersect_genes_symbol)
head(dt1$intersect_genes_symbol)

> head(dt0$intersect_genes_symbol)
[1] "AKR1C1,AKR1C2,AKR1C3"                                             "MT2A,MT1E,MT1F,MT1G,MT1X"                                        
[3] "METRN,FAM173A,CCDC78,HAGHL,NARFL,RPUSD1,CHTF18"                   "RHBDL1,STUB1,JMJD8,WDR24,FBXL16,METRN,FAM173A,CCDC78,HAGHL,NARFL"
[5] "FBXL6,SLC52A2,ADCK5,CPSF1"                                        "MRPS17,GBAS,PSPH,CCT6A,SUMF2"                                    
> head(dt1$intersect_genes_symbol)
[1] "AKR1C1,AKR1C2,AKR1C3"                                             "MT2A,MT1E,MT1F,MT1G,MT1X"                                        
[3] "METRN,FAM173A,CCDC78,HAGHL,NARFL,RPUSD1,CHTF18"                   "RHBDL1,STUB1,JMJD8,WDR24,FBXL16,METRN,FAM173A,CCDC78,HAGHL,NARFL"
[5] "FBXL6,SLC52A2,ADCK5,CPSF1"                                        "MRPS17,GBAS,PSPH,CCT6A,SUMF2"                                    

################################################