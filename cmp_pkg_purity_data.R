dt1 <- get(load("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER_PKG/aran_pkg/CPE/log10/conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
head(dt1[,c("conserved_region", "intersect_genes_symbol")])
# > head(dt1[,c("conserved_region", "intersect_genes_symbol")])
# conserved_region                                           intersect_genes_symbol
# 35   conserved_region_35                                             AKR1C1,AKR1C2,AKR1C3
# 34   conserved_region_34                                         MT2A,MT1E,MT1F,MT1G,MT1X
# 41   conserved_region_41                   METRN,FAM173A,CCDC78,HAGHL,NARFL,RPUSD1,CHTF18
# 78   conserved_region_78 RHBDL1,STUB1,JMJD8,WDR24,FBXL16,METRN,FAM173A,CCDC78,HAGHL,NARFL
# 102 conserved_region_102                                        FBXL6,SLC52A2,ADCK5,CPSF1
# 70   conserved_region_70                                     MRPS17,GBAS,PSPH,CCT6A,SUMF2


dt2 <- get(load("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER/CPE/log10/conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
head(dt2[,c("conserved_region", "intersect_genes_symbol")])
# > head(dt2[,c("conserved_region", "intersect_genes_symbol")])
# conserved_region                                           intersect_genes_symbol
# 34   conserved_region_34                                             AKR1C1,AKR1C2,AKR1C3
# 33   conserved_region_33                                         MT2A,MT1E,MT1F,MT1G,MT1X
# 42   conserved_region_42                   METRN,FAM173A,CCDC78,HAGHL,NARFL,RPUSD1,CHTF18
# 79   conserved_region_79 RHBDL1,STUB1,JMJD8,WDR24,FBXL16,METRN,FAM173A,CCDC78,HAGHL,NARFL
# 104 conserved_region_104                                        FBXL6,SLC52A2,ADCK5,CPSF1
# 71   conserved_region_71                                     MRPS17,GBAS,PSPH,CCT6A,SUMF2


dt1 <- get(load("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER_PKG/aran_pkg/ESTIMATE/log10/conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
head(dt1[,c("conserved_region", "intersect_genes_symbol")])
# conserved_region                                           intersect_genes_symbol
# 55   conserved_region_55                                             AKR1C1,AKR1C2,AKR1C3
# 36   conserved_region_36                                         MT2A,MT1E,MT1F,MT1G,MT1X
# 42   conserved_region_42                   METRN,FAM173A,CCDC78,HAGHL,NARFL,RPUSD1,CHTF18
# 79   conserved_region_79 RHBDL1,STUB1,JMJD8,WDR24,FBXL16,METRN,FAM173A,CCDC78,HAGHL,NARFL
# 102 conserved_region_102                                        FBXL6,SLC52A2,ADCK5,CPSF1
# 71   conserved_region_71                                     MRPS17,GBAS,PSPH,CCT6A,SUMF2

dt2 <- get(load("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER//log10/conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
head(dt2[,c("conserved_region", "intersect_genes_symbol")])
# conserved_region                                           intersect_genes_symbol
# 55   conserved_region_55                                             AKR1C1,AKR1C2,AKR1C3
# 36   conserved_region_36                                         MT2A,MT1E,MT1F,MT1G,MT1X
# 42   conserved_region_42                   METRN,FAM173A,CCDC78,HAGHL,NARFL,RPUSD1,CHTF18
# 79   conserved_region_79 RHBDL1,STUB1,JMJD8,WDR24,FBXL16,METRN,FAM173A,CCDC78,HAGHL,NARFL
# 102 conserved_region_102                                        FBXL6,SLC52A2,ADCK5,CPSF1
# 71   conserved_region_71                                     MRPS17,GBAS,PSPH,CCT6A,SUMF2


x <- get(load("TAD_MATCHING_PURITYTAGGED_ACROSS_HICDS_ALLMATCH_v2_PKG/aran - pkgESTIMATE/log10/conserved_regions_with_genes_purityTagged_tads_0.05_minBpRatio0.8_minInterGenes3.Rdata"))
x$nDS <- as.numeric(sapply(x$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
x <- x[order(x$nDS, decreasing = T),]
conservedPurityTagged <- x[,c("conserved_region", "nDS", "intersect_genes_symbol")]
rownames(conservedPurityTagged) <- NULL
# head(conservedPurityTagged, 20)
# conserved_region nDS     intersect_genes_symbol
# 1  conserved_region_232  53             C1QA,C1QC,C1QB
# 2   conserved_region_18  52       GIMAP7,GIMAP4,GIMAP6
# 3    conserved_region_4  49       LILRA2,LILRB1,LILRB4
# 4  conserved_region_182  43            SNX20,NOD2,CYLD
# 5  conserved_region_235  41      TAP2,PSMB8,TAP1,PSMB9
# 6  conserved_region_196  39           LIPA,IFIT2,IFIT3
# 7  conserved_region_241  36 HLA-DRB5,HLA-DRB6,HLA-DRB1
# 8  conserved_region_207  35         STAMBPL1,ACTA2,FAS
# 9   conserved_region_73  35  TRIM6,TRIM34,TRIM5,TRIM22
# 10 conserved_region_186  33       GIMAP2,GIMAP1,GIMAP5
# 11 conserved_region_152  30        MS4A6A,MS4A4A,MS4A7
# 12  conserved_region_19  30          AMICA1,MPZL2,CD3E
# 13 conserved_region_240  30        CXCL9,CXCL10,CXCL11
# 14  conserved_region_17  30          ANXA6,CCDC69,GM2A
# 15 conserved_region_215  29     CD33,VSIG10L,ETFB,NKG7
# 16 conserved_region_142  28          APOL4,APOL2,APOL1
# 17 conserved_region_208  27       SLFN11,SLFN12,SLFN13
# 18 conserved_region_195  27         SART3,ISCU,TMEM119
# 19  conserved_region_13  27     CCBL2,RBMXL1,GBP3,GBP1
# 20  conserved_region_61  26 APOBEC3C,APOBEC3F,APOBEC3G


x = get(load("TAD_MATCHING_PURITYTAGGED_ACROSS_HICDS_ALLMATCH_v2/aran/log10/conserved_regions_with_genes_purityTagged_tads_0.05_minBpRatio0.8_minInterGenes3.Rdata"))
x$nDS <- as.numeric(sapply(x$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
x <- x[order(x$nDS, decreasing = T),]
conservedPurityTagged <- x[,c("conserved_region", "nDS", "intersect_genes_symbol")]
rownames(conservedPurityTagged) <- NULL
head(conservedPurityTagged, 20)
# conserved_region nDS     intersect_genes_symbol
# 1  conserved_region_242  53             C1QA,C1QC,C1QB
# 2   conserved_region_22  52       GIMAP7,GIMAP4,GIMAP6
# 3    conserved_region_4  49       LILRA2,LILRB1,LILRB4
# 4  conserved_region_191  43            SNX20,NOD2,CYLD
# 5  conserved_region_244  41      TAP2,PSMB8,TAP1,PSMB9
# 6  conserved_region_249  36 HLA-DRB5,HLA-DRB6,HLA-DRB1
# 7  conserved_region_203  36           LIPA,IFIT2,IFIT3
# 8  conserved_region_214  35         STAMBPL1,ACTA2,FAS
# 9  conserved_region_195  33       GIMAP2,GIMAP1,GIMAP5
# 10  conserved_region_81  32  TRIM6,TRIM34,TRIM5,TRIM22
# 11 conserved_region_164  30        MS4A6A,MS4A4A,MS4A7
# 12  conserved_region_21  30          ANXA6,CCDC69,GM2A
# 13 conserved_region_223  29     CD33,VSIG10L,ETFB,NKG7
# 14  conserved_region_11  29            CD3E,CD3D,UBE4A
# 15 conserved_region_151  28          APOL4,APOL2,APOL1
# 16 conserved_region_202  28         SART3,ISCU,TMEM119
# 17 conserved_region_248  28        CXCL9,CXCL10,CXCL11
# 18  conserved_region_24  28          AMICA1,MPZL2,CD3E
# 19 conserved_region_215  27       SLFN11,SLFN12,SLFN13
# 20  conserved_region_16  25     CCBL2,RBMXL1,GBP3,GBP1


x <- get(load("TAD_MATCHING_PURITYTAGGED_ACROSS_HICDS_ALLMATCH_v2_PKG/aran - pkgCPE/log10/conserved_regions_with_genes_purityTagged_tads_0.05_minBpRatio0.8_minInterGenes3.Rdata"))
x$nDS <- as.numeric(sapply(x$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
x <- x[order(x$nDS, decreasing = T),]
conservedPurityTagged <- x[,c("conserved_region", "nDS", "intersect_genes_symbol")]
rownames(conservedPurityTagged) <- NULL
head(conservedPurityTagged, 20)
# conserved_region nDS     intersect_genes_symbol
# 1  conserved_region_216  53             C1QA,C1QC,C1QB
# 2   conserved_region_18  52       GIMAP7,GIMAP4,GIMAP6
# 3    conserved_region_5  50       LILRA2,LILRB1,LILRB4
# 4  conserved_region_166  40            SNX20,NOD2,CYLD
# 5   conserved_region_76  39             ESM1,GZMK,GZMA
# 6  conserved_region_217  38      TAP2,PSMB8,TAP1,PSMB9
# 7  conserved_region_213  38          VSIG10L,ETFB,NKG7
# 8  conserved_region_179  38           LIPA,IFIT2,IFIT3
# 9  conserved_region_221  34 HLA-DRB5,HLA-DRB6,HLA-DRB1
# 10 conserved_region_175  33       GIMAP2,GIMAP1,GIMAP5
# 11 conserved_region_191  33 ANKRD22,STAMBPL1,ACTA2,FAS
# 12  conserved_region_17  31          ANXA6,CCDC69,GM2A
# 13 conserved_region_154  30        MS4A6A,MS4A4A,MS4A7
# 14 conserved_region_112  30       SLC16A14,SP110,SP140
# 15  conserved_region_28  30         SP110,SP140,SP140L
# 16 conserved_region_192  29       SLFN11,SLFN12,SLFN13
# 17 conserved_region_141  29          APOL4,APOL2,APOL1
# 18 conserved_region_167  29        CXCL9,CXCL10,CXCL11
# 19  conserved_region_13  28     CCBL2,RBMXL1,GBP3,GBP1
# 20  conserved_region_67  28             GFM1,LXN,MFSD1


x = get(load("TAD_MATCHING_PURITYTAGGED_ACROSS_HICDS_ALLMATCH_v2/Aran - CPE//log10/conserved_regions_with_genes_purityTagged_tads_0.05_minBpRatio0.8_minInterGenes3.Rdata"))
x$nDS <- as.numeric(sapply(x$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
x <- x[order(x$nDS, decreasing = T),]
conservedPurityTagged <- x[,c("conserved_region", "nDS", "intersect_genes_symbol")]
rownames(conservedPurityTagged) <- NULL
head(conservedPurityTagged, 20)
# conserved_region nDS     intersect_genes_symbol
# 1  conserved_region_212  53             C1QA,C1QC,C1QB
# 2   conserved_region_18  52       GIMAP7,GIMAP4,GIMAP6
# 3    conserved_region_5  50       LILRA2,LILRB1,LILRB4
# 4  conserved_region_171  41            SNX20,NOD2,CYLD
# 5  conserved_region_213  39      TAP2,PSMB8,TAP1,PSMB9
# 6  conserved_region_209  38          VSIG10L,ETFB,NKG7
# 7   conserved_region_27  37     LIPA,IFIT2,IFIT3,IFIT1
# 8  conserved_region_217  35 HLA-DRB5,HLA-DRB6,HLA-DRB1
# 9   conserved_region_26  34         SP110,SP140,SP140L
# 10 conserved_region_113  34       SLC16A14,SP110,SP140
# 11 conserved_region_181  33       GIMAP2,GIMAP1,GIMAP5
# 12 conserved_region_193  32 ANKRD22,STAMBPL1,ACTA2,FAS
# 13 conserved_region_158  30        MS4A6A,MS4A4A,MS4A7
# 14  conserved_region_65  30             GFM1,LXN,MFSD1
# 15 conserved_region_143  29          APOL4,APOL2,APOL1
# 16 conserved_region_172  29        CXCL9,CXCL10,CXCL11
# 17  conserved_region_17  29          ANXA6,CCDC69,GM2A
# 18  conserved_region_13  28     CCBL2,RBMXL1,GBP3,GBP1
# 19 conserved_region_157  27          CD247,CREG1,RCSD1
# 20 conserved_region_194  27       SLFN11,SLFN12,SLFN13



