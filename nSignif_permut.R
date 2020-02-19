final_DT <- get(load(file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")))
final_DT_permut <- get(load(file.path("CREATE_FINAL_TABLE_RANDOM", "all_result_dt.Rdata")))
result_dt <-rbind(final_DT, final_DT_permut)

tad_signif_thresh <- 0.01

dt = aggregate(adjPvalComb~hicds+exprds, data=result_dt, FUN=function(x)sum(x<=tad_signif_thresh))

dt[grepl("NCI", dt$hicds) & dt$exprds == "TCGAluad_norm_luad",]

# hicds             exprds adjPvalComb
# 35                ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad           2
# 36      ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_norm_luad           0
# 37   ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_norm_luad           3
# 38 ENCSR489OCU_NCI-H460_RANDOMNBRGENES_40kb TCGAluad_norm_luad          13
# 39    ENCSR489OCU_NCI-H460_RANDOMSHIFT_40kb TCGAluad_norm_luad           4


tad_signif_thresh <- 0.05

dt = aggregate(adjPvalComb~hicds+exprds, data=result_dt, FUN=function(x)sum(x<=tad_signif_thresh))

dt[grepl("NCI", dt$hicds) & dt$exprds == "TCGAluad_norm_luad",]

# hicds             exprds adjPvalComb
# 35                ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad          67
# 36      ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_norm_luad           0
# 37   ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_norm_luad          53
# 38 ENCSR489OCU_NCI-H460_RANDOMNBRGENES_40kb TCGAluad_norm_luad          94
# 39    ENCSR489OCU_NCI-H460_RANDOMSHIFT_40kb TCGAluad_norm_luad          56


x= result_dt[grepl("NCI", result_dt$hicds) & result_dt$exprds == "TCGAluad_norm_luad",]
x <- x[order(x$adjPvalComb),]

outdt=do.call(rbind, 
  by(x, x$hicds, function(y) y[1:5, c("hicds","region", "region_genes", "adjPvalComb")])
)
rownames(outdt)=NULL
outdt$hicds = gsub("ENCSR489OCU_NCI-H460_","", outdt$hicds)
outdt
# hicds       region                                               region_genes  adjPvalComb
# 1                 40kb chr11_TAD390                                           MMP1,MMP12,MMP13 0.0048984992
# 2                 40kb chr12_TAD187                                DNAJC22,SPATS2,TROAP,TUBA1C 0.0048984992
# 3                 40kb chr10_TAD268                                      BEND3P3,SFTPA1,SFTPA2 0.0102473540
# 4                 40kb chr12_TAD199                                     KRT7,KRT80,KRT81,KRT86 0.0102473540
# 5                 40kb  chr15_TAD39                                       ARHGAP11A,GREM1,SCG5 0.0102473540
# 6       PERMUTG2T_40kb chr10_TAD154                                FAM101B,HSD11B1,ITGA1,STON1 0.6094235284
# 7       PERMUTG2T_40kb chr11_TAD364                                          CASKIN2,PALD1,PGC 0.6094235284
# 8       PERMUTG2T_40kb  chr3_TAD161                                   CD300C,LRRC3,MAP6,SEMA3G 0.6094235284
# 9       PERMUTG2T_40kb   chr6_TAD46                                ARHGEF39,CSTF3,FAM65C,MMP12 0.8818093696
# 10      PERMUTG2T_40kb  chr1_TAD100                    AURKB,KIAA1430,LIF,MLXIPL,RARS2,TMEM127 0.9745803646
# 11   RANDOMMIDPOS_40kb  chr7_TAD555 GIMAP1,GIMAP2,GIMAP4,GIMAP5,GIMAP6,GIMAP7,GIMAP8,LOC728743 0.0001829395
# 12   RANDOMMIDPOS_40kb chr17_TAD242                                          ABCA6,ABCA8,ABCA9 0.0003474006
# 13   RANDOMMIDPOS_40kb chr19_TAD201                         LILRA5,LILRA6,LILRB2,LILRB3,LILRB5 0.0023882761
# 14   RANDOMMIDPOS_40kb  chr1_TAD508                                          FCRL2,FCRL3,FCRL5 0.0152479405
# 15   RANDOMMIDPOS_40kb  chr6_TAD127      HLA-DQA1,HLA-DQA2,HLA-DQB1,HLA-DRB1,HLA-DRB5,HLA-DRB6 0.0152479405
# 16 RANDOMNBRGENES_40kb  chr15_TAD96                                  DUOX1,DUOX2,DUOXA1,DUOXA2 0.0002417633
# 17 RANDOMNBRGENES_40kb chr11_TAD212                                        MS4A14,MS4A4A,MS4A7 0.0019474416
# 18 RANDOMNBRGENES_40kb chr19_TAD348                                             FPR1,FPR2,FPR3 0.0020216491
# 19 RANDOMNBRGENES_40kb  chr1_TAD647                                      DTL,INTS7,LPGAT1,NEK2 0.0027402068
# 20 RANDOMNBRGENES_40kb chr16_TAD181                                        MT1E,MT1L,MT1M,MT2A 0.0027402068
# 21    RANDOMSHIFT_40kb  chr7_TAD555                                       GIMAP4,GIMAP7,GIMAP8 0.0001139450
# 22    RANDOMSHIFT_40kb chr17_TAD242                                          ABCA6,ABCA8,ABCA9 0.0002793643
# 23    RANDOMSHIFT_40kb  chr6_TAD129                                  HLA-DOA,HLA-DPA1,HLA-DPB1 0.0005952197
# 24    RANDOMSHIFT_40kb  chr15_TAD38                                       ARHGAP11A,GREM1,SCG5 0.0068198376
# 25    RANDOMSHIFT_40kb  chr1_TAD336                                             GBP1,GBP2,GBP4 0.0103094278
> 