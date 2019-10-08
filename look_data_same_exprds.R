load("/media/electron/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/0.Scripts/Generate_TADs_for_Marie/bin_comp_uniq_hq_CELL_LINEs_list_debug.Rdata")

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

source("hicds_cancer_annot.R")


inFile <- "TAD_MATCHING_ACROSS_HICDS/all_matching_dt.Rdata"
matching_dt <- get(load(inFile))

# inFile <- "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2/signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"
# matching_dt <- get(load(inFile))


mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))
ds_dt <- unique(final_dt[,c("hicds", "exprds")])
stopifnot(ds_dt$hicds %in% names(hicds_cancer_annot))
ds_dt$hicds_can <- hicds_cancer_annot[paste0(ds_dt$hicds)]
ds_dt <- ds_dt[grepl("_norm_", ds_dt$exprds),]
ds_dt[order(ds_dt$exprds),]
# hicds             exprds hicds_can
# 15726       ENCSR079VIJ_G401_40kb TCGAkich_norm_kich         1
# 15599      ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich         1
# 32936         GSE99051_786_O_40kb TCGAkich_norm_kich         1
# 75927        GSE105381_HepG2_40kb TCGAlihc_norm_lihc         1
# 24165                     LI_40kb TCGAlihc_norm_lihc         0
# 80121       ENCSR444WCZ_A549_40kb TCGAluad_norm_luad         1
# 36322   ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad         1
# 201210                   LG1_40kb TCGAluad_norm_luad         0
# 1618211                  LG2_40kb TCGAluad_norm_luad         0
# 11716       ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc         1
# 114161  ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc         1
# 108162                   LG1_40kb TCGAlusc_norm_lusc         0
# 116163                   LG2_40kb TCGAlusc_norm_lusc         0
# 3938       ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad         1
# 39430        GSE118514_22Rv1_40kb TCGAprad_norm_prad         1
# 38235        GSE118514_RWPE1_40kb TCGAprad_norm_prad         0
xx <-aggregate(hicds_can~exprds, FUN=function(x) unique(x), data=ds_dt)
toKeepExprds <- xx$exprds[lengths(xx$hicds_can) > 1 ]
# exprds hicds_can
# 1 TCGAkich_norm_kich         1
# 2 TCGAlihc_norm_lihc      1, 0
# 3 TCGAluad_norm_luad      1, 0
# 4 TCGAlusc_norm_lusc      1, 0
# 5 TCGAprad_norm_prad      1, 0



# 75927        GSE105381_HepG2_40kb TCGAlihc_norm_lihc         1
# 24165                     LI_40kb TCGAlihc_norm_lihc         0



# 80121       ENCSR444WCZ_A549_40kb TCGAluad_norm_luad         1
# 201210                   LG1_40kb TCGAluad_norm_luad         0

# 80121       ENCSR444WCZ_A549_40kb TCGAluad_norm_luad         1
# 1618211                  LG2_40kb TCGAluad_norm_luad         0


# 80121       ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad         1
# 201210                   LG1_40kb TCGAluad_norm_luad         0

# 80121       ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad         1
# 1618211                  LG2_40kb TCGAluad_norm_luad         0


# 11716       ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc         1
# 108162                   LG1_40kb TCGAlusc_norm_lusc         0


# 11716       ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc         1
# 116163                   LG2_40kb TCGAlusc_norm_lusc         0


# 114161  ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc         1
# 108162                   LG1_40kb TCGAlusc_norm_lusc         0


# 114161  ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc         1
# 116163                   LG2_40kb TCGAlusc_norm_lusc         0

# 3938       ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad         1
# 38235        GSE118514_RWPE1_40kb TCGAprad_norm_prad         0

# 39430        GSE118514_22Rv1_40kb TCGAprad_norm_prad         1
# 38235        GSE118514_RWPE1_40kb TCGAprad_norm_prad         0
# 

# Rscript rankDiff_activDiff.R <hicds_norm> <hicds_tumor> <exprds>
 Rscript rankDiff_activDiff.R LI_40kb GSE105381_HepG2_40kb TCGAlihc_norm_lihc => ok
 Rscript rankDiff_activDiff.R LG1_40kb ENCSR444WCZ_A549_40kb TCGAluad_norm_luad => ok
 Rscript rankDiff_activDiff.R LG2_40kb ENCSR444WCZ_A549_40kb TCGAluad_norm_luad => ok
 Rscript rankDiff_activDiff.R LG1_40kb ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad => ok
 Rscript rankDiff_activDiff.R LG2_40kb ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad => ok
 Rscript rankDiff_activDiff.R LG1_40kb ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc => ok
 Rscript rankDiff_activDiff.R LG2_40kb ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc => ok
 Rscript rankDiff_activDiff.R LG1_40kb ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc => ok
 Rscript rankDiff_activDiff.R LG2_40kb ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc
 Rscript rankDiff_activDiff.R GSE118514_RWPE1_40kb ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad => ok
 Rscript rankDiff_activDiff.R GSE118514_RWPE1_40kb GSE118514_22Rv1_40kb TCGAprad_norm_prad => ok

