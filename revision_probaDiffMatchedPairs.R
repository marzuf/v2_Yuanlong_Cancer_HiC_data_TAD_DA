
outFolder <- file.path("REVISION_PROBADIFFMATCHEDPAIRS_CORRECTED")
dir.create(outFolder, recursive = TRUE)
all_inter_intra1_dt <- get(load("REVISION_INTER_INTRA_PROBA_CORRECTED/all_inter_intra_dt.Rdata"))
all_inter_intra2_dt <- get(load("REVISION_INTER_INTRA_PROBA2_CORRECTED/all_inter_intra_dt.Rdata"))

# outFolder <- file.path("REVISION_PROBADIFFMATCHEDPAIRS_V2_CORRECTED")
# dir.create(outFolder, recursive = TRUE)
# all_inter_intra1_dt <- get(load("REVISION_INTER_INTRA_PROBA_V2_CORRECTED/all_inter_intra_dt.Rdata"))
# all_inter_intra2_dt <- get(load("REVISION_INTER_INTRA_PROBA2_V2_CORRECTED/all_inter_intra_dt.Rdata"))


# Rscript revision_probaDiffMatchedPairs.R


require(ggpubr)
require(ggsci)
require(doMC)
require(foreach)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myWidth <- 500
myWidthGG <- 7
myHeightGG <- 5
myHeight <- 400

plotCex <- 1.2

tadSignifThresh <- 0.01


# 2nd part  => for matching pairs
# refID signif. or not signif.
# 
# similar to normFC and tumorFC -> normInterIntraProba tumor
# 
# diff. inter/intra vs. signif/not signif. norm tumor
# diff. inter/intra vs. norm tumor FC

all_pairs <- c(
  file.path("LI_40kb","GSE105381_HepG2_40kb", "TCGAlihc_norm_lihc"),
  file.path("LG1_40kb" ,"ENCSR444WCZ_A549_40kb", "TCGAluad_norm_luad"),
  file.path("LG2_40kb" ,"ENCSR444WCZ_A549_40kb" ,"TCGAluad_norm_luad"),
  file.path("LG1_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("LG2_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "TCGAprad_norm_prad"),
  file.path("GSE118514_RWPE1_40kb", "GSE118514_22Rv1_40kb", "TCGAprad_norm_prad")
)
all_normal_ds <- as.character(sapply(all_pairs, function(x) dirname(dirname(x))))
all_tumor_ds <-  as.character(sapply(all_pairs, function(x) basename(dirname(x))))

all_col_vars <- c("mean_intraNorm")
col_var = "mean_intraNorm"

###################
### PREPARE SIGNIF DATA
###################

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)

###################
### PREPARE PROBA DIFF DATA
###################
stopifnot(! all_inter_intra1_dt$hicds %in% all_inter_intra2_dt$hicds)
all_inter_intra_dt <- rbind(all_inter_intra1_dt, all_inter_intra2_dt)
stopifnot(final_table_DT$hicds %in% all_inter_intra_dt$hicds)
stopifnot(col_var %in% colnames(all_inter_intra_dt))

all_inter_intra_dt$region_hicdsID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)
stopifnot(!duplicated(all_inter_intra_dt$region_hicdsID))

colVar_values <- setNames(all_inter_intra_dt[,paste0(col_var)], all_inter_intra_dt$region_hicdsID)

###################
### PREPARE THE pairs_data
###################

pairFolder <- file.path("REVISION_RANKDIFF_ACTIVDIFF/")
outFile <- file.path(pairFolder, "matching_data.Rdata")
matching_data <- get(load(file=outFile))

ds1_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["norm_matching_pval_tadRank_dt"]]))
unique(ds1_matching_dt$ref_hicds)
# [1] "LI_40kb"              "LG1_40kb"             "LG2_40kb"             "GSE118514_RWPE1_40kb"
ds2_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["tumor_matching_pval_tadRank_dt"]]))
unique(ds2_matching_dt$ref_hicds)
# [1] "GSE105381_HepG2_40kb"      "ENCSR444WCZ_A549_40kb"     "ENCSR489OCU_NCI-H460_40kb"
# [4] "ENCSR346DCU_LNCaP_40kb"    "GSE118514_22Rv1_40kb"     


matching_withRank_dt <- rbind(ds1_matching_dt, ds2_matching_dt)
rownames(matching_withRank_dt) <- NULL

stopifnot(matching_withRank_dt$matching_exprds == matching_withRank_dt$ref_exprds )

matching_withRank_dt$ref_region_ID <- file.path(matching_withRank_dt$ref_hicds,
                                                matching_withRank_dt$ref_exprds,
                                                matching_withRank_dt$refID
)
matching_withRank_dt$ref_region_hicdsID <- file.path(matching_withRank_dt$ref_hicds,
                                                matching_withRank_dt$refID
)

matching_withRank_dt$matching_region_ID <- file.path(matching_withRank_dt$matching_hicds,
                                                     matching_withRank_dt$matching_exprds,
                                                     matching_withRank_dt$matchingID_maxOverlapBp)

matching_withRank_dt$matching_region_hicdsID <- file.path(matching_withRank_dt$matching_hicds,
                                                     matching_withRank_dt$matchingID_maxOverlapBp)

stopifnot(matching_withRank_dt$ref_region_ID %in% names(regionID_pvals))
matching_withRank_dt$ref_region_pval <- regionID_pvals[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_pval))
stopifnot(round(matching_withRank_dt$ref_region_pval, 6) == round(matching_withRank_dt$adjPval, 6))

stopifnot(matching_withRank_dt$matching_region_ID %in% names(regionID_pvals))
matching_withRank_dt$matching_region_pval <- regionID_pvals[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_pval))

matching_withRank_dt$ref_tadSignif <- ifelse(matching_withRank_dt$adjPval <= tadSignifThresh, "signif.", "not signif.")
my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(matching_withRank_dt$ref_tadSignif))


stopifnot(matching_withRank_dt$ref_region_hicdsID %in% names(colVar_values))
stopifnot(matching_withRank_dt$matching_region_hicdsID %in% names(colVar_values))

matching_withRank_dt_s <- matching_withRank_dt

stopifnot(length(all_col_vars) == 1) # will not work with more variables because colVar_values was built to work with 1 variable
foo <- foreach(col_var=all_col_vars) %dopar% {
  matching_withRank_dt <- matching_withRank_dt_s
  matching_withRank_dt[,paste0("ref_", col_var)] <- colVar_values[matching_withRank_dt$ref_region_hicdsID]
  matching_withRank_dt[,paste0("matching_", col_var)] <- colVar_values[matching_withRank_dt$matching_region_hicdsID]

  matching_withRank_dt <- matching_withRank_dt[
    !is.na(  matching_withRank_dt[,paste0("ref_", col_var)] ) & !is.na( matching_withRank_dt[,paste0("matching_", col_var)] ),
  ]
  
  matching_withRank_dt[,paste0("norm_", col_var)] <- ifelse(matching_withRank_dt$ref_hicds %in% all_normal_ds, 
                                                            matching_withRank_dt[,paste0("ref_", col_var)],
                                                            ifelse(matching_withRank_dt$matching_hicds %in% all_normal_ds, 
                                                                   matching_withRank_dt[,paste0("matching_", col_var)],NA))
  stopifnot(!is.na(matching_withRank_dt[,paste0("norm_", col_var)]))
  
  matching_withRank_dt[,paste0("tumor_", col_var)] <- ifelse(matching_withRank_dt$ref_hicds %in% all_tumor_ds, 
                                                             matching_withRank_dt[,paste0("ref_", col_var)],
                                                             ifelse(matching_withRank_dt$matching_hicds %in% all_tumor_ds, 
                                                                    matching_withRank_dt[,paste0("matching_", col_var)],NA))
  stopifnot(!is.na(matching_withRank_dt[,paste0("tumor_", col_var)]))
  
  
  all_cmps <- unique(file.path(matching_withRank_dt$matching_hicds, matching_withRank_dt$matching_exprds,
                               matching_withRank_dt$ref_hicds, matching_withRank_dt$ref_exprds))
  
  mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_withRank_dt), 
                  " (signif.: ", sum(matching_withRank_dt$adjPval <= tadSignifThresh), ")")
  
  plotTit <- ""
  
  
  # plot
  matching_withRank_dt[,paste0("tumorOverNorm_", col_var)] <- matching_withRank_dt[,paste0("tumor_", col_var)] / matching_withRank_dt[,paste0("norm_", col_var)]  
  
  outFile  <- file.path(outFolder, paste0(col_var, "_ratio_vs_refAdjPvalComb_densplot.", plotType))
  do.call(plotType, list(outFile, height=myWidth, width=myWidth))
  
  
  densplot(
    x=matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  ,
    xlab=paste0("tumorOverNorm_", col_var),
    y=matching_withRank_dt[,paste0("ref_region_pval")] ,
    ylab=paste0("ref_region_pval"),
    cex.main=plotCex,
    cex.axis=plotCex,
    cex.lab=plotCex,
    main=plotTit
  )
  mtext(side=3, text=mySub)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  matching_withRank_dt$tumorMinusNormFCdiff <- matching_withRank_dt$tumorMeanFC - matching_withRank_dt$normMeanFC
  
  outFile  <- file.path(outFolder, paste0(col_var, "_ratio_vs_FC_diff_densplot.", plotType))
  do.call(plotType, list(outFile, height=myWidth, width=myWidth))
  
  densplot(
    x=matching_withRank_dt[,paste0("tumorOverNorm_", col_var)]  ,
    xlab=paste0("tumorOverNorm_", col_var),
    y=matching_withRank_dt[,paste0("tumorMinusNormFCdiff")],
    ylab=paste0("tumorMinusNormFCdiff"),
    cex.main=plotCex,
    cex.axis=plotCex,
    cex.lab=plotCex,
    main=plotTit
  )
  mtext(side=3, text=mySub)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0(col_var, "_matching_withRank_dt.Rdata"))
  save(matching_withRank_dt, file=outFile, version=2) 
  cat(paste0("... written: ", outFile, "\n"))
  
}
