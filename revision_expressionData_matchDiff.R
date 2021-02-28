outFolder <- file.path("REVISION_EXPRESSIONDATA_MATCHDIFF")
dir.create(outFolder, recursive = TRUE)


# Rscript revision_expressionData_matchDiff.R

### CHECK HARD CODED IN REVISION_EXPRESSION DATA!!!
# expr_var <- "FPKM"


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


require(ggpubr)
require(ggsci)
require(doMC)
require(foreach)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myWidth <- 400
myWidthGG <- 7
myHeightGG <- 7
myHeight <- 400

plotCex <- 1.2

tadSignifThresh <- 0.01

# compare the norm tumor FC changes in exprds data and in data retrieved from cell line

ratioAnnotThresh <- 0.75

all_inter_intra1_dt <- get(load("REVISION_INTER_INTRA_PROBA_CORRECTED/all_inter_intra_dt.Rdata"))
all_inter_intra2_dt <- get(load("REVISION_INTER_INTRA_PROBA2_CORRECTED/all_inter_intra_dt.Rdata"))


###################
### PREPARE mrna + SIGNIF DATA
###################

inFolder <- file.path("REVISION_EXPRESSIONDATA")
all_mrna_dt <- get(load(file.path(inFolder, "final_mRNA_dt.Rdata")))

all_mrna_dt$ratioAnnotGenes <- all_mrna_dt$mRNA_nGenes/all_mrna_dt$tad_nGenes

cat(paste0("... discarded: ", sum(all_mrna_dt$ratioAnnotGenes < ratioAnnotThresh), "/", nrow(all_mrna_dt), "\n"))
# ... discarded: 251/15634
filtered_all_mrna_dt <- all_mrna_dt[all_mrna_dt$ratioAnnotGenes >= ratioAnnotThresh,]
nrow(filtered_all_mrna_dt)
# 15383
filtered_all_mrna_dt$regionID2 <- file.path(filtered_all_mrna_dt$hicds, filtered_all_mrna_dt$region)
stopifnot(!duplicated(filtered_all_mrna_dt$regionID2))

# if I have not enough genes -> then set to NA

all_mrna_dt$mRNA_mean[all_mrna_dt$ratioAnnotGenes < ratioAnnotThresh] <- NA
stopifnot(sum(!is.na(all_mrna_dt$mRNA_mean)) == nrow(filtered_all_mrna_dt))
stopifnot(!duplicated(all_mrna_dt$regionID))
mrna_values <- setNames(all_mrna_dt$mRNA_mean, all_mrna_dt$regionID)


###################
### PREPARE SIGNIF DATA
###################

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)
regionID_fcs <- setNames(final_table_DT$meanLogFC, final_table_DT$regionID)


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
stopifnot(matching_withRank_dt$ref_region_ID %in% names(regionID_fcs))
matching_withRank_dt$ref_region_pval <- regionID_pvals[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_pval))
stopifnot(round(matching_withRank_dt$ref_region_pval, 6) == round(matching_withRank_dt$adjPval, 6))

matching_withRank_dt$ref_region_logFC <- regionID_fcs[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_logFC))
stopifnot(! all(matching_withRank_dt$ref_region_logFC == matching_withRank_dt$tumorMeanFC))
stopifnot(matching_withRank_dt$ref_region_logFC == matching_withRank_dt$tumorMeanFC | 
            matching_withRank_dt$ref_region_logFC == matching_withRank_dt$normMeanFC )

stopifnot(matching_withRank_dt$matching_region_ID %in% names(regionID_pvals))
matching_withRank_dt$matching_region_pval <- regionID_pvals[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_pval))

matching_withRank_dt$ref_tadSignif <- ifelse(matching_withRank_dt$adjPval <= tadSignifThresh, "signif.", "not signif.")
my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(matching_withRank_dt$ref_tadSignif))


col_var <- "mRNA_mean"
colVar_values <- mrna_values
newNames <- file.path(dirname(dirname(names(colVar_values))), basename(names(colVar_values)))
stopifnot(!duplicated(newNames))
names(colVar_values) <- newNames

# not true because I have removed those for which not enough gene expression
stopifnot(matching_withRank_dt$ref_region_hicdsID %in% names(colVar_values))
stopifnot(matching_withRank_dt$matching_region_hicdsID %in% names(colVar_values))


matching_withRank_dt[,paste0("ref_", col_var)] <- colVar_values[matching_withRank_dt$ref_region_hicdsID]
matching_withRank_dt[,paste0("matching_", col_var)] <- colVar_values[matching_withRank_dt$matching_region_hicdsID]

# keep only the TADs with enough annotation
# sum(names(colVar_values) %in% filtered_all_mrna_dt$regionID2)
# sum(matching_withRank_dt$ref_region_hicdsID %in% filtered_all_mrna_dt$regionID2)
# this will filter when not enough genes (because set to NA in the introduction)
matching_withRank_dt <- matching_withRank_dt[
  !is.na(  matching_withRank_dt[,paste0("ref_", col_var)] ) & !is.na( matching_withRank_dt[,paste0("matching_", col_var)] ),
]

stopifnot(nrow(matching_withRank_dt) > 0)

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


### NOT NEEDED BECAUSE IN THE EXPRDS THIS ALWAYS TCGA...NORM..TUMOR -> DIRECTION IS ALWAYS THE SAME IRRESPECTIVE OF HICDS
# SO I CAN KEEP matching_withRank_dt$REF_REGION_LOGFC FOR ALL
# matching_withRank_dt[,paste0("norm_vs_tumor_refExprdsFC")] <- ifelse(matching_withRank_dt$ref_hicds %in% all_tumor_ds, 
#                                                            matching_withRank_dt[,paste0("ref_region_logFC")],
#                                                            ifelse(matching_withRank_dt$matching_hicds %in% all_tumor_ds, 
#                                                                   -matching_withRank_dt[,paste0("ref_region_logFC")],NA))
# stopifnot(!is.na(matching_withRank_dt[,paste0("norm_vs_tumor_refExprdsFC")]))

matching_withRank_dt$norm_minus_tumor_cellline_mRNA <- matching_withRank_dt[,paste0("norm_", col_var)] -  
  matching_withRank_dt[,paste0("tumor_", col_var)] 

matching_withRank_dt$norm_over_tumor_cellline_mRNA <- matching_withRank_dt[,paste0("norm_", col_var)] /
  matching_withRank_dt[,paste0("tumor_", col_var)] 


# matching_withRank_dt[matching_withRank_dt$ref_region_ID == "LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12" |
#                        matching_withRank_dt$matching_region_ID == "LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12" ,
#                      c("ref_region_ID", "matching_region_ID","ref_region_logFC", "norm_mRNA_mean", "tumor_mRNA_mean")]
# #                                            ref_region_ID                                 matching_region_ID ref_region_logFC norm_mRNA_mean tumor_mRNA_mean
# 1                  LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12 GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/chr1_TAD10       -0.2967625       12.43167         9.16625
# 10831 GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/chr1_TAD10              LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12       -0.2442317       12.43167         9.16625

all_cmps <- unique(file.path(matching_withRank_dt$matching_hicds, matching_withRank_dt$matching_exprds,
                             matching_withRank_dt$ref_hicds, matching_withRank_dt$ref_exprds))

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_withRank_dt), 
                " (signif.: ", sum(matching_withRank_dt$adjPval <= tadSignifThresh), ")")

plotTit <- "FC from exprds vs. mRNA FPKM diff"

outFile <- file.path(outFolder, "matching_withRank_dt.Rdata")
save(matching_withRank_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

myx <- matching_withRank_dt[,paste0("ref_region_logFC")]  
myy <- log10(matching_withRank_dt$norm_minus_tumor_cellline_mRNA)

toKeep <- !is.na(myx) & !is.na(myy)
myx <- myx[toKeep]
myy <- myy[toKeep]


outFile  <- file.path(outFolder, paste0("exprdsLogFC", "_vs_cellLinesFPKMdiff_densplot.", plotType))
do.call(plotType, list(outFile, height=myWidth, width=myWidth))
densplot(
  x= myx,
  y = myy,
  xlab=paste0("ref_region_logFC"),
  ylab=paste0("norm-tumor mRNA FPKM [log10]"),
  cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
  main=plotTit
)
mtext(side=3, text=mySub)
addCorr(x=myx,y=myy, bty="n", legPos="topright", corMet="spearman")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

myx <- matching_withRank_dt[,paste0("ref_region_logFC")]  
myy <- log2(matching_withRank_dt$norm_over_tumor_cellline_mRNA )
toKeep <- !is.na(myx) & !is.na(myy)
myx <- myx[toKeep]
myy <- myy[toKeep]

plotTit <- "FC from exprds vs. mRNA FPKM ratio"


outFile  <- file.path(outFolder, paste0("exprdsLogFC", "_vs_cellLinesFPKMratio_densplot.", plotType))
do.call(plotType, list(outFile, height=myWidth, width=myWidth))
densplot(
  x= myx,
  y = myy,
  xlab=paste0("ref_region_logFC"),
  ylab=paste0("norm/tumor mRNA FPKM [log2]"),
  cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
  main=plotTit
)
mtext(side=3, text=mySub)
addCorr(x=myx,y=myy, bty="n", legPos="topright", corMet="spearman")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


###################
### PREPARE mrna + proba diff
###################

plot_var <- "mean_intraNorm"

stopifnot(! all_inter_intra1_dt$hicds %in% all_inter_intra2_dt$hicds)
all_inter_intra_dt <- rbind(all_inter_intra1_dt, all_inter_intra2_dt)
stopifnot(final_table_DT$hicds %in% all_inter_intra_dt$hicds)

all_inter_intra_dt$region_hicdsID <- file.path(all_inter_intra_dt$hicds, all_inter_intra_dt$region)
stopifnot(!duplicated(all_inter_intra_dt$region_hicdsID))

plotvar_values <- setNames(all_inter_intra_dt[,paste0(plot_var)], all_inter_intra_dt$region_hicdsID)

stopifnot(matching_withRank_dt$ref_region_hicdsID %in% names(plotvar_values))
stopifnot(matching_withRank_dt$matching_region_hicdsID %in% names(plotvar_values))

matching_withRank_dt[paste0(plot_var, "_ref")] <- plotvar_values[paste0(matching_withRank_dt$ref_region_hicdsID)]
matching_withRank_dt[paste0(plot_var, "_matching")] <- plotvar_values[paste0(matching_withRank_dt$matching_region_hicdsID)]

matching_withRank_dt <- matching_withRank_dt[!is.na(matching_withRank_dt[paste0(plot_var, "_ref")]) &
                                               !is.na(matching_withRank_dt[paste0(plot_var, "_matching")]),]


matching_withRank_dt[,paste0("norm_", plot_var)] <- ifelse(matching_withRank_dt$ref_hicds %in% all_normal_ds, 
                                                          matching_withRank_dt[,paste0(plot_var, "_ref")],
                                                          ifelse(matching_withRank_dt$matching_hicds %in% all_normal_ds, 
                                                                 matching_withRank_dt[,paste0(plot_var,"_matching")],NA))
stopifnot(!is.na(matching_withRank_dt[,paste0("norm_", plot_var)]))

matching_withRank_dt[,paste0("tumor_", plot_var)] <- ifelse(matching_withRank_dt$ref_hicds %in% all_tumor_ds, 
                                                           matching_withRank_dt[,paste0(plot_var,"_ref")],
                                                           ifelse(matching_withRank_dt$matching_hicds %in% all_tumor_ds, 
                                                                  matching_withRank_dt[,paste0(plot_var,"_matching")],NA))
stopifnot(!is.na(matching_withRank_dt[,paste0("tumor_", plot_var)]))

matching_withRank_dt[,paste0("norm_minus_tumor_cellline_", plot_var)] <- matching_withRank_dt[,paste0("norm_", plot_var)] -  
  matching_withRank_dt[,paste0("tumor_", plot_var)] 

matching_withRank_dt[paste0("norm_over_tumor_cellline_", plot_var)] <- matching_withRank_dt[,paste0("norm_", plot_var)] /  
  matching_withRank_dt[,paste0("tumor_", plot_var)] 

myx <- log2(matching_withRank_dt[,paste0("norm_over_tumor_cellline_", plot_var)])
myy <- log2(matching_withRank_dt[,paste0("norm_over_tumor_cellline_mRNA")] )

toKeep <- !is.na(myx) & !is.na(myy)
myx <- myx[toKeep]
myy <- myy[toKeep]

plotTit <- paste0(plot_var, " ratio vs. mRNA FPKM ratio")


outFile  <- file.path(outFolder, paste0("normOverTumor", plot_var,  "_vs_cellLinesFPKMratio_densplot.", plotType))
do.call(plotType, list(outFile, height=myWidth, width=myWidth))
densplot(
  x= myx,
  y = myy,
  xlab=paste0("norm/tumor ", plot_var, " [log2]"),
  ylab=paste0("norm/tumor mRNA FPKM [log2]"),
  cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
  main=plotTit
)
mtext(side=3, text=mySub)
addCorr(x=myx,y=myy, bty="n", legPos="topright", corMet="spearman")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))








  
  
  
  
  
  
  
  
  
  
  
  
  
