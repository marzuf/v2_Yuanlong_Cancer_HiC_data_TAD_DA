require(foreach)
require(doMC)
registerDoMC(40)

# Rscript revision_changes_cptmtLabels.R

setDir <- "/media/electron"
setDir <- ""

outFolder <- "REVISION_CHANGES_CPTMTLABELS"
dir.create(outFolder)

buildTable <- F

hierarchyFolder <- file.path(setDir, 
"/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/2.Results_review/2.Sup_tab_used_for_2nd_submission",
"Supplementary_Data_1_domain_hierarchies_127HiCmaps")
  
matchingFile <- file.path("REVISION_RANKDIFF_ACTIVDIFF", "matching_data.Rdata")

# 
# head /mnt/ndata/Yuanlong/1.Projects/2.PROFILE/2.Results_review/2.Sup_tab_used_for_2nd_submission/Supplementary_Data_1_domain_hierarchies_127HiCmaps/LNCaP_prostate_cancer_binsize=40kb.bed
# #fields: chr, pos_start, pos_end, the full compartment label describing the position in the hierachy, normalized compartment domain rank, ., pos_star
# t, pos_end, color based on eight compartment model (A.1.1 - B.2.2), 1, eight compartment score (A.1.1 to B.2.2: 8 to 1 devided by 8)
# chr1    560001  600000  A.2.1.2.1.2.2.1.1.1.2.2 0.554371002132196       .       560001  600000  #FF9191 1       0.75
# chr1    720001  1080000 A.2.1.2.1.2.2.1.1.1.2.2 0.554371002132196       .       720001  1080000 #FF9191 1       0.75

all_pairs <- c(
  file.path("LI_40kb","GSE105381_HepG2_40kb", "TCGAlihc_norm_lihc"),
  file.path("LG1_40kb" ,"ENCSR444WCZ_A549_40kb", "TCGAluad_norm_luad"),
  file.path("LG2_40kb" ,"ENCSR444WCZ_A549_40kb" ,"TCGAluad_norm_luad"),
  file.path("LG1_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("LG2_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("LG1_40kb", "ENCSR444WCZ_A549_40kb" ,"TCGAlusc_norm_lusc"), 
  file.path("LG2_40kb",  "ENCSR444WCZ_A549_40kb", "TCGAlusc_norm_lusc"), 
  file.path("LG1_40kb","ENCSR489OCU_NCI-H460_40kb", "TCGAlusc_norm_lusc"),
  file.path("LG2_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAlusc_norm_lusc"),
  file.path("GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "TCGAprad_norm_prad"),
  file.path("GSE118514_RWPE1_40kb", "GSE118514_22Rv1_40kb", "TCGAprad_norm_prad")
)
all_normal_ds <- as.character(sapply(all_pairs, function(x) dirname(dirname(x))))
all_tumor_ds <-  as.character(sapply(all_pairs, function(x) basename(dirname(x))))

cptmt_files <-  c(
  "LI_40kb"="LI_liver_binsize=40kb.bed",
    "GSE105381_HepG2_40kb"="HepG2_liver_cancer_binsize=40kb.bed",
"LG1_40kb" ="Lung_tissue_1_binsize=40kb.bed",
  "ENCSR444WCZ_A549_40kb"="A549_alveolar_basal_epithelial_cancer_binsize=40kb.bed",
"LG2_40kb"="Lung_tissue_2_binsize=40kb.bed",
  "ENCSR489OCU_NCI-H460_40kb"="NCI_H460_large_cell_lung_cancer_binsize=40kb.bed",
"GSE118514_RWPE1_40kb"="RWPE1_prostate_epithelium_binsize=40kb.bed",
"ENCSR346DCU_LNCaP_40kb"="LNCaP_prostate_cancer_binsize=40kb.bed",
   "GSE118514_22Rv1_40kb"="22Rv1_prostate_cancer_binsize=40kb.bed"
)


final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_dt <- get(load(final_table_file))
final_table_DT <- final_dt
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)


stopifnot(names(cptmt_files) %in% final_dt$hicds)
sub_final_dt <- final_dt[final_dt$hicds %in% names(cptmt_files),]

hicds = names(cptmt_files)[1]

if(buildTable) {
  tad2cptmt_final_dt <- foreach(hicds = names(cptmt_files), .combine='rbind')%do% {
    
    
    cptmt_dt <- read.delim(file.path(hierarchyFolder, cptmt_files[paste0(hicds)]), header=F, skip=1, stringsAsFactors = FALSE,
                           col.names=c( "chr", "pos_start", "pos_end", "cptmt_label", "normalized_rank",".",
                                        "pos_start", "pos_end", "color", "1", "cptmt_score"))
    stopifnot(cptmt_dt$pos_start.1 == cptmt_dt$pos_start)
    stopifnot(cptmt_dt$pos_end.1 == cptmt_dt$pos_end)
    
    hicds_final_dt <- sub_final_dt[sub_final_dt$hicds == hicds,]  
    stopifnot(nrow(hicds_final_dt) > 0)
    hicds_final_dt$chromo <- gsub("(chr.+)_TAD.+", "\\1", hicds_final_dt$region)
    stopifnot(hicds_final_dt$chromo %in% paste0("chr", 1:22))
    
    hicds_final_dt$midPos <- (hicds_final_dt$start+hicds_final_dt$end)/2
    
    hicds_final_dt$startCptmtLabel <- hicds_final_dt$startCptmtScore <- hicds_final_dt$startCptmtNormRank <- NA
    hicds_final_dt$endCptmtLabel <- hicds_final_dt$endCptmtScore <- hicds_final_dt$endCptmtNormRank <- NA
    hicds_final_dt$midCptmtLabel <- hicds_final_dt$midCptmtScore <- hicds_final_dt$midCptmtNormRank <- NA
    hicds_final_dt$tadCptmtLabel <- hicds_final_dt$tadCptmtScore <- hicds_final_dt$tadCptmtNormRank <- NA
    hicds_final_dt$i_startCptmt <- hicds_final_dt$i_endCptmt <- hicds_final_dt$i_midCptmt <- NA
    i=1
    i=106
    hicds_out_dt <- foreach(i = 1:nrow(hicds_final_dt), .combine='rbind') %dopar% {
      
      chromo <- hicds_final_dt$chromo[i]
      startPos <- hicds_final_dt$start[i]
      endPos <- hicds_final_dt$end[i]
      midPos <- hicds_final_dt$midPos[i]
      
      # !!! there might be gaps !!!
      # 875 chr1 248760001 248920000 B.1.1.2.2.2.2.2.1.1.2.2       0.3274232 .   248760001 248920000
      # 876 chr1 249040001 249240000 B.1.1.2.2.2.2.2.1.1.2.2       0.3274232 .   249040001 249240000
      
      
      # find where start, end, mid fall in
      
      startCptmt <- which(cptmt_dt$chr == chromo & 
                            cptmt_dt$pos_start <= startPos &
                            cptmt_dt$pos_end >= startPos)
      if(length(startCptmt) == 0) startCptmt <- "gap" 
      stopifnot(length(startCptmt) == 1)
      
      midCptmt <- which(cptmt_dt$chr == chromo & 
                          cptmt_dt$pos_start <= midPos &
                          cptmt_dt$pos_end >= midPos)
      if(length(midCptmt) == 0) midCptmt <- "gap" 
      stopifnot(length(midCptmt) == 1)
      
      endCptmt <- which(cptmt_dt$chr == chromo & 
                          cptmt_dt$pos_start <= endPos &
                          cptmt_dt$pos_end >= endPos)
      if(length(endCptmt) == 0) endCptmt <- "gap" 
      stopifnot(length(endCptmt) == 1)
      
      # if start cptmt = mid cptmt -> assign start
      # if end cptmt = mid cptmt -> assign end
      # assign mid otherwise
      tadCptmt <- ifelse(startCptmt == midCptmt, startCptmt,
                         ifelse(endCptmt == midCptmt, endCptmt,midCptmt))
      # for the start
      if(startCptmt == "gap" ) {
        hicds_final_dt$startCptmtLabel[i] <- hicds_final_dt$startCptmtScore[i] <- hicds_final_dt$startCptmtNormRank[i] <- "gap"
        hicds_final_dt$i_startCptmt[i] <- "gap"
      }else {
        hicds_final_dt$startCptmtLabel[i] <- cptmt_dt$cptmt_label[startCptmt] 
        hicds_final_dt$startCptmtScore[i] <- cptmt_dt$cptmt_score[startCptmt] 
        hicds_final_dt$startCptmtNormRank[i] <- cptmt_dt$normalized_rank[startCptmt]
        hicds_final_dt$i_startCptmt[i] <- startCptmt
      }
      # for the end
      if(endCptmt == "gap" ) {
        hicds_final_dt$endCptmtLabel[i] <- hicds_final_dt$endCptmtScore[i] <- hicds_final_dt$endCptmtNormRank[i] <- "gap"
        hicds_final_dt$i_endCptmt[i] <- "gap"
      } else {
        hicds_final_dt$endCptmtLabel[i] <-  cptmt_dt$cptmt_label[endCptmt] 
        hicds_final_dt$endCptmtScore[i] <- cptmt_dt$cptmt_score[endCptmt] 
        hicds_final_dt$endCptmtNormRank[i] <- cptmt_dt$normalized_rank[endCptmt]
        hicds_final_dt$i_endCptmt[i] <- endCptmt
      }
      # for the mid
      if(midCptmt == "gap" ) {
        hicds_final_dt$midCptmtLabel[i] <- hicds_final_dt$midCptmtScore[i] <- hicds_final_dt$midCptmtNormRank[i] <- "gap"
        hicds_final_dt$i_midCptmt[i] <- "gap"
      } else {
        hicds_final_dt$midCptmtLabel[i] <-  cptmt_dt$cptmt_label[midCptmt] 
        hicds_final_dt$midCptmtScore[i] <- cptmt_dt$cptmt_score[midCptmt] 
        hicds_final_dt$midCptmtNormRank[i] <- cptmt_dt$normalized_rank[midCptmt]
        hicds_final_dt$i_midCptmt[i] <- midCptmt
      }
      # for the TAD
      if(tadCptmt == "gap" ) {
        hicds_final_dt$tadCptmtLabel[i] <-  hicds_final_dt$tadCptmtScore[i] <- hicds_final_dt$tadCptmtNormRank[i] <-"gap"
      } else {
        hicds_final_dt$tadCptmtLabel[i] <-  cptmt_dt$cptmt_label[tadCptmt] 
        hicds_final_dt$tadCptmtScore[i] <- cptmt_dt$cptmt_score[tadCptmt] 
        hicds_final_dt$tadCptmtNormRank[i] <- cptmt_dt$normalized_rank[tadCptmt]
      }
      stopifnot(!is.na(hicds_final_dt[i,]))
      hicds_final_dt[i,]
    }
    hicds_out_dt
    
  }
  outFile <- file.path(outFolder, "tad2cptmt_final_dt.Rdata")
  save(tad2cptmt_final_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "tad2cptmt_final_dt.Rdata")
  tad2cptmt_final_dt <- get(load(outFile))
}
# tad2cptmt_final_dt=get(load("REVISION_CHANGES_CPTMTLABELS/tad2cptmt_final_dt.Rdata"))
mean(tad2cptmt_final_dt$tadCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)
# [1] 0.9790391
mean(tad2cptmt_final_dt$tadCptmtLabel == tad2cptmt_final_dt$midCptmtLabel)
# [1] 1
mean(tad2cptmt_final_dt$midCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)
# [1] 0.9790391
mean(tad2cptmt_final_dt$endCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)
# [1] 0.9790391
mean(tad2cptmt_final_dt$endCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)
# [1] 1
mean(tad2cptmt_final_dt$tadCptmtLabel == "gap")
# [1] 0.02096086
mean(tad2cptmt_final_dt$i_endCptmt == tad2cptmt_final_dt$i_startCptmt)
# [1] 0.9390172

tad2cptmt_final_dt$tad_binaryCptmtLab <- substr(tad2cptmt_final_dt$tadCptmtLabel, start=1, stop=1)
tad2cptmt_final_dt$region_ID <- file.path(tad2cptmt_final_dt$hicds, tad2cptmt_final_dt$exprds, tad2cptmt_final_dt$region)
stopifnot(!duplicated(tad2cptmt_final_dt$region_ID))
tad2cptmts <- setNames(tad2cptmt_final_dt$tad_binaryCptmtLab, tad2cptmt_final_dt$region_ID)

# remove the TADs in gaps
tad2cptmt_dt <- tad2cptmt_final_dt[tad2cptmt_final_dt$tadCptmtLabel != "gap",]
stopifnot(tad2cptmt_dt$tad_binaryCptmtLab %in% c("A", "B"))


# matching_data <- get(load("REVISION_RANKDIFF_ACTIVDIFF/matching_data.Rdata"))
matching_data <- get(load(file=matchingFile))

ds1_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["norm_matching_pval_tadRank_dt"]]))
unique(ds1_matching_dt$ref_hicds)
# [1] "LI_40kb"              "LG1_40kb"             "LG2_40kb"             "GSE118514_RWPE1_40kb"
ds2_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["tumor_matching_pval_tadRank_dt"]]))
unique(ds2_matching_dt$ref_hicds)
# [1] "GSE105381_HepG2_40kb"      "ENCSR444WCZ_A549_40kb"     "ENCSR489OCU_NCI-H460_40kb"
# [4] "ENCSR346DCU_LNCaP_40kb"    "GSE118514_22Rv1_40kb"     


matching_withRank_dt <- rbind(ds1_matching_dt, ds2_matching_dt)
rownames(matching_withRank_dt) <- NULL
matching_withRank_dt$ref_chromo <- gsub("(chr.+)_.+", "\\1", matching_withRank_dt$refID )
matching_withRank_dt$matching_chromo <- gsub("(chr.+)_.+", "\\1", matching_withRank_dt$matchingID_maxOverlapBp )
stopifnot( matching_withRank_dt$matching_chromo==matching_withRank_dt$ref_chromo )

stopifnot(matching_withRank_dt$matching_exprds == matching_withRank_dt$ref_exprds )

matching_withRank_dt$ref_region_ID <- file.path(matching_withRank_dt$ref_hicds,
                                                matching_withRank_dt$ref_exprds,
                                                matching_withRank_dt$refID)

matching_withRank_dt$matching_region_ID <- file.path(matching_withRank_dt$matching_hicds,
                                                     matching_withRank_dt$matching_exprds,
                                                     matching_withRank_dt$matchingID_maxOverlapBp)

stopifnot(matching_withRank_dt$ref_region_ID %in% names(regionID_pvals))
matching_withRank_dt$ref_region_pval <- regionID_pvals[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_pval))
stopifnot(round(matching_withRank_dt$ref_region_pval, 6) == round(matching_withRank_dt$adjPval, 6))
stopifnot(matching_withRank_dt$ref_region_ID %in% names(tad2cptmts))
matching_withRank_dt$ref_region_cptmt <- tad2cptmts[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_cptmt))

stopifnot(matching_withRank_dt$matching_region_ID %in% names(regionID_pvals))
matching_withRank_dt$matching_region_pval <- regionID_pvals[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_pval))

stopifnot(matching_withRank_dt$matching_region_ID %in% names(tad2cptmts))
matching_withRank_dt$matching_region_cptmt <- tad2cptmts[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_cptmt))


matching_withRank_dt$ref_tadSignif <- ifelse(matching_withRank_dt$adjPval <= tadSignifThresh, "signif.", "not signif.")

stopifnot(matching_withRank_dt$ref_hicds %in% c(all_tumor_ds, all_normal_ds))
stopifnot(matching_withRank_dt$matching_hicds %in% c(all_tumor_ds, all_normal_ds))

stopifnot(0 == sum(matching_withRank_dt$ref_hicds %in% all_tumor_ds & matching_withRank_dt$matching_hicds %in% all_tumor_ds))
stopifnot(0 == sum(matching_withRank_dt$ref_hicds %in% all_normal_ds & matching_withRank_dt$matching_hicds %in% all_normal_ds))

matching_withRank_dt$normCptmt <- ifelse(matching_withRank_dt$ref_hicds %in% all_normal_ds, matching_withRank_dt$ref_region_cptmt,
                                         ifelse(matching_withRank_dt$matching_hicds %in% all_normal_ds, matching_withRank_dt$matching_region_cptmt,NA))
stopifnot(!is.na(matching_withRank_dt$normCptmt))

matching_withRank_dt$tumorCptmt <- ifelse(matching_withRank_dt$ref_hicds %in% all_tumor_ds, matching_withRank_dt$ref_region_cptmt,
                                         ifelse(matching_withRank_dt$matching_hicds %in% all_tumor_ds, matching_withRank_dt$matching_region_cptmt,NA))
stopifnot(!is.na(matching_withRank_dt$tumorCptmt))

matching_withRank_dt$norm2tumor_cptmtChange <- paste0(matching_withRank_dt$normCptmt, "->", matching_withRank_dt$tumorCptmt)
matching_withRank_dt$tumorMinusNorm_meanLog2FC <- matching_withRank_dt$tumorMeanFC - matching_withRank_dt$normMeanFC

ggboxplot(
  data=matching_withRank_dt, x="norm2tumor_cptmtChange", y="tumorMinusNorm_meanLog2FC"
)

sub_dt <- matching_withRank_dt[matching_withRank_dt$norm2tumor_cptmtChange == "A->B" | 
                                 matching_withRank_dt$norm2tumor_cptmtChange == "B->A",]
ggboxplot(
  data=sub_dt, x="norm2tumor_cptmtChange", y="tumorMinusNorm_meanLog2FC"
)

ggdensity(sub_dt,
          x = paste0("tumorMinusNorm_meanLog2FC"),
          y = "..density..",
          # combine = TRUE,                  # Combine the 3 plots
          xlab = paste0("Rank diff. with matched TAD"),
          # add = "median",                  # Add median line.
          rug = FALSE,                      # Add marginal rug
          color = "norm2tumor_cptmtChange",
          fill = "norm2tumor_cptmtChange",
          palette = "jco"
) 

# same with difference signif. not signif.

# matching_withRank_dt[
#   matching_withRank_dt$matching_region_ID == "LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12" | 
#     matching_withRank_dt$ref_region_ID == "LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12",]
# 
# matchingID_maxOverlapBp      refID            ref_hicds         ref_exprds
# 1                  chr1_TAD10 chr1_TAD12              LI_40kb TCGAlihc_norm_lihc
# 17046              chr1_TAD12 chr1_TAD10 GSE105381_HepG2_40kb TCGAlihc_norm_lihc
# matching_hicds    matching_exprds   adjPval chromo.x start.x   end.x refID_rank
# 1     GSE105381_HepG2_40kb TCGAlihc_norm_lihc 0.4143104     chr1 3400001 3760000  0.6879433
# 17046              LI_40kb TCGAlihc_norm_lihc 0.3482057     chr1 3520001 3800000  0.8377029
# normMeanFC chromo.y start.y   end.y matchingID_rank tumorMeanFC   rankDiff
# 1     -0.2967625     chr1 3520001 3800000       0.8377029  -0.2442317 -0.1497596
# 17046 -0.2967625     chr1 3400001 3760000       0.6879433  -0.2442317  0.1497596
# ref_region_ID
# 1                  LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12
# 17046 GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/chr1_TAD10
# matching_region_ID ref_region_pval matching_region_pval
# 1     GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/chr1_TAD10       0.4143104            0.3482057
# 17046              LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12       0.3482057            0.4143104
# ref_tadSignif
# 1       not signif.
# 17046   not signif.
# # normMeanFC always the same -> correct ! 



my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(matching_withRank_dt$ref_tadSignif))

legTitle <- ""


tumor_ds