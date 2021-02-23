require(foreach)
require(doMC)
registerDoMC(40)

# Rscript revision_changes_cptmtLabels.R

setDir <- "/media/electron"
setDir <- ""

outFolder <- "REVISION_CHANGES_CPTMTLABELS"
dir.create(outFolder)

buildTable <- T

hierarchyFolder <- file.path(setDir, 
"/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/2.Results_review/2.Sup_tab_used_for_2nd_submission",
"Supplementary_Data_1_domain_hierarchies_127HiCmaps")
  
# 
# head /mnt/ndata/Yuanlong/1.Projects/2.PROFILE/2.Results_review/2.Sup_tab_used_for_2nd_submission/Supplementary_Data_1_domain_hierarchies_127HiCmaps/LNCaP_prostate_cancer_binsize=40kb.bed
# #fields: chr, pos_start, pos_end, the full compartment label describing the position in the hierachy, normalized compartment domain rank, ., pos_star
# t, pos_end, color based on eight compartment model (A.1.1 - B.2.2), 1, eight compartment score (A.1.1 to B.2.2: 8 to 1 devided by 8)
# chr1    560001  600000  A.2.1.2.1.2.2.1.1.1.2.2 0.554371002132196       .       560001  600000  #FF9191 1       0.75
# chr1    720001  1080000 A.2.1.2.1.2.2.1.1.1.2.2 0.554371002132196       .       720001  1080000 #FF9191 1       0.75

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
[1] 1
mean(tad2cptmt_final_dt$tadCptmtLabel == "gap")
[1] 0.02096086

# remove the TADs in gaps
tad2cptmt_dt <- tad2cptmt_final_dt[tad2cptmt_final_dt$tadCptmtLabel != "gap",]
tad2cptmt_dt$tad_binaryCptmtLab <- substr(tad2cptmt_dt$tadCptmtLabel, start=1, stop=1)
stopifnot(tad2cptmt_dt$tad_binaryCptmtLab %in% c("A", "B"))



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
