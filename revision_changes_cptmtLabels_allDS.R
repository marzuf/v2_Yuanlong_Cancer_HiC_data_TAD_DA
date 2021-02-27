require(foreach)
require(doMC)
registerDoMC(40)
require(ggsci)
require(ggpubr)
require(ggplot2)

# Rscript revision_changes_cptmtLabels_allDS.R

setDir <- "/media/electron"
setDir <- ""

plotType <- "svg"
myHeightGG <- 5
myWidthGG <- 7


mytheme <-     theme(
  # text = element_text(family=fontFamily),
  panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
  panel.background = element_rect(fill = "transparent"),
  panel.grid.major.x =  element_blank(),
  panel.grid.minor.x =  element_blank(),
  axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
  axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
  axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
  plot.title = element_text(hjust=0.5, size = 16, face="bold"),
  plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
  legend.title = element_text(face="bold")
) 

outFolder <- "REVISION_CHANGES_CPTMTLABELS_ALLDS"
dir.create(outFolder)

buildTable <- TRUE

hierarchyFolder <- file.path(setDir, 
"/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/2.Results_review/2.Sup_tab_used_for_2nd_submission",
"Supplementary_Data_1_domain_hierarchies_127HiCmaps")
  
matchingFile <- file.path("REVISION_RANKDIFF_ACTIVDIFF", "matching_data.Rdata")

tadSignifThresh <- 0.01

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
   "GSE118514_22Rv1_40kb"="22Rv1_prostate_cancer_binsize=40kb.bed",
"Barutcu_MCF-10A_40kb"="MCF_10A_breast_epithelium_binsize=40kb.bed",
"Barutcu_MCF-7_40kb"="MCF_7_breast_cancer_binsize=40kb.bed" ,
"ENCSR079VIJ_G401_40kb"="G401_rhabdoid_tumor_kidney_binsize=40kb.bed",
"ENCSR312KHQ_SK-MEL-5_40kb"="SK_MEL_5_melanoma_binsize=40kb.bed",      
"ENCSR401TBQ_Caki2_40kb"="Caki2_renal_cell_cancer_binsize=40kb.bed",
"ENCSR504OTV_transverse_colon_40kb" ="Transverse_colon_binsize=40kb.bed",
"ENCSR549MGQ_T47D_40kb"="T47D_breast_cancer_binsize=40kb.bed"    ,
"ENCSR862OGI_RPMI-7951_40kb"="RPMI_7951_melanoma_binsize=40kb.bed",
"GSE105194_cerebellum_40kb"="Cerebellum_binsize=40kb.bed",
 "GSE105194_spinal_cord_40kb"="Spinal_cord_binsize=40kb.bed",
"GSE105318_DLD1_40kb"="DLD1_colorectal_cancer_binsize=40kb.bed",
"GSE109229_BT474_40kb"="BT474_breast_cancer_binsize=40kb.bed",
"GSE109229_SKBR3_40kb"="SKBR3_breast_cancer_binsize=40kb.bed",
"GSE118588_Panc_beta_40kb"= "Pancreas_beta_binsize=40kb.bed",
 "GSE99051_786_O_40kb"="786_O_renal_cancer_binsize=40kb.bed",
"HMEC_40kb"="HMEC_mammary_epithelium_binsize=40kb.bed",
 "K562_40kb"="K562_ML_binsize=40kb.bed",
 "PA2_40kb" ="Pancreas_1_binsize=40kb.bed",
"PA3_40kb"  ="Pancreas_2_binsize=40kb.bed",    
 "Panc1_rep12_40kb" ="Pancreatic_cancer_binsize=40kb.bed",
"Rao_HCT-116_2017_40kb"="HCT_116_colon_cancer_binsize=40kb.bed"
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
sum(tad2cptmt_final_dt$startCptmtLabel == "gap")
# [1] 0
sum(tad2cptmt_final_dt$endCptmtLabel == "gap")
# [1] 0

stopifnot(sum(tad2cptmt_final_dt$startCptmtLabel == "gap") == 0)
stopifnot(sum(tad2cptmt_final_dt$endCptmtLabel == "gap") == 0)
stopifnot(tad2cptmt_final_dt$endCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)

### UPDATE HERE THE TADS ARE THE COMPARTMENTS... ONLY MIDCPTMT MIGHT BE IN GAP
# retaint the start (= end ) compartment
tad2cptmt_final_dt$tadCptmtLabel <- tad2cptmt_final_dt$startCptmtLabel
tad2cptmt_final_dt$tadCptmtNormRank <- tad2cptmt_final_dt$startCptmtNormRank

# from A.1.1 to B.2.2
tad2cptmt_final_dt$tad_eightCptmtLab <- substr(tad2cptmt_final_dt$tadCptmtLabel, start=1, stop=5)

tad2cptmt_final_dt$tad_binaryCptmtLab <- substr(tad2cptmt_final_dt$tadCptmtLabel, start=1, stop=1)
tad2cptmt_final_dt$region_ID <- file.path(tad2cptmt_final_dt$hicds, tad2cptmt_final_dt$exprds, tad2cptmt_final_dt$region)
stopifnot(!duplicated(tad2cptmt_final_dt$region_ID))
tad2cptmts <- setNames(tad2cptmt_final_dt$tad_binaryCptmtLab, tad2cptmt_final_dt$region_ID)
tad2cptmts_full <- setNames(tad2cptmt_final_dt$tadCptmtLabel, tad2cptmt_final_dt$region_ID)
tad2cptmts_eight <- setNames(tad2cptmt_final_dt$tad_eightCptmtLab, tad2cptmt_final_dt$region_ID)

# remove the TADs in gaps
tad2cptmt_dt <- tad2cptmt_final_dt[tad2cptmt_final_dt$tadCptmtLabel != "gap",]
stopifnot(tad2cptmt_dt$tad_binaryCptmtLab %in% c("A", "B"))
### UPDATE HERE THE TADS ARE THE COMPARTMENTS... ONLY MIDCPTMT MIGHT BE IN GAP
# START CPTMT = ENDCPTMT, AND NEITHER START NOR END ARE GAP

stopifnot(nrow(tad2cptmt_dt) == nrow(tad2cptmt_final_dt))
outFile <- file.path(outFolder, "tad2cptmt_dt.Rdata")
save(tad2cptmt_dt, file=outFile, version=2)
stop("-ok \n")


######################## PLOT ALL DATASETS
# this can be done for all datasets -> cptmt assignment

all_ds_dt <- tad2cptmt_final_dt
all_ds_dt$tadSignif <- ifelse(all_ds_dt$adjPvalComb <= tadSignifThresh, "signif.", "not signif.")

my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(all_ds_dt$tadSignif))
legTitle <- ""

stopifnot(!duplicated(all_ds_dt$region_ID))

all_cmps <- unique(file.path(all_ds_dt$hicds, all_ds_dt$exprds))

mySub <- paste0("# DS = ", length(all_cmps), "; # TADs = ", nrow(all_ds_dt), 
                " (signif.: ", sum(all_ds_dt$adjPvalComb <= tadSignifThresh), ")")

all_ds_dt$region_cptmt <- all_ds_dt$tad_binaryCptmtLab
all_ds_dt$region_cptmt_eight <- all_ds_dt$tad_eightCptmtLab

  
for(suffix in c("", "_eight")) {
  ################# BARPLOT DISTRIBUTION CPTMT CATEGORY SIGNIF/NOT SIGNIF - BINARY
  
  # barplot: ratio of b->a ou a->b changes in signif. and not signif
  agg_dt <- aggregate(as.formula(paste0("region_ID ~ region_cptmt", suffix, "+tadSignif")),
                      data=all_ds_dt, FUN=length)
  colnames(agg_dt)[colnames(agg_dt)=="region_ID"] <- "nTADs"
  tmp_dt <- aggregate(region_ID~tadSignif, data=all_ds_dt, FUN=length)
  nTot <- setNames(tmp_dt$region_ID, tmp_dt$tadSignif)
  agg_dt$ratioTADs <- agg_dt$nTADs/nTot[paste0(agg_dt$tadSignif)]
  plotTit <- paste0("ratio TADs by cptmt", suffix, " (all DS)")
  ggbar_p <-  ggbarplot(agg_dt, 
                        y="ratioTADs",
                        x="tadSignif", 
                        fill=paste0("region_cptmt", suffix)) +
    ggtitle(plotTit, subtitle=mySub)+
    mytheme +
    labs(x="" , y ="ratio of TADs", color=paste0(legTitle),fill=paste0(legTitle)) + 
    theme(
      axis.text.x = element_text(hjust=1, vjust=0.5,size=10,angle=90)
    )
  outFile <- file.path(outFolder,paste0("cptmt", suffix, "_signif_notSignif_allDS_barplot.", plotType))
  ggsave(ggbar_p, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ################# BARPLOT DISTRIBUTION CPTMT CHANGE CATEGORY SIGNIF/NOT SIGNIF - BINARY
  
  # barplot: ratio of b->a ou a->b changes in signif. and not signif
  agg_dt <- aggregate(as.formula(paste0("region_ID ~ region_cptmt", suffix, "+tadSignif")),
                      data=all_ds_dt, FUN=length)
  colnames(agg_dt)[colnames(agg_dt)=="region_ID"] <- "nTADs"
  tmp_dt <- aggregate(as.formula(paste0("region_ID~region_cptmt", suffix)), data=all_ds_dt, FUN=length)
  nTot <- setNames(tmp_dt$region_ID, tmp_dt[, paste0("region_cptmt", suffix)])
  agg_dt$ratioTADs <- agg_dt$nTADs/nTot[paste0(agg_dt[,paste0("region_cptmt", suffix)])]
  plotTit <- paste0("ratio TADs by signif. cptmt", suffix, " (all DS)")
  ggbar_p <-  ggbarplot(agg_dt, 
                        y="ratioTADs",
                        x=paste0("region_cptmt", suffix),
                        fill=paste0("tadSignif")) +
    ggtitle(plotTit, subtitle=mySub) +
    mytheme +
    labs(x="" , y ="ratio of TADs", color=paste0(legTitle),fill=paste0(legTitle)) + 
    theme(
      axis.text.x = element_text(hjust=1, vjust=0.5,size=10,angle=90)
    )
  outFile <- file.path(outFolder,paste0("signif_cptmt", suffix, "_allDS_barplot.", plotType))
  ggsave(ggbar_p, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
}








