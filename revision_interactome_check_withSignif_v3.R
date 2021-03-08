inFolder <- "REVISION_INTERACTOME"
ref_hicds <- "GSE118514_RWPE1_40kb"

plotType <- "png"
myWidth <- 400
myHeight <- 400

plotCex <- 1.2

# Rscript revision_interactome_check_withSignif_v3.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- "REVISION_INTERACTOME_CHECK_WITHSIGNIF_V3"
dir.create(outFolder, recursive=TRUE)

all_pairs_cols <- c("ENCSR346DCU_LNCaP_40kb" = "LNCaP",
                    "GSE118514_RWPE1_40kb" = "RWPE1",
                    "GSE118514_22Rv1_40kb"="22Rv1")

ref_lab <- as.character(all_pairs_cols[ref_hicds])

ref_dt <- get(load(file.path(inFolder, paste0(ref_hicds, "_all_interactome_dt.Rdata"))))

ref_dt$coord <- file.path(ref_dt$chromo, ref_dt$tad_start, ref_dt$tad_end)
stopifnot(!duplicated(ref_dt$coord))

tadSignifThresh <- 0.01
final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$region_ID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
final_table_DT <- final_table_DT[final_table_DT$hicds == ref_hicds & final_table_DT$exprds == "TCGAprad_norm_prad",]
stopifnot(nrow(final_table_DT) > 0)
final_table_DT$chromo <- gsub("(chr.+)_TAD.+", "\\1",final_table_DT$region  )
signif_coords <- file.path(final_table_DT$chromo[final_table_DT$adjPvalComb <= tadSignifThresh],
                           final_table_DT$start[final_table_DT$adjPvalComb <= tadSignifThresh],
                           final_table_DT$end[final_table_DT$adjPvalComb <= tadSignifThresh])


# REVISION_INTERACTOME/ENCSR346DCU_LNCaP_40kb_all_interactome_dt.Rdata
# REVISION_INTERACTOME/GSE118514_22Rv1_40kb_all_interactome_dt.Rdata
# REVISION_INTERACTOME/GSE118514_RWPE1_40kb_all_interactome_dt.Rdata

cmp_col1 <- as.character(all_pairs_cols[all_pairs_cols != ref_lab][1])
cmp_col2 <- as.character(all_pairs_cols[all_pairs_cols != ref_lab][2])

cmp_ds1 <-  as.character(names(all_pairs_cols)[all_pairs_cols != ref_lab][1])
cmp_ds2 <-  as.character(names(all_pairs_cols)[all_pairs_cols != ref_lab][2])

all_cols <- c("sig_ratio_diff", "mean_p_diff", paste0("nMatchingEntries_", ref_lab),
              paste0("nRatioMatchEntries_", ref_lab), paste0("nSameEntries_", ref_lab))
plotcol = all_cols[1]

ref_dt[,paste0("sig_ratio_diff_",cmp_col1 )] <- ref_dt[,paste0("sig_ratio_",cmp_col1 )]-
                                                      ref_dt[,paste0("sig_ratio_",ref_lab )]
  
ref_dt[,paste0("mean_p_diff_",cmp_col1 )] <- ref_dt[,paste0("mean_p_",cmp_col1 )]-
  ref_dt[,paste0("mean_p_",ref_lab )]

ref_dt[,paste0("sig_ratio_diff_",cmp_col2 )] <- ref_dt[,paste0("sig_ratio_",cmp_col2 )]-
  ref_dt[,paste0("sig_ratio_",ref_lab )]


ref_dt[,paste0("mean_p_diff_",cmp_col2 )] <- ref_dt[,paste0("mean_p_",cmp_col2)]-
  ref_dt[,paste0("mean_p_",ref_lab )]

dot_dt <- ref_dt[ref_dt$coord %in% signif_coords,]

for(plotcol in all_cols) {
  
  
  myxlab <- paste0(plotcol, " - ", cmp_col1)
  myylab <- paste0(plotcol, " - ", cmp_col2)
  
  plotTit <- paste0(plotcol, "; ref: ", ref_lab)
  
  mySub <- paste0("# TADs = ", nrow(ref_dt))
  
  myy <- ref_dt[,paste0(plotcol, "_", cmp_col1)]
  myx <- ref_dt[,paste0(plotcol, "_", cmp_col2)]
  
  points_x <- dot_dt[,paste0(plotcol, "_", cmp_col1)]
  points_y <- dot_dt[,paste0(plotcol, "_", cmp_col2)]
  
  outFile <- file.path(outFolder, paste0(ref_lab,  "_", plotcol, "_", cmp_col1, "_", cmp_col2, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(x=myx,
           y=myy,
           xlab = myxlab,
           ylab =myylab,
           main =plotTit,
           cex.main=plotCex,
           cex.lab=plotCex,
           cex.axis=plotCex)
  mtext(side=3, text = plotTit)
  addCorr(x=myx,y=myy, legPos="topleft", bty="n")
  points(x=points_x,y=points_y, col="darkgreen", size=5)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}
######################################################### CHECK
other1_dt <- get(load(file.path(inFolder, paste0(cmp_ds1, "_all_interactome_dt.Rdata"))))
other2_dt <- get(load(file.path(inFolder, paste0(cmp_ds2, "_all_interactome_dt.Rdata"))))

all_merge_dt <- merge(ref_dt, 
                   merge(other1_dt,other2_dt, by=c("chromo", "tad_start", "tad_end")),
                   by=c("chromo", "tad_start", "tad_end"))

stopifnot(all_merge_dt$mean_p_LNCaP == all_merge_dt$mean_p_LNCaP.x)
stopifnot(all_merge_dt$mean_p_LNCaP == all_merge_dt$mean_p_LNCaP.y)
stopifnot(all_merge_dt$sig_ratio_LNCaP == all_merge_dt$sig_ratio_LNCaP.x)
stopifnot(all_merge_dt$sig_ratio_LNCaP == all_merge_dt$sig_ratio_LNCaP.y)
stopifnot(all_merge_dt$nMatchingEntries_LNCaP_22Rv1 == 
            all_merge_dt$nMatchingEntries_LNCaP_RWPE1)
stopifnot(all_merge_dt$nSameEntries_RWPE1_22Rv1 == all_merge_dt$nSameEntries_22Rv1_RWPE1)
stopifnot(all_merge_dt$nSameEntries_RWPE1_LNCaP == all_merge_dt$nSameEntries_LNCaP_RWPE1)



dt1 = get(load("REVISION_CHANGES_CPTMTLABELS_ALLDS_ALLTADS/tad2cptmt_dt.Rdata"))
dt2 = get(load("REVISION_CHANGES_CPTMTLABELS_ALLDS/tad2cptmt_dt.Rdata"))
dt = merge(dt1, dt2, by=c("chromo", "region", "start", "end", "region_ID"))
stopifnot(dt$startCptmtScore.x ==dt$startCptmtScore.y )
stopifnot(dt$adjPvalComb.x ==dt$adjPvalComb.y )
stopifnot(dt$tad_eightCptmtLab.x ==dt$tad_eightCptmtLab.y )

