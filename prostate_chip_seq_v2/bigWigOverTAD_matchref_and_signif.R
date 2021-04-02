plotType <- "png"
myHeightGG <- 6
myWidthGG <- 7

outFolder <- "BIGWIGOVERTAD_MATCHREF_AND_SIGNIF"
dir.create(outFolder, recursive = TRUE)

# Rscript bigWigOverTAD_matchref_and_signif.R

require(ggsci)
require(ggpubr)
require(ggplot2)

pvalthresh <- 0.05
runFolder <- ".."
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
resultData$region_id <- file.path(resultData$dataset, resultData$region)
resultData$signif_lab <- ifelse(resultData$adjPvalComb <= pvalthresh, "signif", "notsignif")
resultData$direction_lab <- ifelse(resultData$meanLogFC <0, "down", 
                                   ifelse(resultData$meanLogFC >0, "up",NA))
stopifnot(!is.na(resultData$direction_lab))
resultData$signif_lab <- paste0(resultData$signif_lab, ".", resultData$direction_lab)
signif_labs <- setNames(resultData$signif_lab, resultData$region_id)

#   name - name field from bed, which should be unique
#   size - size of bed (sum of exon sizes
#   covered - # bases within exons covered by bigWig
#   sum - sum of values over all bases covered
#   mean0 - average over bases with non-covered bases counting as zeroes
#   mean - average over just covered bases

exprds <- "TCGAprad_norm_prad"


all_files <- c("chip_22Rv1_cover_22Rv1_TADs.bed",
               "chip_RWPE1_cover_RWPE1_TADs.bed")

for(infile in all_files) {
  
  cat(paste0("... start ", infile, "\n"))
  
  
  hicds <- ifelse(grepl("22Rv1", infile), "GSE118514_22Rv1_40kb",
                  ifelse(grepl("RWPE1", infile), "GSE118514_RWPE1_40kb", NA))
  stopifnot(!is.na(hicds))
  
  match_hicds <- ifelse(grepl("RWPE1", infile), "GSE118514_22Rv1_40kb",
                  ifelse(grepl("22Rv1", infile), "GSE118514_RWPE1_40kb", NA))
  stopifnot(!is.na(match_hicds))
  
  stopifnot(hicds != match_hicds)
  cat("match_hicds = ", match_hicds,"\n")
  
  match_hicds_lab <- gsub(".+_(.+)_40kb","\\1", match_hicds)
  
  ref_hicds_lab <- gsub(".+_(.+)_40kb","\\1", hicds)
  
  
  cat("match_hicds_lab = ", match_hicds_lab,"\n")
  
  cat("ref_hicds_lab = ", ref_hicds_lab,"\n")
  
  
  # ref is 22rv1 -> i want rwpe1 over 22rv1
  matchfile <- ifelse(grepl("RWPE1", match_hicds), "chip_RWPE1_cover_22Rv1_TADs.bed",
                      ifelse(grepl("22Rv1", match_hicds), "chip_22Rv1_cover_RWPE1_TADs.bed", NA))
                     
  stopifnot(!is.na(matchfile))
  
  
  cat("ref file ", infile,"\n")
  cat("matchfile  ", matchfile,"\n")
  
  
  bw_over_tad_dt <- read.delim(infile,
                               header=F,
                               col.names=c("region", "size_ref", "covered_ref", "sum_ref", "mean0_ref", "mean_ref"),
                               stringsAsFactors = FALSE)
  # bw_over_tad_dt$region_id <- file.path(hicds, exprds, bw_over_tad_dt$region)
  
  
  match_bw_over_tad_dt <- read.delim(matchfile,
                               header=F,
                               col.names=c("region", "size_match", "covered_match", "sum_match", "mean0_match", "mean_match"),
                               stringsAsFactors = FALSE)
  # match_bw_over_tad_dt$region_id <- file.path(hicds, exprds, match_bw_over_tad_dt$region)
  
  
  cat("nrow bw ", nrow(bw_over_tad_dt), "\n")
  cat("nrow match_bw_over_tad_dt ", nrow(match_bw_over_tad_dt), "\n")
  
  stopifnot(setequal(bw_over_tad_dt$region, match_bw_over_tad_dt$region))
  
  all_over_tad_dt <- merge(bw_over_tad_dt, match_bw_over_tad_dt,
                           by =c("region"), all=FALSE)
  
  stopifnot(all_over_tad_dt$size_ref == all_over_tad_dt$size_match)
  
  cat("nrow all_over_tad_dt ", nrow(all_over_tad_dt), "\n")
  
  
  # stopifnot(nrow(all_over_tad_dt) == nrow(match_bw_over_tad_dt))
  # stopifnot(nrow(all_over_tad_dt) == nrow(bw_over_tad_dt))
  stopifnot(!is.na(all_over_tad_dt))
  
  all_over_tad_dt$region_id <- file.path(hicds, exprds, all_over_tad_dt$region)
  
  
  sub_result_dt <- resultData[resultData$hicds == hicds & resultData$exprds == exprds,]
  
  
  merge_dt <- merge(sub_result_dt, all_over_tad_dt, by="region_id", all.x=T, all.y=F)
  stopifnot(!is.na(merge_dt))
  
  plotcols <- c( "sum", "mean0", "mean")
  # stopifnot(plotcols %in% colnames(bw_over_tad_dt))
  
  merge_dt$signif_lab2 <- gsub("notsignif.+", "notsignif", merge_dt$signif_lab)
  
  save(merge_dt, file="merge_dt.Rdata", version=2)
  
  for(toplot in plotcols) {
    
    cat(paste0("... start ", infile, " - ", toplot, "\n"))
    
    newcol <- paste0(toplot, "_ratio")
    
    merge_dt[,newcol] <-  merge_dt[,paste0(toplot, "_match")]/ merge_dt[,paste0(toplot, "_ref")]
    
    mysub2 <- paste0(names(table(merge_dt$signif_lab2)), "=", as.numeric(table(merge_dt$signif_lab2)), 
                     collapse=";")
    
    plotTit <- gsub("\\.bed", "", infile)
    legTitle <- ""
    mySub <- paste0(toplot, " - ", mysub2)
    
    p3 <-  ggboxplot(merge_dt,
                     y = paste0(newcol),
                     outlier.shape=NA,
                     add="jitter",
                     rug = FALSE,  
                     ylab = paste0(toplot, ": ", match_hicds_lab, "/", ref_hicds_lab),
                     xlab ="",
                     x="signif_lab2",                   # fill = paste0("tad_signif"),
                     color = "signif_lab2",
                     palette = "jco") +
      ggtitle(plotTit, subtitle = mySub) +
      labs(color=paste0(legTitle),fill=paste0(legTitle)) +
      guides(color=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    
    
    outFile <- file.path(outFolder, paste0(plotTit, "_", toplot, "_match_", match_hicds_lab,
                                           "_over_ref_", ref_hicds_lab, "_boxplot.", plotType))
    ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    merge_dt[,"newcol_log2"] <-  log2(merge_dt[,paste0(newcol)])
    
    p3 <-  ggboxplot(merge_dt,
                     y = paste0("newcol_log2"),
                     outlier.shape=NA,
                     add="jitter",
                     rug = FALSE,  
                     ylab = paste0(toplot, ": ", match_hicds_lab, "/", ref_hicds_lab, " [log2]"),
                     xlab ="",
                     x="signif_lab2",                   # fill = paste0("tad_signif"),
                     color = "signif_lab2",
                     palette = "jco") +
      ggtitle(plotTit, subtitle = mySub) +
      labs(color=paste0(legTitle),fill=paste0(legTitle)) +
      guides(color=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    
    
    outFile <- file.path(outFolder, paste0(plotTit, "_", toplot, "_match_", match_hicds_lab,
                                           "_over_ref_", ref_hicds_lab, "_log2_boxplot.", plotType))
    ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
    
  }
  
  
}
















# # take the 22rv1 tad cover by 22rv1
# 
# bw_over_tad_dt <- read.delim("chip_22Rv1_cover_22Rv1_TADs.bed", 
#                              header=F,
#                              col.names=c("region", "size", "covered", "sum", "mean0", "mean"),
#                              stringsAsFactors = FALSE)
# bw_over_tad_dt$region_id <- file.path(hicds, exprds, bw_over_tad_dt$region)
# 
# 
# # take the 22rv1 tad cover by rwpe1
# 
# merge on region id
# take the ratio of the coverage 
# 
# # -> compare if the methylation is different between normal and cancer for the signif ones
# 
# 
# ------------------------------------------------------------------------------------------------------------
#   bigWigAverageOverBed v2 - Compute average score of big wig over each bed, which may have introns.
# usage:
#   bigWigAverageOverBed in.bw in.bed out.tab
# The output columns are:
#   name - name field from bed, which should be unique
# size - size of bed (sum of exon sizes
#                     covered - # bases within exons covered by bigWig
#                       sum - sum of values over all bases covered
#                     mean0 - average over bases with non-covered bases counting as zeroes
#                     mean - average over just covered bases
#                     Options:
#                       -stats=stats.ra - Output a collection of overall statistics to stat.ra file
#                     -bedOut=out.bed - Make output bed that is echo of input bed but with mean column appended
#                     -sampleAroundCenter=N - Take sample at region N bases wide centered around bed item, rather
#                     than the usual sample in the bed item.
#                     -minMax - include two additional columns containing the min and max observed in the area.