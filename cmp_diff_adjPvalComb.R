
# Rscript cmp_diff_adjPvalComb.R

outFolder <- "CMP_DIFF_ADJPVALCOMB"
dir.create(outFolder, recursive = TRUE)
plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.4
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.4


v0_dt <- get(load(file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")))
v0_rank_dt <- do.call(rbind, by(v0_dt, list(v0_dt$hicds, v0_dt$exprds), function(sub_dt) {
  sub_dt$tadRank <- rank(sub_dt$adjPvalComb, ties="min")
  sub_dt
}))

bothSameNbr_dt <- get(load(file.path("CREATE_FINAL_TABLE_BOTHSAMENBR/all_result_dt.Rdata")))
bothSameNbr_rank_dt <- do.call(rbind, by(bothSameNbr_dt, list(bothSameNbr_dt$hicds, bothSameNbr_dt$exprds), function(sub_dt) {
  sub_dt$tadRank <- rank(sub_dt$adjPvalComb, ties="min")
  sub_dt
}))

bothSameNbrDouble_dt <- get(load(file.path("CREATE_FINAL_TABLE_BOTHSAMENBRDOUBLE/all_result_dt.Rdata")))
bothSameNbrDouble_rank_dt <- do.call(rbind, by(bothSameNbrDouble_dt, list(bothSameNbrDouble_dt$hicds, bothSameNbrDouble_dt$exprds), function(sub_dt) {
  sub_dt$tadRank <- rank(sub_dt$adjPvalComb, ties="min")
  sub_dt
}))


bothRandom_dt <- get(load(file.path("CREATE_FINAL_TABLE_BOTHRANDOM//all_result_dt.Rdata")))
bothRandom_rank_dt <- do.call(rbind, by(bothRandom_dt, list(bothRandom_dt$hicds, bothRandom_dt$exprds), function(sub_dt) {
  sub_dt$tadRank <- rank(sub_dt$adjPvalComb, ties="min")
  sub_dt
}))

bothRandomRescaled_dt <- get(load(file.path("CREATE_FINAL_TABLE_BOTHRANDOMRESCALED//all_result_dt.Rdata")))
bothRandomRescaled_rank_dt <- do.call(rbind, by(bothRandomRescaled_dt, list(bothRandomRescaled_dt$hicds, bothRandomRescaled_dt$exprds), function(sub_dt) {
  sub_dt$tadRank <- rank(sub_dt$adjPvalComb, ties="min")
  sub_dt
}))

merged_dt <- merge(merge(v0_rank_dt[,c("hicds", "exprds", "region", "adjPvalComb", "tadRank")], 
                   bothSameNbr_rank_dt[,c("hicds", "exprds", "region", "adjPvalComb", "tadRank")], 
                   by=c("hicds", "exprds", "region"), all=TRUE, suffixes = c("_v0", "_bothSameNbr")),
                   bothRandom_rank_dt[,c("hicds", "exprds", "region", "adjPvalComb", "tadRank")],
                   by=c("hicds", "exprds", "region"), all=TRUE)
colnames(merged_dt)[colnames(merged_dt) == "adjPvalComb"] <- "adjPvalComb_bothRandom"
colnames(merged_dt)[colnames(merged_dt) == "tadRank"] <- "tadRank_bothRandom"

merged_dt <- merge( merged_dt,bothRandomRescaled_rank_dt[,c("hicds", "exprds", "region", "adjPvalComb", "tadRank")],
                   by=c("hicds", "exprds", "region"), all=TRUE)
colnames(merged_dt)[colnames(merged_dt) == "adjPvalComb"] <- "adjPvalComb_bothRandomRescaled"
colnames(merged_dt)[colnames(merged_dt) == "tadRank"] <- "tadRank_bothRandomRescaled"

merged_dt <- merge( merged_dt,bothSameNbrDouble_rank_dt[,c("hicds", "exprds", "region", "adjPvalComb", "tadRank")],
                   by=c("hicds", "exprds", "region"), all=TRUE)
colnames(merged_dt)[colnames(merged_dt) == "adjPvalComb"] <- "adjPvalComb_bothSameNbrDouble"
colnames(merged_dt)[colnames(merged_dt) == "tadRank"] <- "tadRank_bothSameNbrDouble"



stopifnot(!is.na(merged_dt))

tadSignifThresh <- 0.05




agg_v0_dt <- aggregate(adjPvalComb_v0 ~ hicds+exprds, data=merged_dt, function(x) sum(x <= tadSignifThresh))
agg_random_dt <- aggregate(adjPvalComb_bothRandom ~ hicds+exprds, data=merged_dt, function(x) sum(x <= tadSignifThresh))
agg_sameNbr_dt <- aggregate(adjPvalComb_bothSameNbr ~ hicds+exprds, data=merged_dt, function(x) sum(x <= tadSignifThresh))

nMerged_dt <- merge(agg_sameNbr_dt, merge(agg_v0_dt, agg_random_dt, all=TRUE, by=c("hicds", "exprds")), all=TRUE, by=c("hicds", "exprds"))

agg_randomRescaled_dt <- aggregate(adjPvalComb_bothRandomRescaled ~ hicds+exprds, data=merged_dt, function(x) sum(x <= tadSignifThresh))
nMerged_dt <- merge(nMerged_dt,agg_randomRescaled_dt, all=TRUE, by=c("hicds", "exprds"))

agg_sameNbrDouble_dt <- aggregate(adjPvalComb_bothSameNbrDouble ~ hicds+exprds, data=merged_dt, function(x) sum(x <= tadSignifThresh))
nMerged_dt <- merge(nMerged_dt,agg_sameNbrDouble_dt, all=TRUE, by=c("hicds", "exprds"))



nMerged_dt$labcols <- all_cols[all_cmps[nMerged_dt$exprds]]

stopifnot(!is.na(nMerged_dt))

save(nMerged_dt, file ="nMerged_dt.Rdata", version=2)

for(patt in c("bothSameNbrDouble")) {
#for(patt in c("bothRandom", "bothSameNbr", "bothRandomRescaled")) {
  
  my_x <- nMerged_dt$adjPvalComb_v0
  my_xlab <- "# signif v0"
  my_y <- nMerged_dt[,paste0("adjPvalComb_", patt)]
  my_ylab <- paste0("# signif ", patt)
  
  outFile <- file.path(outFolder, paste0("nbr_", patt, "_vs_nbr_v0_signif_", tadSignifThresh, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x=my_x,
    y=my_y,
    xlab=my_xlab,
    ylab=my_ylab,
    pch=16,
    col=nMerged_dt$labcols,
    main=paste0("# signif. ", patt, " vs. v0"),
    cex.axis=plotCex,
    cex.lab=plotCex,
    cex.main=plotCex
  )
  mtext(side=3, text=paste0("all DS - n=",nrow(nMerged_dt), "; p-val <=", tadSignifThresh))
  # addCorr(x=my_x, y  = my_y, bty="n", legPos = "topleft")
  legend("topleft",
         legend=names(all_cols),
         pch=16,
         col = all_cols, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}





i_col=1
j_col=2

nDS <- length(unique(file.path(merged_dt$hicds, merged_dt$exprds)))

init_merged_dt <- merged_dt
merged_dt <- merged_dt[,grepl("adjPval", colnames(merged_dt))]

for(i_col in 1:(ncol(merged_dt)-1)) {
  
  my_x <- -log10(merged_dt[,i_col])
  my_xlab <- paste0(colnames(merged_dt)[i_col])
  cmp_x <- gsub("adjPvalComb_", "", my_xlab)
  my_xlab <- paste0(my_xlab, " [-log10]")
  
  for(j_col in (i_col+1):ncol(merged_dt)){
    
    my_y <- -log10(merged_dt[,j_col])
    my_ylab <- paste0(colnames(merged_dt)[j_col])
    cmp_y <- gsub("adjPvalComb_", "", my_ylab)
    my_ylab <- paste0(my_ylab, " [-log10]")
    
    outFile <- file.path(outFolder, paste0(cmp_y, "_vs_", cmp_x, "_adjPvalComb_log10.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x=my_x,
      y=my_y,
      xlab=my_xlab,
      ylab=my_ylab,
      main=paste0(cmp_y, " vs. ", cmp_x, " adjPvalComb"),
      cex.axis=plotCex,
      cex.lab=plotCex,
      cex.main=plotCex
    )
    curve(1*x, add=TRUE, lty=2, col="grey")
    mtext(side=3, text=paste0("all DS - n=", nDS, "; # TADs=", nrow(merged_dt)))
    addCorr(x=my_x, y  = my_y, bty="n", legPos = "topleft")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
  }
  
  
}

merged_dt <- init_merged_dt[,grepl("tadRank", colnames(init_merged_dt))]

for(i_col in 1:(ncol(merged_dt)-1)) {
  
  my_x <- log10(merged_dt[,i_col])
  my_xlab <- paste0(colnames(merged_dt)[i_col])
  cmp_x <- gsub("tadRank_", "", my_xlab)
  my_xlab <- paste0(my_xlab, " [log10]")
  
  for(j_col in (i_col+1):ncol(merged_dt)){
    
    my_y <- log10(merged_dt[,j_col])
    my_ylab <- paste0(colnames(merged_dt)[j_col])
    cmp_y <- gsub("tadRank_", "", my_ylab)
    my_ylab <- paste0(my_ylab, " [log10]")
    
    outFile <- file.path(outFolder, paste0(cmp_y, "_vs_", cmp_x, "_tadRank_log10.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x=my_x,
      y=my_y,
      xlab=my_xlab,
      ylab=my_ylab,
      main=paste0(cmp_y, " vs. ", cmp_x, " tadRank"),
      cex.axis=plotCex,
      cex.lab=plotCex,
      cex.main=plotCex
    )
    curve(1*x, add=TRUE, lty=2, col="grey")
    mtext(side=3, text=paste0("all DS - n=", nDS, "; # TADs=", nrow(merged_dt)))
    addCorr(x=my_x, y  = my_y, bty="n", legPos = "topleft")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
  }
  
  
}
for(i_col in 1:(ncol(merged_dt)-1)) {
  
  my_x <- (merged_dt[,i_col])
  my_xlab <- paste0(colnames(merged_dt)[i_col])
  cmp_x <- gsub("tadRank_", "", my_xlab)
  my_xlab <- paste0(my_xlab, "")
  
  for(j_col in (i_col+1):ncol(merged_dt)){
    
    my_y <- (merged_dt[,j_col])
    my_ylab <- paste0(colnames(merged_dt)[j_col])
    cmp_y <- gsub("tadRank_", "", my_ylab)
    my_ylab <- paste0(my_ylab, "")
    
    outFile <- file.path(outFolder, paste0(cmp_y, "_vs_", cmp_x, "_tadRank.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x=my_x,
      y=my_y,
      xlab=my_xlab,
      ylab=my_ylab,
      main=paste0(cmp_y, " vs. ", cmp_x, " tadRank"),
      cex.axis=plotCex,
      cex.lab=plotCex,
      cex.main=plotCex
    )
    curve(1*x, add=TRUE, lty=2, col="grey")
    mtext(side=3, text=paste0("all DS - n=", nDS, "; # TADs=", nrow(merged_dt)))
    addCorr(x=my_x, y  = my_y, bty="n", legPos = "topleft")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
  }
  
  
}

