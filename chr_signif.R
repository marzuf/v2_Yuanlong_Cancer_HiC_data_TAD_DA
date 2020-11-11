result_dt <- get(load("DATA_TMP/all_result_dt.Rdata"))

result_dt$chromo <- gsub("(chr.+)_.+", "\\1", result_dt$region)

signif_thresh <- 0.01

agg_dt <- aggregate(adjPvalComb~hicds + exprds + chromo, data=result_dt, FUN =function(x)sum(x<= signif_thresh))

agg_dt$chromo <- factor(agg_dt$chromo, levels=paste0("chr", 1:22))
stopifnot(!is.na(agg_dt$chromo))

ggboxplot(agg_dt, y = paste0("adjPvalComb"), x ="chromo", add="jitter",
          xlab ="") + 
  ggtitle(paste0(""), subtitle=paste0("")) + 
  theme(
    axis.text.x=element_text(angle=90)
  )
