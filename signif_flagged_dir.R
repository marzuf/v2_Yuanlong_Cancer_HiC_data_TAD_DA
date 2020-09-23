require(ggplot2)
require(reshape2)
require(ggsci)

# Rscript signif_flagged_dir.R

outFolder <- "SIGNIF_FLAGGED_DIR"
dir.create(outFolder, recursive = TRUE)

plotType <- "svg"
myWidth <- 8
myHeight <- 5

curr_hicds <- "ENCSR489OCU_NCI-H460_40kb"
curr_exprds <- "TCGAluad_mutKRAS_mutEGFR"  

puritycorr_dt <- get(load("ALLTADS_AND_PURITY_FINAL/aran/CPE/log10/all_ds_corrPurity_dt.Rdata"))
agg_pc_dt <- aggregate(purityCorr~dataset+region, data=puritycorr_dt, FUN=mean)

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
final_dt$dataset <- file.path(final_dt$hicds, final_dt$exprds)

signifThresh <- 0.01

merge_dt <- merge(final_dt, agg_pc_dt, by=c("dataset", "region"), all=F)

merge_dt$signif_tad <- merge_dt$adjPvalComb <= signifThresh
merge_dt$logFC_sign <- sign(merge_dt$meanLogFC)

corrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif_tad], probs = 0.05))

merge_dt$purity_flagged <- merge_dt$purityCorr <= corrThresh

signif_dt <- merge_dt[merge_dt$signif_tad,]

tmp <- aggregate(purity_flagged~dataset+logFC_sign, FUN=sum, data=signif_dt)
tmp[tmp$dataset == "ENCSR489OCU_NCI-H460_40kb/TCGAluad_mutKRAS_mutEGFR",]

all_ds_dt <- do.call(rbind, by(signif_dt, signif_dt$dataset, function(sub_dt) {
  tmp <- aggregate(purity_flagged~dataset+logFC_sign, FUN=sum, data=sub_dt)
  tmp$purity_flagged_ratio <- tmp$purity_flagged/sum(tmp$purity_flagged)
  tmp
}))
rownames(all_ds_dt) <- NULL

order_dt <- aggregate(purity_flagged~dataset, data=all_ds_dt, FUN=sum)
order_dt <- order_dt[order(order_dt$purity_flagged, decreasing = TRUE),]

plot_dt <- melt(all_ds_dt, id=c("dataset", "logFC_sign"))

plot_dt$dataset <- factor(plot_dt$dataset, levels = order_dt$dataset)

plot_dt$logFC_lab <- ifelse(plot_dt$logFC_sign == 1, "higher expr. cond2",
                            ifelse(plot_dt$logFC_sign == -1,"higher expr. cond1", NA ))
stopifnot(!is.na(plot_dt$logFC_lab))

add_theme <- function(p) {
  p+ 
    scale_fill_aaas()+
    # scale_y_continuous(expand=c(0,0))+
    theme(
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_blank(),#text(size=12, hjust=0.5, vjust=0.5),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold"),
      legend.text=element_text(size=11)
    )
}

plotTit <- "Purity-flagged signif. TADs by cond. of higher expr."
subTit <- paste0("# DS ", length(unique(plot_dt$dataset)))

nbr_p <- ggplot(plot_dt[plot_dt$variable=="purity_flagged",], aes(x=dataset, y=value, fill=logFC_lab))+
  ggtitle(plotTit, subtitle=subTit)+
  labs(x="all datasets", y="# purity-flagged signif. TADs", fill="")+
  geom_bar(stat="identity", position="stack") 
nbr_p <- add_theme(nbr_p)
outFile <- file.path(outFolder, paste0("nbr_purityFlagged_signif_by_cond_higherExpr_barplot.", plotType))
ggsave(nbr_p, filename = outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile,"\n"))

ratio_p <- ggplot(plot_dt[plot_dt$variable=="purity_flagged_ratio",], aes(x=dataset, y=value, fill=logFC_lab))+
  ggtitle(plotTit, subtitle=subTit)+
  labs(x="all datasets", y="ratio purity-flagged signif. TADs", fill="")+
  geom_bar(stat="identity", position="stack")
ratio_p <- add_theme(ratio_p)
outFile <- file.path(outFolder, paste0("ratio_purityFlagged_signif_by_cond_higherExpr_barplot.", plotType))
ggsave(ratio_p, filename = outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile,"\n"))

purity_byCond_data <- get(load("PURITY_BY_COND_FINAL/aran/CPE/log10/all_ds_corrPurity_data.Rdata"))
purity_level_dt <- do.call(rbind, lapply(purity_byCond_data, function(x)x[["ds_purity_dt"]]))
purity_level_dt$cond1 <- gsub("TCGA.+_(.+)_.+", "\\1", basename(purity_level_dt$dataset))
purity_level_dt$cond2 <- gsub("TCGA.+_.+_(.+)", "\\1", basename(purity_level_dt$dataset))
purity_level_dt$cond_id <- purity_level_dt$cond
purity_level_dt$cond <- ifelse(purity_level_dt$cond_id == purity_level_dt$cond1, "cond1", 
                               ifelse(purity_level_dt$cond_id == purity_level_dt$cond2, "cond2", NA))
stopifnot(!is.na(purity_level_dt$cond))
agg_cond_pl_dt <- aggregate(purity~cond+dataset, FUN=mean, data=purity_level_dt)

ratio_cond_pl_dt <- do.call(rbind, by(agg_cond_pl_dt, agg_cond_pl_dt$dataset, function(sub_dt){
  sub_dt$purity[sub_dt$cond == "cond1"]/sub_dt$purity[sub_dt$cond == "cond2"]
}))
ratio_cond_pl_dt <- data.frame(ratio_cond_pl_dt)
ratio_cond_pl_dt$dataset <- rownames(ratio_cond_pl_dt)
rownames(ratio_cond_pl_dt) <- NULL


# logFC < 0 => higher expression in cond1
# logFC > 0 => higher expression in cond2

ratio_cond_signif_flagged_dt <- all_ds_dt[all_ds_dt$logFC_sign == -1, c("dataset", "purity_flagged_ratio")]

ratio_dt <- merge(ratio_cond_pl_dt, ratio_cond_signif_flagged_dt, by=c("dataset"))

outFile <- file.path(outFolder, paste0("ratio_purityLevel_purityFlagged_dotplot.", plotType))
do.call(plotType,list(outFile, height=myHeight, width=myHeight))
plot(
  x=ratio_dt$ratio_cond_pl_dt,
  y=ratio_dt$purity_flagged_ratio,
  xlab="ratio cond1/cond2 mean purity level",
  ylab="ratio purity-flagged signif. TADs with higher expr. in cond1",
  pch=16,
  main= "Ratio purity level and ratio purity-flagged signif. TADs"
)
mtext(side=3, text=paste0("# DS = ", length(unique(ratio_dt$dataset))))
points(
  x=ratio_dt$ratio_cond_pl_dt[ratio_dt$dataset == file.path(curr_hicds, curr_exprds)],
  y=ratio_dt$purity_flagged_ratio[ratio_dt$dataset == file.path(curr_hicds, curr_exprds)],
  col="red",
  pch=16,
  cex=2
)
legend("topright",
       col="red", 
       pch=16, pt.cex=2,
       legend=paste0(curr_hicds, " -\n" , curr_exprds), 
       bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))

