
outFolder <- "TOPTAD_TOPGENE"
dir.create(outFolder, recursive = TRUE)

# Rscript topTAD_topGene.R

require(foreach)
require(doMC)
registerDoMC(40)

require(reshape2)
require(ggpubr)
require(ggsci)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "svg"
myHeight <- 7
myWidth <- 12

plotCex  <- 1.2

signif_dt <- get(load(file.path("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))


signif_dt$ds <- file.path(signif_dt$hicds, signif_dt$exprds)

all_ds <- unique(signif_dt$ds)

buildTable <- F

pipFolder <- file.path("PIPELINE/OUTPUT_FOLDER")

if(buildTable) {
  
  all_withfc_dt <- foreach(ds =all_ds, .combine='rbind') %dopar% {
    
    sub_dt <- signif_dt[signif_dt$ds == ds,]
    
    # 0_prepGeneData/pipeline_geneList.Rdata 
    # 1_runGeneDE/DE_topTable.Rdata
    geneList <- get(load(file.path(pipFolder, dirname(ds), basename(ds), "0_prepGeneData/pipeline_geneList.Rdata")))
    geneDE_dt <- get(load(file.path(pipFolder, dirname(ds), basename(ds), "1_runGeneDE/DE_topTable.Rdata")))
    geneDE_dt$genes <- as.character(geneDE_dt$genes)
    stopifnot(names(geneList) %in% geneDE_dt$genes)
    geneDE_dt <- geneDE_dt[geneDE_dt$genes %in% names(geneList),]
    geneDE_dt$entrezID <- geneList[geneDE_dt$genes]
    stopifnot(!is.na(geneDE_dt$entrezID))
    stopifnot(setequal(geneList, geneDE_dt$entrezID))  
    stopifnot(setequal(geneList, sub_dt$entrezID))  
    gene2fc <- setNames(geneDE_dt$logFC, geneDE_dt$entrezID)
    stopifnot(setequal(names(gene2fc), sub_dt$entrezID))  
    sub_dt$entrezID <- as.character(sub_dt$entrezID)
    sub_dt$geneFC <- gene2fc[sub_dt$entrezID]
    stopifnot(!is.na(sub_dt$geneFC))
    sub_dt
  }
  
  
  outFile <- file.path(outFolder, "all_withfc_dt.Rdata")
  save(all_withfc_dt, file=outFile, version=2)
  cat(paste0("written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_withfc_dt.Rdata")
  all_withfc_dt <- get(load(outFile))
  
}

abs_fc_thresh <- 1
signif_tad_thresh <- 0.01
signif_gene_thresh <- 0.01

abs_fc_qtThresh <- 0.95

normQt <- 0.95

all_rel_dt <- do.call(rbind, by(all_withfc_dt, all_withfc_dt$ds, function(sub_dt) {
  
  sub_dt$tad_rank_rel <- sub_dt$tad_rank/max(sub_dt$tad_rank)
  sub_dt$gene_rank_rel <- sub_dt$gene_rank/max(sub_dt$gene_rank)
  
  sub_dt$gene_fcRank <- rank(-abs(sub_dt$logFC), ties.method="min")
  sub_dt$gene_pvalRank <- rank(sub_dt$adj.P.Val, ties.method="min")
  
  sub_dt$gene_fcRank_rel <- sub_dt$gene_fcRank/max(sub_dt$gene_fcRank)
  
  sub_dt$gene_absFC_zscore <- as.numeric(scale(abs(sub_dt$logFC), center=T, scale=T))
  
  genefc_qt <- as.numeric(quantile(abs(sub_dt$logFC), probs=normQt))
  
  sub_dt$abs_gene_fc_qtThresh <- abs(sub_dt$logFC)/genefc_qt
  sub_dt
}))
stopifnot(nrow(all_rel_dt) == nrow(all_withfc_dt))
stopifnot(all_rel_dt$gene_pvalRank == all_rel_dt$gene_rank)

tmp_dt <- signif_dt[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tmp_dt <- unique(tmp_dt)
tmp_dt2 <- aggregate(tad_adjCombPval~hicds+exprds, FUN=function(x) sum(x <= signif_tad_thresh), data = tmp_dt)
tmp_dt2 <- tmp_dt2[order(tmp_dt2$tad_adjCombPval, decreasing=T),]
tmp_dt2$ds <- paste0(tmp_dt2$hicds, "\n", tmp_dt2$exprds)


tad_dt <- all_rel_dt[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tad_dt <- unique(tad_dt)
sum(tad_dt$tad_adjCombPval<=0.01)

zMinAgg_dt <- aggregate(gene_absFC_zscore~hicds+exprds+region+tad_adjCombPval, data=all_rel_dt, FUN=min)
zMinAgg_dt$ds <- paste0(zMinAgg_dt$hicds,"\n",zMinAgg_dt$exprds)
zMinAgg_dt$ds <- factor(zMinAgg_dt$ds, levels=tmp_dt2$ds)

signif_zMinAgg_dt <- zMinAgg_dt[zMinAgg_dt$tad_adjCombPval <= signif_tad_thresh,]

pBox <- ggboxplot(signif_zMinAgg_dt, y = paste0("gene_absFC_zscore"), x ="ds", add="jitter",
                  xlab ="Datasets ranked by decreasing # signif. TADs",
                  ylab="min z-score TAD gene abs FC") + 
  ggtitle(paste0("min gene FC z-score in signif. TADs"), subtitle=paste0("TAD p-val. signif. thresh. <= ", signif_tad_thresh)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme(
    # text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
    
  )

outFile <- file.path(outFolder, paste0("minGeneAbsFCzscoreInSignifTADs_boxplot.", "svg"))
ggsave(pBox, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))


zMeanAgg_dt <- aggregate(gene_absFC_zscore~hicds+exprds+region+tad_adjCombPval, data=all_rel_dt, FUN=mean)
zMeanAgg_dt$ds <- paste0(zMeanAgg_dt$hicds,"\n",zMeanAgg_dt$exprds)
zMeanAgg_dt$ds <- factor(zMeanAgg_dt$ds, levels=tmp_dt2$ds)

signif_zMeanAgg_dt <- zMeanAgg_dt[zMeanAgg_dt$tad_adjCombPval <= signif_tad_thresh,]

pBox <- ggboxplot(signif_zMeanAgg_dt, y = paste0("gene_absFC_zscore"), x ="ds", add="jitter",
                  xlab ="Datasets ranked by decreasing # signif. TADs",
                  ylab="mean z-score TAD gene abs FC") + 
  ggtitle(paste0("mean gene FC z-score in signif. TADs"), subtitle=paste0("TAD p-val. signif. thresh. <= ", signif_tad_thresh)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme(
    # text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
    
  )

outFile <- file.path(outFolder, paste0("meanGeneAbsFCzscoreInSignifTADs_boxplot.", "svg"))
ggsave(pBox, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))

stop("ok\n")



outFile <- file.path(outFolder, "relGeneRankPval_relTADrank.png")
do.call("png", list(outFile, height=400, width=400))
densplot(y=all_rel_dt$gene_rank_rel, 
           x=all_rel_dt$tad_rank_rel,
         xlab = "rel. TAD rank",
         ylab="rel. gene rank (pval)",
         cex.lab=plotCex,
         cex.axis=plotCex,
         cex.main=plotCex,
         main=paste0("gene rank vs. TAD rank")
         )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "relGeneRankFC_relTADrank.png")
do.call("png", list(outFile, height=400, width=400))
densplot(y=all_rel_dt$gene_fcRank_rel,
         x=all_rel_dt$tad_rank_rel,
         xlab = "rel. TAD rank",
         ylab="rel. gene rank (abs FC)",
         cex.lab=plotCex,
         cex.axis=plotCex,
         cex.main=plotCex,
         main=paste0("gene rank vs. TAD rank")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "relGeneQtNorm_relTADrank.png")
do.call("png", list(outFile, height=400, width=400))
densplot(y=all_rel_dt$abs_gene_fc_qtThresh, 
         x=all_rel_dt$tad_rank_rel,
         xlab = "rel. TAD rank",
         ylab=paste0("gene FC (",normQt, "-qt norm)"),
         cex.lab=plotCex,
         cex.axis=plotCex,
         cex.main=plotCex,
         main=paste0("gene rank vs. TAD rank")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


agg_pvalRank_dt <- aggregate(gene_rank_rel ~ hicds + exprds + region+tad_rank_rel, data=all_rel_dt, FUN=min)
agg_fcRank_dt <- aggregate(gene_fcRank_rel ~ hicds + exprds + region+tad_rank_rel, data=all_rel_dt, FUN=min)
agg_maxPvalRank_dt <- aggregate(gene_rank_rel ~ hicds + exprds + region+tad_rank_rel, data=all_rel_dt, FUN=max)
agg_maxFcRank_dt <- aggregate(gene_fcRank_rel ~ hicds + exprds + region+tad_rank_rel, data=all_rel_dt, FUN=max)

agg_minAbsFC_dt <- aggregate(geneFC ~ hicds + exprds + region+tad_rank_rel, data=all_rel_dt, FUN=function(x) min(abs(x)))
agg_maxAbsFC_dt <- aggregate(geneFC ~ hicds + exprds + region+tad_rank_rel, data=all_rel_dt, FUN=function(x) max(abs(x)))

stopifnot(nrow(agg_pvalRank_dt) == nrow(tad_dt))
stopifnot(nrow(agg_fcRank_dt) == nrow(tad_dt))
stopifnot(nrow(agg_maxAbsFC_dt) == nrow(tad_dt))
stopifnot(nrow(agg_minAbsFC_dt) == nrow(tad_dt))

outFile <- file.path(outFolder, "minRelGeneRankFC_relTADrank.png")
do.call("png", list(outFile, height=400, width=400))
densplot(y=agg_fcRank_dt$gene_fcRank_rel,
         x=agg_fcRank_dt$tad_rank_rel,
         xlab = "rel. TAD rank",
         ylab="min gene rank (abs FC)",
         cex.lab=plotCex,
         cex.axis=plotCex,
         cex.main=plotCex,
         main=paste0("gene rank vs. TAD rank")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "maxRelGeneRankFC_relTADrank.png")
do.call("png", list(outFile, height=400, width=400))
densplot(y=agg_maxFcRank_dt$gene_fcRank_rel,
         x=agg_maxFcRank_dt$tad_rank_rel,
         xlab = "rel. TAD rank",
         ylab="max gene rank (abs FC)",
         cex.lab=plotCex,
         cex.axis=plotCex,
         cex.main=plotCex,
         main=paste0("gene rank vs. TAD rank")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, "minRelGeneRankPval_relTADrank.png")
do.call("png", list(outFile, height=400, width=400))
densplot(y=agg_pvalRank_dt$gene_rank_rel,
         x=agg_pvalRank_dt$tad_rank_rel,
         xlab = "rel. TAD rank",
         ylab="min gene rank (pval)",
         cex.lab=plotCex,
         cex.axis=plotCex,
         cex.main=plotCex,
         main=paste0("gene rank vs. TAD rank")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "maxRelGeneRankPval_relTADrank.png")
do.call("png", list(outFile, height=400, width=400))
densplot(y=agg_maxPvalRank_dt$gene_rank_rel,
         x=agg_maxPvalRank_dt$tad_rank_rel,
         xlab = "rel. TAD rank",
         ylab="max gene rank (pval)",
         cex.lab=plotCex,
         cex.axis=plotCex,
         cex.main=plotCex,
         main=paste0("gene rank vs. TAD rank")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, "minGeneAbsFC_relTADrank.png")
do.call("png", list(outFile, height=400, width=400))
densplot(y=agg_minAbsFC_dt$geneFC,
         x=agg_minAbsFC_dt$tad_rank_rel,
         xlab = "rel. TAD rank",
         ylab="min gene abs FC",
         cex.lab=plotCex,
         cex.axis=plotCex,
         cex.main=plotCex,
         main=paste0("min gene FC vs. TAD rank")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "maxGeneAbsFC_relTADrank.png")
do.call("png", list(outFile, height=400, width=400))
densplot(y=agg_maxAbsFC_dt$geneFC,
         x=agg_maxAbsFC_dt$tad_rank_rel,
         xlab = "rel. TAD rank",
         ylab="max gene abs FC",
         cex.lab=plotCex,
         cex.axis=plotCex,
         cex.main=plotCex,
         main=paste0("max gene FC vs. TAD rank")
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

stop("-ok \n")



signifOnly_dt <- signif_dt[signif_dt$tad_adjCombPval <= signif_tad_thresh,]

minRankBySignifTAD_dt <- aggregate(gene_rank ~hicds+exprds+region, data=signifOnly_dt, min)
minRankBySignifTAD_dt$ds <- paste0(minRankBySignifTAD_dt$hicds,"\n", minRankBySignifTAD_dt$exprds)

minRankBySignifTAD_dt$ds <- factor(minRankBySignifTAD_dt$ds, levels = as.character(tmp_dt2$ds))
stopifnot(!is.na(minRankBySignifTAD_dt$ds))

minRankBySignifTAD_dt$min_gene_rank_log10 <- log10(minRankBySignifTAD_dt$gene_rank)

pBox <- ggboxplot(minRankBySignifTAD_dt, y = paste0("gene_rank"), x ="ds", add="jitter",
                  # xlab ="", 
                  xlab ="Datasets ranked by decreasing # signif. TADs",
                  ylab="min gene rank") + 
  ggtitle(paste0("min gene rank in signif. TADs"), subtitle=paste0("TAD p-val. signif. thresh. <= ", signif_tad_thresh)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme(
    # text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("minGeneRankInSignifTADs_boxplot.", "svg"))
ggsave(pBox, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))

pBox <- ggboxplot(minRankBySignifTAD_dt, y = paste0("min_gene_rank_log10"), x ="ds", add="jitter",
                  xlab ="Datasets ranked by decreasing # signif. TADs",
                  ylab="min gene rank [log10]") + 
  ggtitle(paste0("min gene rank in signif. TADs"), subtitle=paste0("TAD p-val. signif. thresh. <= ", signif_tad_thresh)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme(
    # text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
    
  )

outFile <- file.path(outFolder, paste0("minGeneRankInSignifTADs_log10_boxplot.", "svg"))
ggsave(pBox, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))



stop("-ok \n")



all_withfc_dt <- do.call(rbind, by(all_withfc_dt, all_withfc_dt$ds, function(sub_dt){
  fc_qt <- quantile(abs(sub_dt$geneFC), probs=abs_fc_qtThresh)
  sub_dt$highQt_geneFC <-  abs(sub_dt$logFC) >= fc_qt
  sub_dt
}))

fcQt_dt <- do.call(rbind, by(all_withfc_dt, all_withfc_dt$ds, function(sub_dt){
  fc_qt <- quantile(abs(sub_dt$geneFC), probs=abs_fc_qtThresh)
  data.frame(
    ds = unique(sub_dt$ds), 
    fc_qt = fc_qt,
    stringsAsFactors = FALSE
  )
}))
plotTit <- paste0("thresh. abs. FC ", abs_fc_qtThresh)
mySub <- paste0("# DS = ", length(fcQt_dt$ds))
legTitle <- ""

fcQt_dt <- fcQt_dt[order(fcQt_dt$fc_qt, decreasing = TRUE),]
fcQt_dt$ds_lab <- paste0(dirname(fcQt_dt$ds), "\n", basename(fcQt_dt$ds))
fcQt_dt$ds_lab <- factor(fcQt_dt$ds_lab, levels=fcQt_dt$ds_lab)
  
pQt <- ggbarplot(fcQt_dt,
                y = "fc_qt",
                x = "ds_lab",
                # y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                ylab = paste0("abs. FC ", abs_fc_qtThresh),
                xlab ="all datasets",
                # add = "median",                  # Add median line.
                rug = FALSE,                      # Add marginal rug
                color = "grey",
                fill = "grey",
                palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  # scale_color_manual(values=my_cols)+
  # scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand=c(0,0))+
  # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(
    # text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_text(size=6, angle=90, hjust=1, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  ) 

outFile <- file.path(outFolder, paste0("abs_fc_qt_thresh.", "svg"))
ggsave(pQt, file=outFile, height=myHeight*1.6, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))
# stop("-ok")

all_withfc_dt$signif_tad <- all_withfc_dt$tad_adjCombPval <= signif_tad_thresh
all_withfc_dt$signif_gene <- all_withfc_dt$adj.P.Val <= signif_gene_thresh
all_withfc_dt$high_geneFC <- abs(all_withfc_dt$logFC) >= abs_fc_thresh

any_highqtfc_gene_dt <- aggregate(highQt_geneFC ~ hicds + exprds + region, data=all_withfc_dt, any)
colnames(any_highqtfc_gene_dt)[colnames(any_highqtfc_gene_dt) == "highQt_geneFC"] <- "any_highQt_geneFC"


any_signif_gene_dt <- aggregate(signif_gene ~ hicds + exprds + region, data=all_withfc_dt, any)
colnames(any_signif_gene_dt)[colnames(any_signif_gene_dt) == "signif_gene"] <- "any_signif_gene"

any_highfc_gene_dt <- aggregate(high_geneFC ~ hicds + exprds + region, data=all_withfc_dt, any)
colnames(any_highfc_gene_dt)[colnames(any_highfc_gene_dt) == "high_geneFC"] <- "any_high_geneFC"

signif_dt <- all_withfc_dt[,c("hicds", "exprds", "region", "signif_tad")]
signif_dt <- unique(signif_dt)
stopifnot(!duplicated(file.path(signif_dt$hicds, signif_dt$exprds, signif_dt$region)))

all_dt <- merge(merge(merge(signif_dt, any_signif_gene_dt, by=c("hicds", "exprds", "region"), all=T),
                 any_highfc_gene_dt, by=c("hicds", "exprds", "region"), all=T), any_highqtfc_gene_dt, by=c("hicds", "exprds", "region"), all=T)

all_dt$any_signif_and_signif <- all_dt$any_signif_gene & all_dt$signif_tad
all_dt$any_high_and_signif <- all_dt$any_high_geneFC & all_dt$signif_tad
all_dt$any_highQt_and_signif <- all_dt$any_highQt_geneFC & all_dt$signif_tad

only_signif_dt <- all_dt[all_dt$signif_tad,]

only_signif_dt$anySignif_anyHigh <- only_signif_dt$any_high_geneFC & only_signif_dt$any_signif_gene
only_signif_dt$anySignif_noHigh <- !only_signif_dt$any_high_geneFC & only_signif_dt$any_signif_gene
only_signif_dt$noSignif_anyHigh <- only_signif_dt$any_high_geneFC & !only_signif_dt$any_signif_gene
only_signif_dt$noSignif_noHigh <- !only_signif_dt$any_high_geneFC & !only_signif_dt$any_signif_gene

only_signif_dt$anySignif_anyHighQt <- only_signif_dt$any_highQt_geneFC & only_signif_dt$any_signif_gene
only_signif_dt$anySignif_noHighQt <- !only_signif_dt$any_highQt_geneFC & only_signif_dt$any_signif_gene
only_signif_dt$noSignif_anyHighQt <- only_signif_dt$any_highQt_geneFC & !only_signif_dt$any_signif_gene
only_signif_dt$noSignif_noHighQt <- !only_signif_dt$any_highQt_geneFC & !only_signif_dt$any_signif_gene


only_signif_dt$region <- NULL
agg_onlysignif_dt <- aggregate(.~hicds+exprds, FUN=sum, data=only_signif_dt)

agg_onlysignif_dt$none_signif_gene <- agg_onlysignif_dt$signif_tad - agg_onlysignif_dt$any_signif_gene
agg_onlysignif_dt$none_high_geneFC <- agg_onlysignif_dt$signif_tad - agg_onlysignif_dt$any_high_geneFC
agg_onlysignif_dt$none_highQt_geneFC <- agg_onlysignif_dt$signif_tad - agg_onlysignif_dt$any_highQt_geneFC

stopifnot(agg_onlysignif_dt$none_highfc_gene >= 0)
stopifnot(agg_onlysignif_dt$none_signif_gene >= 0)

stopifnot(agg_onlysignif_dt$signif_tad == agg_onlysignif_dt$anySignif_noHigh + agg_onlysignif_dt$noSignif_noHigh +
                                                agg_onlysignif_dt$anySignif_anyHigh + agg_onlysignif_dt$noSignif_anyHigh )
stopifnot(agg_onlysignif_dt$any_signif_gene == agg_onlysignif_dt$anySignif_noHigh + agg_onlysignif_dt$anySignif_anyHigh)
stopifnot(agg_onlysignif_dt$any_high_gene == agg_onlysignif_dt$anySignif_anyHigh + agg_onlysignif_dt$noSignif_anyHigh)

stopifnot(agg_onlysignif_dt$signif_tad == agg_onlysignif_dt$anySignif_noHighQt + agg_onlysignif_dt$noSignif_noHighQt +
            agg_onlysignif_dt$anySignif_anyHighQt + agg_onlysignif_dt$noSignif_anyHighQt )
stopifnot(agg_onlysignif_dt$any_signif_gene == agg_onlysignif_dt$anySignif_noHighQt + agg_onlysignif_dt$anySignif_anyHighQt)
stopifnot(agg_onlysignif_dt$any_HighQt_gene == agg_onlysignif_dt$anySignif_anyHighQt + agg_onlysignif_dt$noSignif_anyHighQt)



tmp_dt <- agg_onlysignif_dt[order(agg_onlysignif_dt$signif_tad, decreasing = TRUE),]
tmp_dt$ds <- file.path(tmp_dt$hicds, "\n",tmp_dt$exprds)

plot_dt <- melt(agg_onlysignif_dt, id=c("hicds", "exprds"))
plot_dt$ds <- file.path(plot_dt$hicds, "\n", plot_dt$exprds)

plot_dt$ds <- factor(plot_dt$ds, levels= as.character(tmp_dt$ds))
stopifnot(!is.na(plot_dt$ds))

all_plot_vars <- c("high_geneFC", "highQt_geneFC", "signif_gene")

all_subtits <- setNames(c(paste0("abs. FC >= ", abs_fc_thresh),
                          paste0("abs. FC >= ", abs_fc_qtThresh, "-qt." ),
                          paste0("p-val. <= ", signif_gene_thresh)),
                        all_plot_vars)



foo <- foreach(plot_var = all_plot_vars) %dopar% {


  plotTit <- paste0(plot_var)
  subTit <- paste0("TAD: adj. p-val. signif. <= ", signif_tad_thresh, "; gene ", all_subtits[plot_var])
  myxlab <- paste0("All datasets (n=", length(unique(plot_dt$ds)), ")")
  myylab <- "# signif. TADs"


  p_withLab <- ggbarplot(
    data=plot_dt[as.character(plot_dt$variable) %in% c(paste0("none_", plot_var), paste0("any_", plot_var)), ],
    y = "value", fill = "variable", x ="ds", palette="jco"
  ) +
    ggtitle(paste0(plotTit), subtitle=paste0(subTit)) +
    # scale_color_manual(values="my_cols")+
    # scale_fill_manual(values=my_cols)  +
    # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    labs(x=myxlab, y=myylab, fill="") +
    # guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand=c(0,2))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      # text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      # axis.text.x = element_blank(),
      axis.text.x = element_text(size=6, angle=90, hjust=1, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    )

  outFile <- file.path(outFolder, paste0("nSignifTADs_withAny_none_", plot_var, "_allDS_barplot_withLabs.", "svg"))
  ggsave(p_withLab, file=outFile, height=myHeight*1.6, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))


  p_noLab <- p_withLab +  theme(axis.text.x = element_blank())

  outFile <- file.path(outFolder, paste0("nSignifTADs_withAny_none_", plot_var, "_allDS_barplot_noLabs.", "svg"))
  ggsave(p_noLab, file=outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))

}


all_plot_vars2 <- c("High", "HighQt")

plot_var <- "High"

all_subtits <- setNames(c(paste0("abs. FC >= ", abs_fc_thresh),
                          paste0("abs. FC >= ", abs_fc_qtThresh, "-qt." ),
                          paste0("p-val. <= ", signif_gene_thresh)),
                        all_plot_vars2)


save(plot_dt, file="plot_dt.Rdata", version=2)

plot_var = all_plot_vars2[1]
foo <- foreach(plot_var = all_plot_vars2) %dopar% {
  
  
  all_cols <- c(
    paste0("anySignif_any", plot_var),paste0("anySignif_no", plot_var), paste0("noSignif_any", plot_var), paste0("noSignif_no", plot_var))
  
  stopifnot(all_cols %in% as.character(plot_dt$variable))
  
  plotTit <- paste0(plot_var)
  subTit <- paste0("TAD: adj. p-val. signif. <= ", signif_tad_thresh, "; gene: p-val. signif. <= ", signif_gene_thresh , " and FC ", all_subtits[plot_var])
  myxlab <- paste0("All datasets (n=", length(unique(plot_dt$ds)), ")")
  myylab <- "# signif. TADs"
  
  
  p_withLab <- ggbarplot(
    data=plot_dt[as.character(plot_dt$variable)%in%all_cols,],
    y = "value", fill = "variable", x ="ds", palette="jco"
  ) +
    ggtitle(paste0(plotTit), subtitle=paste0(subTit)) +
    # scale_color_manual(values="my_cols")+
    # scale_fill_manual(values=my_cols)  +
    # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    labs(x=myxlab, y=myylab, fill="") +
    # guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand=c(0,2))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      # text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      # axis.text.x = element_blank(),
      axis.text.x = element_text(size=6, angle=90, hjust=1, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    ) 
  
  outFile <- file.path(outFolder, paste0("nSignifTADs_withAny_none_cmbSignifFC_", plot_var, "_allDS_barplot_withLabs.", "svg"))
  ggsave(p_withLab, file=outFile, height=myHeight*1.6, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  p_noLab <- p_withLab +  theme(axis.text.x = element_blank())
  
  outFile <- file.path(outFolder, paste0("nSignifTADs_withAny_none_cmbSignifFC_", plot_var, "_allDS_barplot_noLabs.", "svg"))
  ggsave(p_noLab, file=outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
}




# plot_var1 <- "signif_gene"
# plot_var2 <- "highQt_geneFC"
# 
# plot_dt$mygroup <- ifelse(grepl("signif_gene", plot_dt$variable) , "signif", "fc")
# 
# 
# ggplot(  plot_dt[as.character(plot_dt$variable) %in% c(paste0("none_", plot_var1), paste0("any_", plot_var1)) |
#                    as.character(plot_dt$variable) %in% c(paste0("none_", plot_var2), paste0("any_", plot_var2)),],
#          aes(y = value, x = mygroup, group=factor(ds), color=variable))+
#   geom_bar(stat="identity")
# 
# # ggplot(  plot_dt[as.character(plot_dt$variable) %in% c(paste0("none_", plot_var1), paste0("any_", plot_var1)) | 
# #                    as.character(plot_dt$variable) %in% c(paste0("none_", plot_var2), paste0("any_", plot_var2)),],
# #          aes(y = value, x = mygroup, fill=variable))+
# #   facet_wrap(~ ds,nrow=1) + #coord_flip()+
# #   geom_bar(stat="identity")
# # + 
# # theme(panel.spacing.y = unit(0,'npc'))
# 
# # set.seed(1234)
# # data <- data.frame(
#   animal = sample(c('bear','tiger','lion'), 50, replace=T),
#   color = sample(c('black','brown','orange'), 50, replace=T),
#   period = sample(c('first','second','third'), 50, replace=T),
#   value = sample(1:100, 50, replace=T))
# ggplot(data, aes(x=period, y=value, fill=color, group=animal, color=animal)) +
#   geom_bar(stat="identity", position="dodge")



# ggbarplot(
#   data=plot_dt[as.character(plot_dt$variable) %in% c(paste0("none_", plot_var1), paste0("any_", plot_var1)) | 
#                  as.character(plot_dt$variable) %in% c(paste0("none_", plot_var2), paste0("any_", plot_var2)),],
#   y = "value", fill = "variable", x ="ds", palette="jco"
# ) +
#   ggtitle(paste0(plotTit), subtitle=paste0(subTit)) +
#   # scale_color_manual(values="my_cols")+
#   # scale_fill_manual(values=my_cols)  +
#   # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
#   labs(x=myxlab, y=myylab, fill="") +
#   # guides(color=FALSE)+
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand=c(0,2))+
#   # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
#   theme(
#     # text = element_text(family=fontFamily),
#     panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
#     panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
#     panel.background = element_rect(fill = "transparent"),
#     panel.grid.major.x =  element_blank(),
#     panel.grid.minor.x =  element_blank(),
#     axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
#     axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
#     axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
#     # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
#     # axis.text.x = element_blank(),
#     axis.text.x = element_text(size=6, angle=90, hjust=1, vjust=0.5),
#     plot.title = element_text(hjust=0.5, size = 16, face="bold"),
#     plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
#     legend.title = element_text(face="bold")
#   ) 


