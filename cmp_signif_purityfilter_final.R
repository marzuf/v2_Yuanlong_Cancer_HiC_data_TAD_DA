outFolder <- "CMP_SIGNIF_PURITYFILTER_FINAL"
dir.create(outFolder, recursive = T)
# Rscript cmp_signif_purityfilter_final.R

init_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
pf_dt <- get(load("CREATE_FINAL_TABLE_PURITYFILTER/all_result_dt.Rdata"))

signifThresh <- 0.01

init_dt$signif_tad <- init_dt$adjPvalComb <= signifThresh
pf_dt$signif_tad <- pf_dt$adjPvalComb <= signifThresh

agg_init_dt <- aggregate(signif_tad ~ hicds+exprds, FUN=sum, data=init_dt)
agg_pf_dt <- aggregate(signif_tad ~ hicds+exprds, FUN=sum, data=pf_dt)


plot_dt <- merge(agg_pf_dt, agg_init_dt, by=c("hicds", "exprds"), suffixes = c("_pF", "_init"), all.x=T, all.y=F)
stopifnot(!is.na(plot_dt))

plotTit <- paste0("# signif. TADs (# ds = ", nrow(plot_dt), ")")
subTit <- paste0("TAD adj. p-val <= ", signifThresh)

my_x <- plot_dt$signif_tad_init
my_y <- plot_dt$signif_tad_pF

plotType <- "svg"
myHeight <- myWidth <- 7
plotCex <- 1.2

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path(outFolder,  paste0("allDS_nbrSignifTADs_pF_vs_init_scatterplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
plot(x=my_x,y=my_y, main=plotTit, 
     pch=16, cex=0.7,
     cex.main=plotCex,
     cex.axis=plotCex,
     cex.lab=plotCex,
     xlab=paste0("# signif. TADs init. data"),
     ylab=paste0("# signif. TADs purity-filtered data"))
curve(1*x, col="grey", add=TRUE)
addCorr(x=my_x, y=my_y, legPos="topleft", bty="n")
mtext(side=3, text = subTit, font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))

