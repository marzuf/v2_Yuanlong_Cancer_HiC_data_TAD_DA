
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
# Rscript cmp_meanCorr_permut.R

plotType <- "svg"
myHeight <- myWidth <- 7
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "CMP_MEANCORR_PERMUT"
dir.create(outFolder, recursive = TRUE)
  
all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
hicds = all_hicds[1]
all_hicds <- all_hicds[grep("NCI-H460", all_hicds)]

myHicds <- "ENCSR489OCU_NCI-H460_"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_ds_meanCorr <- foreach(hicds = all_hicds) %dopar% {
 get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata")))
}
names(all_ds_meanCorr) <- gsub(myHicds, "", all_hicds)

outFile <- file.path(outFolder, paste0("allDS_meanCorr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  all_ds_meanCorr, 
  plotTit = paste0(myHicds, " - n=", length(all_hicds), " - TAD meanCorr"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("allDS_meanCorr_density_noG2t.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  all_ds_meanCorr[!grepl("PERMUTG2T", names(all_ds_meanCorr))], 
  plotTit = paste0(myHicds, " - n=", length(all_hicds), " - TAD meanCorr"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

tadSignifThresh <- 0.05

plot_dt  <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  mean_corr <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata")))
  empPval <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
  empPval <- p.adjust(empPval, method="BH")
  data.frame(
    hicds= hicds,
    meanCorr = as.numeric(mean_corr),
    adjCombPval = empPval[names(mean_corr)],
    region = names(mean_corr),
    stringsAsFactors = FALSE
  )
}
stopifnot(!is.na(plot_dt))
plot_dt$signif <- ifelse(plot_dt$adjCombPval <= tadSignifThresh, "signif.", "not signif.")
plot_dt$hicds_lab <- gsub(myHicds, "", plot_dt$hicds)

box_meanCorr <- ggboxplot(data= plot_dt, x="signif", y= "meanCorr", xlab="") +
  ggtitle("IntraTAD corr.", subtitle = paste0(myHicds, "; TAD signif. thresh <= ", tadSignifThresh))+
  facet_grid(~hicds_lab, switch="x") + 
  scale_y_continuous(name=paste0("TAD meanCorr"),
                     breaks = scales::pretty_breaks(n = 10))+
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"))

outFile <- file.path(outFolder, paste0("allDS_meanCorr_signifNotSginfi_boxplot.", plotType))
ggsave(box_meanCorr, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

#################################################################################
################################################################################# emp pval
#################################################################################

all_ds_empPval <- foreach(hicds = all_hicds) %dopar% {
 empPval <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
 empPval <- p.adjust(empPval, method="BH")
 -log10(empPval)
}
names(all_ds_empPval) <- gsub(myHicds, "", all_hicds)

outFile <- file.path(outFolder, paste0("allDS_empPval_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  all_ds_empPval, 
  plotTit = paste0(myHicds, "- n=", length(all_hicds), " - adj. emp. p-val [-log10]"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("allDS_empPval_density_noG2t.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  all_ds_empPval[!grepl("PERMUTG2T", names(all_ds_meanCorr))], 
  plotTit = paste0(myHicds, " - n=", length(all_hicds), " - adj. emp. p-val [-log10]"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))










