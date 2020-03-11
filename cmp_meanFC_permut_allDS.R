
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
# Rscript cmp_meanFC_permut_allDS.R

plotType <- "svg"
myHeight <- myWidth <- 7
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "CMP_MEANFC_PERMUT_ALLDS"
dir.create(outFolder, recursive = TRUE)
  
all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))


# 
rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T", "RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT")



source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_ds_meanFC_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    stopifnot(file.exists(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "3_runMeanTADLogFC", "all_meanLogFC_TAD.Rdata")))
     data.frame(
   hicds = hicds,
   exprds = exprds,
   
   meanFC = get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "3_runMeanTADLogFC", "all_meanLogFC_TAD.Rdata"))),
   stringsAsFactors = FALSE)
  }
  exprds_dt
}

all_ds_meanFC_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                                   gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                   gsub(".+RANDOMMIDPOSSTRICT_40kb", "RANDOMMIDPOSSTRICT",
                               gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                    gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                               gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_ds_meanFC_dt$hicds))))))
all_ds_meanFC_dt$hicds_lab[! all_ds_meanFC_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"



outFile <- file.path(outFolder, paste0("allDS_meanFC_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_ds_meanFC_dt$meanFC, all_ds_meanFC_dt$hicds_lab),
  plotTit = paste0("all DS - n=", length(unique(file.path(all_ds_meanFC_dt$hicds, all_ds_meanFC_dt$exprds))), " - TAD meanLogFC"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("allDS_meanFC_density_noG2t.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_ds_meanFC_dt$meanFC[!grepl("PERMUTG2T", all_ds_meanFC_dt$hicds_lab)], all_ds_meanFC_dt$hicds_lab[!grepl("PERMUTG2T", all_ds_meanFC_dt$hicds_lab)]),
  plotTit = paste0("all DS - n=", length(unique(file.path(all_ds_meanFC_dt$hicds, all_ds_meanFC_dt$exprds))), " - TAD meanLogFC"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

tadSignifThresh <- 0.01

plot_dt  <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    mean_fc <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "3_runMeanTADLogFC", "all_meanLogFC_TAD.Rdata")))
    empPval <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
    empPval <- p.adjust(empPval, method="BH")
    data.frame(
      hicds= hicds,
      exprds=exprds,
      meanFC = as.numeric(mean_fc),
      adjCombPval = empPval,
      region = names(mean_fc),
      stringsAsFactors = FALSE
    )
  }
  exprds_dt
}

plot_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                          gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                     gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                          gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                               gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", plot_dt$hicds)))))
plot_dt$hicds_lab[! plot_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"


plot_dt$signif <- ifelse(plot_dt$adjCombPval <= tadSignifThresh, "signif.", "not signif.")

box_meanFC <- ggboxplot(data= plot_dt, x="signif", y= "meanFC", xlab="") +
  ggtitle("TAD meanLogFC", subtitle = paste0("all ds - n=", length(unique(plot_dt$hicds))))+
  facet_grid(~hicds_lab, switch="x") + 
  scale_y_continuous(name=paste0("TAD meanFC"),
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

outFile <- file.path(outFolder, paste0("allDS_meanFC_signifNotSginfi_boxplot.", plotType))
ggsave(box_meanFC, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



