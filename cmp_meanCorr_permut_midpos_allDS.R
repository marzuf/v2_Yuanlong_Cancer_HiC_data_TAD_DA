
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
# Rscript cmp_meanCorr_permut_midpos_allDS.R

plotType <- "svg"
myHeight <- myWidth <- 7
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "CMP_MEANCORR_PERMUT_MIDPOS_ALLDS"
dir.create(outFolder, recursive = TRUE)
  
all_hicds_init <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds1 <- all_hicds_init[ ! (grepl("RANDOM", all_hicds_init) | grepl("PERMUT", all_hicds_init)) ]
all_hicds2 <- all_hicds_init[grepl("RANDOMMIDPOS", all_hicds_init)]
all_hicds <- c(all_hicds1, all_hicds2)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))


# 
rd_patterns <- c("RANDOMMIDPOS","RANDOMMIDPOSDISC", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T" )



source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_ds_meanCorr_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
     data.frame(
   hicds = hicds,
   exprds = exprds,
   meanCorr = get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata"))),
   stringsAsFactors = FALSE)
  }
  exprds_dt
}

all_ds_meanCorr_dt$hicds_lab <-  gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                      gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                               gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                    gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                               gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_ds_meanCorr_dt$hicds)))))
all_ds_meanCorr_dt$hicds_lab[! all_ds_meanCorr_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"



outFile <- file.path(outFolder, paste0("allDS_meanCorr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_ds_meanCorr_dt$meanCorr, all_ds_meanCorr_dt$hicds_lab),
  plotTit = paste0("all DS - n=", length(unique(file.path(all_ds_meanCorr_dt$hicds, all_ds_meanCorr_dt$exprds))), " - TAD meanCorr"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("allDS_meanCorr_density_noG2t.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_ds_meanCorr_dt$meanCorr[!grepl("PERMUTG2T", all_ds_meanCorr_dt$hicds_lab)], all_ds_meanCorr_dt$hicds_lab[!grepl("PERMUTG2T", all_ds_meanCorr_dt$hicds_lab)]),
  plotTit = paste0("all DS - n=", length(unique(file.path(all_ds_meanCorr_dt$hicds, all_ds_meanCorr_dt$exprds))), " - TAD meanCorr"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

tadSignifThresh <- 0.01

plot_dt  <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    mean_corr <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata")))
    empPval <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
    empPval <- p.adjust(empPval, method="BH")
    data.frame(
      hicds= hicds,
      exprds=exprds,
      meanCorr = as.numeric(mean_corr),
      adjCombPval = empPval,
      region = names(mean_corr),
      stringsAsFactors = FALSE
    )
  }
  exprds_dt
}

plot_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                                     gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                          gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                               gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", plot_dt$hicds))))
plot_dt$hicds_lab[! plot_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"


plot_dt$signif <- ifelse(plot_dt$adjCombPval <= tadSignifThresh, "signif.", "not signif.")

box_meanCorr <- ggboxplot(data= plot_dt, x="signif", y= "meanCorr", xlab="") +
  ggtitle("IntraTAD corr.", subtitle = paste0("all ds - n=", length(unique(file.path(plot_dt$hicds, plot_dt$exprds)))))+
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

all_ds_empPval_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
     empPval <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")))
     empPval <- p.adjust(empPval, method="BH")
     data.frame(
       hicds= hicds,
       adjCombPval_log10=-log10(empPval),
       stringsAsFactors = FALSE
     )
  }
  exprds_dt
}
all_ds_empPval_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                          gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                               gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                    gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_ds_empPval_dt$hicds))))
all_ds_empPval_dt$hicds_lab[! all_ds_empPval_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"



outFile <- file.path(outFolder, paste0("allDS_empPval_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_ds_empPval_dt$adjCombPval_log10, all_ds_empPval_dt$hicds_lab),
  plotTit = paste0("all DS - n=", length(unique(all_ds_empPval_dt$hicds)), " - adj. emp. p-val [-log10]"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("allDS_empPval_density_noG2t.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_ds_empPval_dt$adjCombPval_log10[!grepl("PERMUTG2T", all_ds_empPval_dt$hicds_lab)], all_ds_empPval_dt$hicds_lab),
  plotTit = paste0( all DS, " - n=", length(unique(all_ds_empPval_dt$hicds)), " - adj. emp. p-val [-log10]"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))










