
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)

# Rscript cmp_FCC_permut_allDS.R

plotType <- "svg"
myHeight <- myWidth <- 7
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "CMP_FCC_PERMUT_ALLDS"
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")

all_hicds <- all_hicds[!( grepl("RANDOMNBRGENES", all_hicds) | grepl("PERMUTG2T", all_hicds) | grepl("RANDOMSHIFT", all_hicds) | grepl("RANDOMMIDPOSSTRICT", all_hicds) )]

all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))



# 
rd_patterns <- c("RANDOMMIDPOS","RANDOMMIDPOSDISC" , "RANDOMMIDPOSSTRICT", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T")
plot_patt <-  c("RANDOMMIDPOS","RANDOMMIDPOSDISC" )


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_ds_fcc_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    fcc_file <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown", "all_obs_prodSignedRatio.Rdata")
    if(!file.exists(fcc_file)) return(NULL)
    data.frame(
      hicds = hicds,
      exprds = exprds,
      FCC = get(load(fcc_file)),
      stringsAsFactors = FALSE)
  }
  exprds_dt
}

all_ds_fcc_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                                     gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                          gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                               gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                                    gsub(".+RANDOMMIDPOSSTRICT_40kb", "RANDOMMIDPOSSTRICT",
                                                         gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_ds_fcc_dt$hicds))))))
all_ds_fcc_dt$hicds_lab[! all_ds_fcc_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"

save(all_ds_fcc_dt, file="all_ds_fcc_dt.Rdata", version=2)

all_ds_fcc_dt$hicds_lab <- factor(all_ds_fcc_dt$hicds_lab, levels = c("OBSERVED", rd_patterns))

box_fcc <- ggboxplot(data= all_ds_fcc_dt, x="hicds_lab", y= "FCC", xlab="") +
  ggtitle("TAD FCC score", subtitle = paste0("all ds - n=", length(unique(file.path(all_ds_fcc_dt$hicds, all_ds_fcc_dt$exprds)))))+
  # facet_grid(~hicds_lab, switch="x") + 
  scale_y_continuous(name=paste0("FCC"),
                     breaks = scales::pretty_breaks(n = 10))+
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"))

outFile <- file.path(outFolder, paste0("allDS_FCC_boxplot.", plotType))
ggsave(box_fcc, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

all_ds_fcc_dt$hicds_lab <- as.character(all_ds_fcc_dt$hicds_lab)

outFile <- file.path(outFolder, paste0("allDS_FCC_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_ds_fcc_dt$FCC, all_ds_fcc_dt$hicds_lab), legPos = "topleft",
  plotTit = paste0("all DS - FCC scores"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



sub_dt <- all_ds_fcc_dt[as.character(all_ds_fcc_dt$hicds_lab) %in% c("OBSERVED",plot_patt ),]
sub_dt$hicds_lab <- factor(sub_dt$hicds_lab, levels = c("OBSERVED",plot_patt ))


save(sub_dt, version=2, file="sub_dt.Rdata")

box_fcc <- ggboxplot(data= sub_dt, x="hicds_lab", y= "FCC", xlab="") +
  ggtitle("TAD FCC score", subtitle = paste0("all ds - n=", length(unique(file.path(sub_dt$hicds, sub_dt$exprds)))))+
  # facet_grid(~hicds_lab, switch="x") + 
  scale_y_continuous(name=paste0("FCC"),
                     breaks = scales::pretty_breaks(n = 10))+
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"))

outFile <- file.path(outFolder, paste0("allDS_FCC_boxplot_sub.", plotType))
ggsave(box_fcc, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

sub_dt$hicds_lab <- as.character(sub_dt$hicds_lab)


outFile <- file.path(outFolder, paste0("allDS_FCC_density_sub.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(sub_dt$FCC, sub_dt$hicds_lab),
  plotTit = paste0("all DS - FCC scores"),  legPos = "topleft")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))









