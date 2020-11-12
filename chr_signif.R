outFolder <- file.path("CHR_SIGNIF")
dir.create(outFolder)

# Rscript chr_signif.R

require(ggpubr)

result_dt <- get(load("CREATE_FINAL_TABLE//all_result_dt.Rdata"))

result_dt$chromo <- gsub("(chr.+)_.+", "\\1", result_dt$region)

signif_thresh <- 0.01

plotType <- "svg"
myHeight <- 6
myWidth <- 8

agg_dt <- aggregate(adjPvalComb~hicds + exprds + chromo, data=result_dt, FUN =function(x)sum(x<= signif_thresh))

agg_dt$chromo <- factor(agg_dt$chromo, levels=paste0("chr", 1:22))
stopifnot(!is.na(agg_dt$chromo))

size_dt <- read.delim("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes", 
                      header=F, stringsAsFactors = FALSE, sep="\t", col.names = c("chromo", "size"))

chromoSize <- setNames(size_dt$size, size_dt$chromo)

signif_dt <- result_dt[result_dt$adjPvalComb <= signif_thresh,]
signif_dt$midPos <- 0.5*(signif_dt$start + signif_dt$end)
signif_dt$midPos_rel <- signif_dt$midPos/chromoSize[signif_dt$chromo]
stopifnot(!is.na(signif_dt$midPos_rel))

signif_dt$chromo <- factor(signif_dt$chromo, levels=paste0("chr", 1:22))

signif_dt$chromo_num <- as.numeric(signif_dt$chromo)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path(outFolder, paste0("signifTADs_dist_by_chromo.", "svg"))
do.call(plotType, list(outFile, height=myHeight, width=myHeight))
densplot(
  main = "Distribution of signif. TADs",
  x =   signif_dt$chromo_num,
  y = signif_dt$midPos_rel,
  xlab ="chromosome",
  ylab ="relative position"
)
mtext(side=3, text = paste0("all DS - n=", length(unique(file.path(signif_dt$hicds, signif_dt$exprds))), "; signif. p-val. thresh <= ", signif_thresh ))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
# stop("-ok")

pBox <- ggboxplot(agg_dt, y = paste0("adjPvalComb"), x ="chromo", add="jitter",
          xlab ="", ylab="# signif. TADs") + 
  ggtitle(paste0("# signif. TADs by chromo (all DS; n=", length(unique(file.path(agg_dt$hicds, agg_dt$exprds))) , ")"), subtitle=paste0("TAD p-val. signif. thresh. <= ", signif_thresh)) + 
  theme(
    axis.text.x=element_text(angle=90),
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
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("nbrSignifTADs_by_chromo.", "svg"))
ggsave(pBox, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))
# stop("-ok")
