require(foreach)
require(doMC)
registerDoMC(40)

require(ggplot2)
require(ggpubr)
require(ggsci)

source("revision_settings.R")

plotType <- "png"
myHeightGG <- 5
myWidthGG <- 6

# Rscript revision_aggChipPeaks_data_fc.R

outFolder <- "REVISION_AGGCHIPPEAKS_DATA_FC"
dir.create(outFolder, recursive = TRUE)
pvalthresh <- 0.01

inFiles <- list(
  c("subset_chr12_54160001_54440000_22Rv1_ENCFF286BKT_fcOverControl.bedGraph",
"subset_chr12_54160001_54440000_RWPE1_ENCFF039XYU_fcOverControl.bedGraph", "chr12_54160001_54440000"),
c(
"subset_chr17_46720001_46880000_22Rv1_ENCFF286BKT_fcOverControl.bedGraph",
"subset_chr17_46720001_46880000_RWPE1_ENCFF039XYU_fcOverControl.bedGraph", "chr17_46720001_46880000"),
c(
"subset_chr7_116080001_116320000_22Rv1_ENCFF286BKT_fcOverControl.bedGraph",
"subset_chr7_116080001_116320000_RWPE1_ENCFF039XYU_fcOverControl.bedGraph", "chr7_116080001_116320000"))

inFolder <- "prostate_chip_seq/SUBSET_FOR_IGV_BIGWIG"

for(i in seq_along(inFiles)){
  
  plot_tit <- inFiles[[i]][3]
  
  dt_22rv1 <- read.delim(file.path(inFolder, inFiles[[i]][1]), header=F,col.names=c("chromo", "start", "end", "FC"))
  dt_RWPE1 <- read.delim(file.path(inFolder, inFiles[[i]][2]), header=F,col.names=c("chromo", "start", "end", "FC"))
  
  dt_22rv1$dataset <- "22Rv1"
  dt_RWPE1$dataset <- "RWPE1"
  
  plot_dt <- rbind(dt_22rv1,dt_RWPE1)
  plot_dt$FC_log10 <- log10(plot_dt$FC)
  
  fillvar <- "dataset"
  legTitle <- paste0("")
  # my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], sort(unique(sub_plot_dt[,paste0("tad_signif")])))
  my_cols <- setNames(pal_jama()(5)[c(3, 2,4, 1)], sort(unique(plot_dt[,paste0(fillvar)])))
  # plotTit <- paste0(icol)
  plotTit <- paste0(plot_tit)
  
  mySub1 <- paste0("")
  mySub2 <-  paste0("")
  
  
  colplot <- "FC"
  legTitle <- ""
  for(colplot in c("FC", "FC_log10")) {
    myylab <- paste0(colplot)
    mySub <- paste0("FC over control (", colplot, ")")
    
    p3 <- ggboxplot(plot_dt,
                    y = paste0(colplot),
                    x = paste0(fillvar),
                    # combine = TRUE,    
                    add="jitter",
                    outlier.shape=NA,
                    xlab = "",
                    ylab = paste0(myylab),
                    # add = "median",                  # Add median line.
                    rug = FALSE,                      # Add marginal rug
                    color = paste0(fillvar),
                    # fill = paste0(fillvar),
                    # color = paste0("tad_signif"),
                    # fill = paste0("tad_signif"),
                    palette = "jco"
    ) +
      # geom_hline(yintercept=1)+
      ggtitle(plotTit, subtitle = mySub)+
      # scale_color_manual(values=my_cols)+
      # scale_fill_manual(values=my_cols)  +
      labs(color=paste0(legTitle),fill=paste0(legTitle)) +
      guides(color=FALSE, fill=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
      # scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
      mytheme
    
    outFile <- file.path(outFolder, paste0(plot_tit, "_", colplot,"_", paste0(unique(plot_dt[,fillvar]), collapse="_"), "_boxplot.", plotType))
    cat(paste0("... written: ", outFile, "\n"))
    
    ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
  }
  
  
  
}
