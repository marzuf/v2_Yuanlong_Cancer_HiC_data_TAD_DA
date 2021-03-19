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
"subset_chr12_54160001_54440000_RWPE1_ENCFF039XYU_fcOverControl.bedGraph"),
c(
"subset_chr17_46720001_46880000_22Rv1_ENCFF286BKT_fcOverControl.bedGraph",
"subset_chr17_46720001_46880000_RWPE1_ENCFF039XYU_fcOverControl.bedGraph"),
c(
"subset_chr7_116080001_116320000_22Rv1_ENCFF286BKT_fcOverControl.bedGraph",
"subset_chr7_116080001_116320000_RWPE1_ENCFF039XYU_fcOverControl.bedGraph"))

inFolder <- "prostate_chip_seq/SUBSET_FOR_IGV_BIGWIG"

for(i in seq_along(inFiles)){
  
  dt_22rv1 <- read.delim(file.path(inFolder, inFiles[[i]][1]), header=F,col.names=c("chromo", "start", "end", "FC"))
  dt_RWPE1 <- read.delim(file.path(inFolder, inFiles[[i]][2]), header=F,col.names=c("chromo", "start", "end", "FC"))
  
  dt_22rv1$dataset <- "22Rv1"
  dt_RWPE1$dataset <- "RWPE1"
  
  plot_dt <- rbind(dt_22rv1,dt_RWPE1)
  
  fillvar <- "dataset"
  legTitle <- paste0("")
  # my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], sort(unique(sub_plot_dt[,paste0("tad_signif")])))
  my_cols <- setNames(pal_jama()(5)[c(3, 2,4, 1)], sort(unique(plot_dt[,paste0(fillvar)])))
  # plotTit <- paste0(icol)
  plotTit <- paste0("FC over control")
  
  mySub1 <- paste0("")
  mySub2 <-  paste0("")
  
  mySub <- paste0(mySub1, "; ", mySub2)
  
  colplot <- "FC"
  myylab <- paste0(colplot)
  p3 <- ggboxplot(plot_dt,
                  y = paste0(colplot),
                  x = paste0(fillvar),
                  # combine = TRUE,    
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
    geom_hline(yintercept=1)+
    ggtitle(plotTit, subtitle = mySub)+
    scale_color_manual(values=my_cols)+
    scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle)) +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
    mytheme
  
  outFile <- file.path(outFolder, paste0(plotvar, "_ratio_", ref_lab, "_vs_",  match_lab, "_boxplot.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
}


runFolder <- "."
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
resultData$region_id <- file.path(resultData$dataset, resultData$region)
resultData$signif_lab1 <- ifelse(resultData$adjPvalComb <= pvalthresh, "signif", "notsignif")
resultData$direction_lab <- ifelse(resultData$meanLogFC <0, "down", 
                                   ifelse(resultData$meanLogFC >0, "up",NA))
stopifnot(!is.na(resultData$direction_lab))
resultData$signif_lab2 <- paste0(resultData$signif_lab1, ".", resultData$direction_lab)
resultData$signif_lab3 <- resultData$signif_lab2
resultData$signif_lab3[resultData$signif_lab1 == "notsignif"] <- "notsignif"
resultData$chromo <- gsub("(chr.+)_TAD.+", "\\1", resultData$region)
stopifnot(resultData$chromo %in% paste0("chr", 1:22))

exprds <- "TCGAprad_norm_prad"

ref_hicds <- "GSE118514_RWPE1_40kb"
match_hicds <- "GSE118514_22Rv1_40kb"

# match_hicds<- "GSE118514_RWPE1_40kb"
# ref_hicds <- "GSE118514_22Rv1_40kb"

ref_lab <- gsub(".+_(.+)_40kb", "\\1", ref_hicds)
match_lab <- gsub(".+_(.+)_40kb", "\\1", match_hicds)

hicds_chips <- c( "GSE118514_RWPE1_40kb" = "RWPE1_ENCFF149SDU",
                "GSE118514_22Rv1_40kb" = "22Rv1_ENCFF655WXZ")

peakFolder <- "prostate_chip_seq"
ref_peaks_file <- file.path(peakFolder, paste0(hicds_chips[paste0(ref_hicds)], "_repPeaks.bed"))
ref_peaks_dt <- read.delim(ref_peaks_file, header=FALSE, col.names = 
                             c("chrom" ,"chromStart", "chromEnd" ,"name","score","strand" ,"signalValue","pValue","qValue","peak"),
                           stringsAsFactors = FALSE)
stopifnot(is.numeric(ref_peaks_dt$chromStart))
stopifnot(is.numeric(ref_peaks_dt$chromEnd))
stopifnot(resultData$chromo %in% ref_peaks_dt$chrom )
# take only signif. peaks
sum(ref_peaks_dt$qValue < -log10(0.05)) # 18
ref_peaks_dt <- ref_peaks_dt[ref_peaks_dt$qValue >= -log10(0.05),]
# [1] 95495

match_peaks_file <- file.path(peakFolder, paste0(hicds_chips[paste0(match_hicds)], "_repPeaks.bed"))
match_peaks_dt <- read.delim(match_peaks_file, header=FALSE, col.names = 
                             c("chrom" ,"chromStart", "chromEnd" ,"name","score","strand" ,"signalValue","pValue","qValue","peak"),
                             stringsAsFactors = FALSE)
stopifnot(is.numeric(match_peaks_dt$chromStart))
stopifnot(is.numeric(match_peaks_dt$chromEnd))
stopifnot(resultData$chromo %in% match_peaks_dt$chrom )
# take only signif. peaks
sum(match_peaks_dt$qValue < -log10(0.05)) # 21
match_peaks_dt <- match_peaks_dt[match_peaks_dt$qValue >= -log10(0.05),]
# [1] 45955
# 

ref_final_dt <- resultData[resultData$hicds == ref_hicds & 
                             resultData$exprds == exprds, ]
stopifnot(nrow(ref_final_dt)> 0)

ref_final_dt$match_meanQvalue <-ref_final_dt$match_nbrPeaks <- ref_final_dt$ref_meanQvalue <- ref_final_dt$ref_nbrPeaks <- NA



i_tad=1
tad_withPeaks_dt <- foreach(i_tad = 1:nrow(ref_final_dt), .combine='rbind') %dopar% {
  
  if(i_tad==4) save(ref_final_dt, file="ref_final_dt.Rdata", version=2)
  if(i_tad==4) save(match_peaks_dt, file="match_peaks_dt.Rdata", version=2)
  if(i_tad==4) save(ref_peaks_dt, file="ref_peaks_dt.Rdata", version=2)
  
  curr_chromo="chr12"	;curr_start=54160001;	curr_end=54440000
  
  curr_chromo <- ref_final_dt$chromo[i_tad]
  curr_start <- ref_final_dt$start[i_tad]
  curr_end <- ref_final_dt$end[i_tad]

  sub_match_dt <- match_peaks_dt[match_peaks_dt$chrom == curr_chromo &
                               match_peaks_dt$chromStart >= curr_start &
                               match_peaks_dt$chromEnd <= curr_end ,
                               ]
  # if(nrow(sub_match_dt) > 0 ) {
    ref_final_dt$match_meanQvalue[i_tad] <- mean(sub_match_dt$qValue)
    ref_final_dt$match_nbrPeaks[i_tad] <- nrow(sub_match_dt)
    
  # } else {
  #   ref_final_dt$match_meanQvalue[i_tad] <- NA
  #   ref_final_dt$match_nbrPeaks[i_tad] <- NA
  #   
  # }
  
  
  sub_ref_dt <- ref_peaks_dt[ref_peaks_dt$chrom == curr_chromo &
                               ref_peaks_dt$chromStart >= curr_start &
                               ref_peaks_dt$chromEnd <= curr_end ,
  ]
  
  # if(nrow(sub_match_dt) > 0 ) {
    ref_final_dt$ref_meanQvalue[i_tad] <- mean(sub_ref_dt$qValue)
    ref_final_dt$ref_nbrPeaks[i_tad] <- nrow(sub_ref_dt)
  #   
  #   
  # } else {
  #   ref_final_dt$ref_meanQvalue[i_tad] <- NA
  #   ref_final_dt$ref_nbrPeaks[i_tad] <- NA
  #   
  # }
  # 
  
  
  ref_final_dt[i_tad,]
}

# tad_withPeaks_dt <- ref_final_dt
outFile <- file.path(outFolder, paste0(ref_hicds, "_TADs_match_", match_hicds, "_", "tad_withPeaks_dt.Rdata"))
save(tad_withPeaks_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

plotvar <- "nbrPeaks"
fillvar <- "signif_lab3"

for(plotvar in c("nbrPeaks", "meanQvalue")) {
  
  tad_withPeaks_dt[,paste0("ratio_refMatch_", plotvar)] <- tad_withPeaks_dt[,paste0("ref_", plotvar)]/
    tad_withPeaks_dt[,paste0("match_", plotvar)]
  
  colplot <- paste0("ratio_refMatch_", plotvar)
  
  tad_withPeaks_dt[,colplot] <- tad_withPeaks_dt[,paste0("ref_", plotvar)]/
    tad_withPeaks_dt[,paste0("match_", plotvar)]
  
  plot_dt <- na.omit(tad_withPeaks_dt[,c(colplot, fillvar)])
  
  legTitle <- paste0("")
  # my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], sort(unique(sub_plot_dt[,paste0("tad_signif")])))
  my_cols <- setNames(pal_jama()(5)[c(3, 2,4, 1)], sort(unique(plot_dt[,paste0(fillvar)])))
  # plotTit <- paste0(icol)
  plotTit <- paste0(plotvar, " ratio ", ref_lab, " (ref.)/", match_lab)
  
  mySub1 <- paste0(nrow(plot_dt), "/", nrow(tad_withPeaks_dt))
  mySub2 <-  paste0(names(table(plot_dt[,fillvar])),"=", as.numeric(table(plot_dt[,fillvar])), collapse="; ") 
  
  mySub <- paste0(mySub1, "; ", mySub2)
  
  myylab <- paste0(colplot)
  
  # p3 <- ggdensity(plot_dt,
  #                 x = paste0(colplot),
  #                 y = "..density..",
  #                 # combine = TRUE,                  # Combine the 3 plots
  #                 xlab = paste0(myxlab),
  #                 # add = "median",                  # Add median line.
  #                 rug = FALSE,                      # Add marginal rug
  #                 color = paste0(fillvar),
  #                 fill = paste0(fillvar),
  #                 # color = paste0("tad_signif"),
  #                 # fill = paste0("tad_signif"),
  #                 palette = "jco"
  # ) +
  #   ggtitle(plotTit, subtitle = mySub)+
  #   scale_color_manual(values=my_cols)+
  #   scale_fill_manual(values=my_cols)  +
  #   labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  #   guides(color=FALSE)+
  #   scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  #   scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  #   mytheme
  # 
  # outFile <- file.path(outFolder, paste0(icol, "_", ref_hicds_col, "_vs_", match_hicds_col, "_", plotvar, "_density.", plotType))
  # ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  # cat(paste0("... written: ", outFile, "\n"))
  
  
  p3 <- ggboxplot(plot_dt,
                  y = paste0(colplot),
                  x = paste0(fillvar),
                  # combine = TRUE,    
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
    geom_hline(yintercept=1)+
    ggtitle(plotTit, subtitle = mySub)+
    scale_color_manual(values=my_cols)+
    scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle)) +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
    mytheme
  
  outFile <- file.path(outFolder, paste0(plotvar, "_ratio_", ref_lab, "_vs_",  match_lab, "_boxplot.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}
for(plotvar in c("nbrPeaks", "meanQvalue")) {
  

  colplot <- paste0("ref_", plotvar)
  stopifnot(colplot %in% colnames(tad_withPeaks_dt))
  

  plot_dt <- na.omit(tad_withPeaks_dt[,c(colplot, fillvar)])
  
  legTitle <- paste0("")
  # my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], sort(unique(sub_plot_dt[,paste0("tad_signif")])))
  my_cols <- setNames(pal_jama()(5)[c(3, 2,4, 1)], sort(unique(plot_dt[,paste0(fillvar)])))
  # plotTit <- paste0(icol)
  plotTit <- paste0(plotvar, " ", ref_lab, " (ref.)")
  
  mySub1 <- paste0(nrow(plot_dt), "/", nrow(tad_withPeaks_dt))
  mySub2 <-  paste0(names(table(plot_dt[,fillvar])),"=", as.numeric(table(plot_dt[,fillvar])), collapse="; ") 
  
  mySub <- paste0(mySub1, "; ", mySub2)
  
  myylab <- paste0(colplot)
  
  # p3 <- ggdensity(plot_dt,
  #                 x = paste0(colplot),
  #                 y = "..density..",
  #                 # combine = TRUE,                  # Combine the 3 plots
  #                 xlab = paste0(myxlab),
  #                 # add = "median",                  # Add median line.
  #                 rug = FALSE,                      # Add marginal rug
  #                 color = paste0(fillvar),
  #                 fill = paste0(fillvar),
  #                 # color = paste0("tad_signif"),
  #                 # fill = paste0("tad_signif"),
  #                 palette = "jco"
  # ) +
  #   ggtitle(plotTit, subtitle = mySub)+
  #   scale_color_manual(values=my_cols)+
  #   scale_fill_manual(values=my_cols)  +
  #   labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  #   guides(color=FALSE)+
  #   scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  #   scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  #   mytheme
  # 
  # outFile <- file.path(outFolder, paste0(icol, "_", ref_hicds_col, "_vs_", match_hicds_col, "_", plotvar, "_density.", plotType))
  # ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  # cat(paste0("... written: ", outFile, "\n"))
  
  
  p3 <- ggboxplot(plot_dt,
                  y = paste0(colplot),
                  x = paste0(fillvar),
                  # combine = TRUE,    
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
    geom_hline(yintercept=1)+
    ggtitle(plotTit, subtitle = mySub)+
    scale_color_manual(values=my_cols)+
    scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle)) +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
    mytheme
  
  outFile <- file.path(outFolder, paste0(plotvar, "_", ref_lab, "_boxplot.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}


