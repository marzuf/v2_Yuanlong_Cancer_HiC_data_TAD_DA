plotType <- "png"
myHeightGG <- 6
myWidthGG <- 7

outFolder <- "BIGWIGOVERTAD_AND_SIGNIF"
dir.create(outFolder, recursive = TRUE)

# Rscript bigWigOverTAD_and_signif.R

require(ggsci)
require(ggpubr)
require(ggplot2)

pvalthresh <- 0.05
runFolder <- ".."
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
resultData$region_id <- file.path(resultData$dataset, resultData$region)
resultData$signif_lab <- ifelse(resultData$adjPvalComb <= pvalthresh, "signif", "notsignif")
resultData$direction_lab <- ifelse(resultData$meanLogFC <0, "down", 
                                   ifelse(resultData$meanLogFC >0, "up",NA))
stopifnot(!is.na(resultData$direction_lab))
resultData$signif_lab <- paste0(resultData$signif_lab, ".", resultData$direction_lab)
signif_labs <- setNames(resultData$signif_lab, resultData$region_id)

#   name - name field from bed, which should be unique
#   size - size of bed (sum of exon sizes
#   covered - # bases within exons covered by bigWig
#   sum - sum of values over all bases covered
#   mean0 - average over bases with non-covered bases counting as zeroes
#   mean - average over just covered bases

exprds <- "TCGAprad_norm_prad"


all_files <- c("chip_22Rv1_cover_22Rv1_TADs.bed",
               "chip_RWPE1_cover_RWPE1_TADs.bed")

for(infile in all_files) {
  
  cat(paste0("... start ", infile, "\n"))
  
  
  hicds <- ifelse(grepl("22Rv1", infile), "GSE118514_22Rv1_40kb",
                  ifelse(grepl("RWPE1", infile), "GSE118514_RWPE1_40kb", NA))
  stopifnot(!is.na(hicds))
  
  sub_result_dt <- resultData[resultData$hicds == hicds & resultData$exprds == exprds,]
  
  bw_over_tad_dt <- read.delim(infile,
                               header=F,
                               col.names=c("region", "size", "covered", "sum", "mean0", "mean"),
                               stringsAsFactors = FALSE)
  bw_over_tad_dt$region_id <- file.path(hicds, exprds, bw_over_tad_dt$region)
  
  merge_dt <- merge(sub_result_dt, bw_over_tad_dt, by="region_id", all.x=T, all.y=F)
  stopifnot(!is.na(merge_dt))
  
  plotcols <- c( "sum", "mean0", "mean")
  stopifnot(plotcols %in% colnames(bw_over_tad_dt))
  
  merge_dt$signif_lab2 <- gsub("notsignif.+", "notsignif", merge_dt$signif_lab)
  
  for(toplot in plotcols) {
    
    cat(paste0("... start ", infile, " - ", toplot, "\n"))
    
    
    mysub2 <- paste0(names(table(merge_dt$signif_lab2)), "=", as.numeric(table(merge_dt$signif_lab2)), 
                     collapse=";")
    
    plotTit <- gsub("\\.bed", "", infile)
    legTitle <- ""
    mySub <- paste0(toplot, " - ", mysub2)
    
   p3 <-  ggboxplot(merge_dt,
              y = paste0(toplot),
              outlier.shape=NA,
              add="jitter",
              rug = FALSE,  
              ylab = paste0(toplot),
              xlab ="",
              x="signif_lab2",                   # fill = paste0("tad_signif"),
              color = "signif_lab2",
              # fill = "signif_lab2",
              palette = "jco") +
     ggtitle(plotTit, subtitle = mySub) +
     labs(color=paste0(legTitle),fill=paste0(legTitle)) +
     guides(color=FALSE, fill=FALSE)+
     scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
   
   
   outFile <- file.path(outFolder, paste0(plotTit, "_", toplot, "_boxplot.", plotType))
   ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
   cat(paste0("... written: ", outFile, "\n"))
   
   merge_dt$toplot_log10 <- log10(merge_dt[,toplot])
   
   p3 <-  ggboxplot(merge_dt,
                    y = paste0("toplot_log10"),
                    outlier.shape=NA,
                    add="jitter",
                    rug = FALSE,  
                    ylab = paste0(toplot, " [log10]"),
                    xlab ="",
                    x="signif_lab2",                   # fill = paste0("tad_signif"),
                    color = "signif_lab2",
                    # fill = "signif_lab2",
                    palette = "jco") +
     ggtitle(plotTit, subtitle = mySub) +
     labs(color=paste0(legTitle),fill=paste0(legTitle)) +
     guides(color=FALSE, fill=FALSE)+
     scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
   
   
   outFile <- file.path(outFolder, paste0(plotTit, "_", toplot, "_log10_boxplot.", plotType))
   ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
   cat(paste0("... written: ", outFile, "\n"))
   
   
  }
  
  
}


