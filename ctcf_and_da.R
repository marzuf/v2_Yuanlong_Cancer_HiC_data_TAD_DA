
# Rscript ctcf_and_da.R

library("readxl")
library(doMC)
library(foreach)
library(stringr)


registerDoMC(40)
# runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878" #PIPELINE/OUTPUT_FOLDER/GM12878_40kb/TCGAluad_norm_luad/11sameNbr_runEmpPvalCombined/"
# hicds <- "GM12878_40kb"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.2



do_densplot_withCorr <- function(xvar, yvar, plot_dt) {
  my_x <- plot_dt[,paste0(xvar)]
  my_y <- plot_dt[,paste0(yvar)]
  densplot(
    x=my_x,
    y=my_y,
    xlab=paste0(xvar),
    ylab=paste0(yvar),
    cex.main=plotCex,
    cex.axis=plotCex,
    cex.lab = plotCex,
    pch=16
  )
  addCorr(x=my_x, y=my_y, bty="n")
}



runFolder <- "." 
hicds <- "ENCSR444WCZ_A549_40kb"
exprds <- "TCGAluad_mutKRAS_mutEGFR"

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
ds_final_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds, ]
stopifnot(nrow(ds_final_dt) > 0)

buildTable <- TRUE

outFolder <- file.path("CTCF_AND_DA", hicds, exprds)
dir.create(outFolder, recursive = TRUE)

tad_dt <- read.delim(file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt"), stringsAsFactors = FALSE, 
                     header=F, col.names = c("chromo", "region", "start", "end"))

### KEEP ONLY TAD REGIONS
tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
stopifnot(!grepl("BOUND", tad_dt$region))

ctcf_dt <- read_excel("13059_2020_2108_MOESM2_ESM.xlsx", sheet="CTCFs")
ctcf_dt <- as.data.frame(ctcf_dt)
ctcf_dt <- ctcf_dt[, 1:7]

# assign ctcf BS to tads
ctcf_dt$chr <- as.character(ctcf_dt$chr)
ctcf_dt <- ctcf_dt[ctcf_dt$chr %in% tad_dt$chromo,]
stopifnot(nrow(ctcf_dt) > 0)
stopifnot(is.numeric(ctcf_dt$start))
stopifnot(is.numeric(ctcf_dt$end))
stopifnot(ctcf_dt$start <= ctcf_dt$end)
ctcf_dt$region <- NA

i=1
if(buildTable){
  ctcf2tad_dt <- foreach(i = 1:nrow(ctcf_dt), .combine='rbind') %dopar% {
    
    chr <- ctcf_dt$chr[i]
    ctcf_start <- ctcf_dt$start[i]
    ctcf_end <- ctcf_dt$end[i]
    
    subtad_dt <- tad_dt[tad_dt$chromo == chr,]
    stopifnot(nrow(subtad_dt) > 0)
    
    # assign if start after tad start and end before tad end
    test1 <- which(ctcf_start >= subtad_dt$start & ctcf_end <= subtad_dt$end)
    test2 <- which(subtad_dt$start <= ctcf_start & subtad_dt$end >= ctcf_end)
    stopifnot(test1==test2)
    stopifnot(length(test1) == 0 | length(test1) == 1)
    if(length(test1) == 1) {
      ctcf_dt$region[i] <- subtad_dt$region[test1]
    } else {
      ctcf_dt$region[i] <- NA
    }
    ctcf_dt[i,]  
  }
  outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
  save(ctcf2tad_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,"\n"))
} else {
  outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
  ctcf2tad_dt <- get(load(outFile))
}

###################### todo: TO CAHNGE KEEP HICDS + EXPRDS + KEEP THE 0s !!!


# load("CTCF_AND_DA/ENCSR444WCZ_A549_40kb/TCGAluad_mutKRAS_mutEGFR/ctcf2tad_dt.Rdata")
cat(paste0("# of CTCF BS:\t",nrow(ctcf2tad_dt), "\n"))
cat(paste0("# of CTCF BS in TADs:\t",sum(!is.na(ctcf2tad_dt$region)), "\n"))
cat(paste0("# of CTCF BS out of TADs:\t",sum(is.na(ctcf2tad_dt$region)), "\n"))

colnames(ds_final_dt)[colnames(ds_final_dt) == "start"] <- "tad_start"
colnames(ds_final_dt)[colnames(ds_final_dt) == "end"] <- "tad_end"

merged_dt <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanLogFC", "meanCorr", "adjPvalComb", "tad_start", "tad_end")],
                   ctcf2tad_dt,
                   by="region", all=FALSE)
tmp <- merged_dt$region
tmp <- gsub("(.+)_.+", "\\1", tmp)
stopifnot(tmp == merged_dt$chr)

merged_dt <- merged_dt[order(merged_dt$chr, merged_dt$start, merged_dt$end ),]
# stopifnot(diff(merged_dt$start) >= 0) # not true because multiple chromo
merged_dt$region <- as.character(merged_dt$region)

#############################################################################
# presence of motifs
#############################################################################

aggByOrientation_dt <- aggregate(chr ~ orientation + region + meanCorr + meanLogFC + adjPvalComb, FUN=length, data=merged_dt)
colnames(aggByOrientation_dt)[colnames(aggByOrientation_dt) == "chr"] <- "CTCF_count"

aggByOrientation_dt$orientation_lab <- ifelse(aggByOrientation_dt$orientation == ">", "forward", 
                                              ifelse(aggByOrientation_dt$orientation == "<", "reverse", NA))
stopifnot(!is.na(aggByOrientation_dt$orientation_lab))

wide_aggByOrientation_dt <- reshape(aggByOrientation_dt[,c("region", "orientation_lab","CTCF_count")], 
                                    idvar="region", direction="wide", timevar = "orientation_lab")



agg_dt <- aggregate(chr ~ region + meanCorr + meanLogFC + adjPvalComb, FUN=length, data=merged_dt)
colnames(agg_dt)[colnames(agg_dt) == "chr"] <- "CTCF_totCount"

agg_merged_dt <- merge(wide_aggByOrientation_dt, agg_dt, by="region", all=TRUE)

agg_merged_dt$CTCF_count.forward[is.na(agg_merged_dt$CTCF_count.forward)] <- 0
agg_merged_dt$CTCF_count.reverse[is.na(agg_merged_dt$CTCF_count.reverse)] <- 0
stopifnot(!is.na(agg_merged_dt))


stopifnot(agg_merged_dt$CTCF_count.reverse + agg_merged_dt$CTCF_count.forward == agg_merged_dt$CTCF_totCount)

agg_merged_dt$adjPvalComb_log10 <- -log10(agg_merged_dt$adjPvalComb)

ntot <- nrow(agg_merged_dt)
nDS <- length(unique(file.path(agg_merged_dt$hicds, agg_merged_dt$exprds)))


all_yvars <-  c("meanCorr", "meanLogFC", "adjPvalComb_log10")
xvar <- "CTCF_totCount"
yvar <- "adjPvalComb_log10"
for(yvar in all_yvars){
  outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_densplot_withCorr(xvar, yvar, agg_merged_dt)
  mtext(side=3, text = paste0("n=", ntot, "; nDS=", nDS))
  title(main = paste0(yvar, " vs. ", xvar))
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
}

#############################################################################
# look by cluster
#############################################################################

if(buildTable) {
  clustByTAD_dt <- foreach(region = unique(merged_dt$region), .combine='rbind') %dopar% {
    
    tad_agg_dt <- merged_dt[merged_dt$region == region,]  
    stopifnot(diff(tad_agg_dt$start) >= 0) 
    
    nConvergent <- str_count(paste0(tad_agg_dt$orientation, collapse=""), "><")
    
    data.frame(
      region  = region,
      orientation = rle(tad_agg_dt$orientation)$values,
      nInClust = rle(tad_agg_dt$orientation)$lengths,
      nConvergent = nConvergent,
      stringsAsFactors=FALSE
    )
  }
  outFile <- file.path(outFolder, "clustByTAD_dt.Rdata")
  save(clustByTAD_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,"\n"))
} else {
  outFile <- file.path(outFolder, "clustByTAD_dt.Rdata")
  clustByTAD_dt <- get(load(outFile))
}
# load("CTCF_AND_DA/ENCSR444WCZ_A549_40kb/TCGAluad_mutKRAS_mutEGFR/clustByTAD_dt.Rdata")

clust_merged_dt <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanLogFC", "meanCorr", "adjPvalComb", "tad_start", "tad_end")],
                   clustByTAD_dt,
                   by="region", all.x=TRUE, all.y=FALSE)
clust_merged_dt$chr <- gsub("(.+)_.+", "\\1", as.character(clust_merged_dt$region))
stopifnot(clust_merged_dt$chr %in% tad_dt$chromo)

tad_conv_dt <- clust_merged_dt
tad_conv_dt$nInClust <- tad_conv_dt$orientation <- NULL
tad_conv_dt <- unique(tad_conv_dt)
stopifnot(!duplicated(tad_conv_dt$region))
tad_conv_dt$nConvergent[is.na(tad_conv_dt$nConvergent)] <- 0


tad_conv_dt$adjPvalComb_log10 <- -log10(tad_conv_dt$adjPvalComb)

xvar <- "nConvergent"
all_yvars <-  c("meanCorr", "meanLogFC", "adjPvalComb_log10")
yvar <- "adjPvalComb_log10"
for(yvar in all_yvars){
  outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_densplot_withCorr(xvar, yvar,tad_conv_dt )
  mtext(side=3, text = paste0("n=", ntot, "; nDS=", nDS))
  title(main = paste0(yvar, " vs. ", xvar))
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
}





ntot <- nrow(maxClust_tad_dt)
nDS <- length(unique(file.path(maxClust_tad_dt$hicds, maxClust_tad_dt$exprds)))

tmp_dt <- clust_merged_dt
tmp_dt$orientation <- tmp_dt$nConvergent <- NULL
stopifnot(is.numeric(tmp_dt$nInClust))
maxClust_tad_dt <- aggregate(nInClust~ ., data = tmp_dt, FUN=max)
colnames(maxClust_tad_dt) [colnames(maxClust_tad_dt) == "nInClust"] <- "max_nInClust"

maxClust_tad_dt$adjPvalComb_log10 <- -log10(maxClust_tad_dt$adjPvalComb)

xvar <- "max_nInClust"
all_yvars <-  c("meanCorr", "meanLogFC", "adjPvalComb_log10")
yvar <- "meanCorr"
for(yvar in all_yvars){
  outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_densplot_withCorr(xvar, yvar,maxClust_tad_dt )
  mtext(side=3, text = paste0("n=", ntot, "; nDS=", nDS))
  title(main = paste0(yvar, " vs. ", xvar))
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
}

# just for check
tmp1 <- aggregate(nInClust ~ region + hicds + exprds, data=clust_merged_dt, FUN=sum)
x1 <- setNames(tmp1$nInClust, tmp1$region)
x2 <- setNames(agg_merged_dt$CTCF_totCount, agg_merged_dt$region)
stopifnot(setequal(x1,x2))
stopifnot(setequal(names(x1),names(x2)))

ratioMax_dt <- merge(maxClust_tad_dt,agg_merged_dt, by=intersect(colnames(maxClust_tad_dt), colnames(agg_merged_dt)), all=T )
ratioMax_dt$ratioMaxClust <- ratioMax_dt$max_nInClust/ ratioMax_dt$CTCF_totCount
stopifnot(ratioMax_dt$ratioMaxClust <= 1)
  
xvar <- "ratioMaxClust"
all_yvars <-  c("meanCorr", "meanLogFC", "adjPvalComb_log10")
yvar <- "adjPvalComb_log10"
for(yvar in all_yvars){
  outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_densplot_withCorr(xvar, yvar,ratioMax_dt )
  mtext(side=3, text = paste0("n=", ntot, "; nDS=", nDS))
  title(main = paste0(yvar, " vs. ", xvar))
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
}
