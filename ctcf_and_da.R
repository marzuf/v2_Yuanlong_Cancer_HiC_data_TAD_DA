
# Rscript ctcf_and_da.R

library("readxl")
library(doMC)
library(foreach)
library(stringr)


require(ggpubr)
require(ggsci)


pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_obs_exprds <- sapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))

registerDoMC(40)
# runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878" #PIPELINE/OUTPUT_FOLDER/GM12878_40kb/TCGAluad_norm_luad/11sameNbr_runEmpPvalCombined/"
# hicds <- "GM12878_40kb"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
plotType <- "png"
myHeight <- myWidth <- 400
plotCex <- 1.2

plotTypeGG <- "svg"
ggHeight <- 6
ggWidth <- 5

fontFamily <- "Hershey"

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


plot_density <- function(p) {
  p2 <- p+  
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    )
  return(p2)
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


### RETRIEVE THE REGION WITH 0 CTCF
agg_merged_dt2 <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanCorr", "meanLogFC", "adjPvalComb")], 
                        agg_merged_dt, all.x=T, all.y=T, by=c("region", "meanCorr", "meanLogFC", "adjPvalComb"))
agg_merged_dt2$CTCF_count.forward[is.na(agg_merged_dt2$CTCF_count.forward)] <- 0
agg_merged_dt2$CTCF_count.reverse[is.na(agg_merged_dt2$CTCF_count.reverse)] <- 0
agg_merged_dt2$CTCF_totCount[is.na(agg_merged_dt2$CTCF_totCount)] <- 0
stopifnot(!is.na(agg_merged_dt2))
stopifnot(agg_merged_dt2$CTCF_count.reverse + agg_merged_dt2$CTCF_count.forward == agg_merged_dt2$CTCF_totCount)

agg_merged_dt <- agg_merged_dt2

agg_merged_dt$adjPvalComb_log10 <- -log10(agg_merged_dt$adjPvalComb)


ntot <- nrow(agg_merged_dt)
nDS <- length(unique(file.path(agg_merged_dt$hicds, agg_merged_dt$exprds)))

################### DENSPLOT CTCF count

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

signifThresh <- 0.01

################### BARPLOT CTCF count

agg_merged_dt$signif_lab <- ifelse(agg_merged_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")

agg_merged_dt$CTCF_totCount_lab <- ifelse(agg_merged_dt$CTCF_totCount == 0, "0",
                                          ifelse(agg_merged_dt$CTCF_totCount > 0 & agg_merged_dt$CTCF_totCount <= 5, ">0 & <=5",
                                                 ifelse(agg_merged_dt$CTCF_totCount > 5 & agg_merged_dt$CTCF_totCount <= 10, ">5 & <=10",
                                                        ifelse(agg_merged_dt$CTCF_totCount > 10 & agg_merged_dt$CTCF_totCount <= 15, ">10 & <=15",
                                                               ifelse(agg_merged_dt$CTCF_totCount > 15 & agg_merged_dt$CTCF_totCount <= 20, ">15 & <=20",
                                                                      ifelse(agg_merged_dt$CTCF_totCount > 20, ">20", NA))))))

stopifnot(!is.na(agg_merged_dt$CTCF_totCount_lab))

CTCFcount_levels <- c("0", ">0 & <=5", ">5 & <=10", ">10 & <=15", ">15 & <=20", ">20")

plot_dt <- aggregate(region ~ CTCF_totCount_lab+signif_lab, data=agg_merged_dt, FUN=length)
plot_dt$CTCF_totCount_lab <- factor(plot_dt$CTCF_totCount_lab, levels=rev(CTCFcount_levels))
stopifnot(!is.na(plot_dt$CTCF_totCount_lab))

tmp <- aggregate(region~signif_lab, data=plot_dt, FUN=sum)

totSignif <- setNames(tmp$region, tmp$signif_lab)

stopifnot(sum(ds_final_dt$adjPvalComb <= signifThresh) == totSignif[ "signif."])
stopifnot(sum(ds_final_dt$adjPvalComb > signifThresh) == totSignif["not signif."])

plot_dt$signif_lab <- as.character(plot_dt$signif_lab)
plot_dt$region_ratio <- plot_dt$region/totSignif[plot_dt$signif_lab]

stopifnot(!is.na(plot_dt$region_ratio))


plotTit <- "Ratio domains by # of CTCT_totCount"
subTit <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")

p <- ggbarplot(plot_dt, x="signif_lab", y="region_ratio", fill="CTCF_totCount_lab", 
                    xlab = "", ylab = "Ratio of domains")+
  scale_fill_nejm() + 
  labs(fill="CTCT_totCount") + 
  ggtitle(plotTit, subtitle=subTit)+
  theme(
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  )



outFile <- file.path(outFolder, paste0("CTCF_totCount_byTAD_bySignif_ratio_barplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))

################### DENSITY CTCFtotcount
count_dt <- agg_merged_dt
count_dt$signif_lab <- ifelse(count_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")

totSignif <- table(count_dt$signif_lab)

plot_var <- "CTCF_totCount"
plotTit <- paste0(plot_var, " dist.")
mySub <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
legTitle <- ""

p <- plot_density(ggdensity(count_dt,
                            x = plot_var,
                            y = "..density..",
                            # combine = TRUE,                  # Combine the 3 plots
                            xlab = plot_var,
                            # add = "median",                  # Add median line.
                            rug = FALSE,                      # Add marginal rug
                            color = "signif_lab",
                            fill = "signif_lab",
                            palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density"))

outFile <- file.path(outFolder, paste0(paste0(plot_var, "_signif_notSignif_dist_density."), plotTypeGG))
ggsave(p, file=outFile, height=ggHeight, width=ggWidth*1.5)
cat(paste0("... written: ", outFile, "\n"))



################### BARPLOT CTCF count by orientation

tmp <- merged_dt
tmp$signif_lab <- ifelse(tmp$adjPvalComb <= signifThresh, "signif.", "not signif.")
class_agg_dt_tmp <- aggregate(chr ~ signif_lab + Triplet_class, FUN=length, data=tmp)
colnames(class_agg_dt_tmp)[colnames(class_agg_dt_tmp) == "chr"] <- "CTCF_Count"

tmp2 <- aggregate(region~signif_lab, data=tmp, FUN=length)
totSignif <- setNames(tmp2$region, tmp2$signif_lab)

plot_dt <- class_agg_dt_tmp
plot_dt$count_ratio <- plot_dt$CTCF_Count/totSignif[plot_dt$signif_lab]
stopifnot(plot_dt$count_ratio <= 1 & plot_dt$count_ratio >= 0)
stopifnot(!is.na(plot_dt$count_ratio))

plotTit <- "Ratio motifs by class and signif."
subTit <- paste0("# in", names(totSignif), "=", totSignif, collapse="; ")

p <- ggbarplot(plot_dt, x="signif_lab", y="count_ratio", fill="Triplet_class", 
               xlab = "", ylab = "Ratio of motifs")+
  scale_fill_nejm() + 
  labs(fill="Triplet_class") + 
  ggtitle(plotTit, subtitle=subTit)+
  theme(
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  )

outFile <- file.path(outFolder, paste0("CTCF_totCount_byTripletClass_bySignif_ratio_barplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))

################### DENSITY CTCF TAD dist
dist_dt <- merged_dt
dist_dt$midPosDist <- abs(((merged_dt$start+merged_dt$end)/2) - ((merged_dt$tad_start+merged_dt$tad_end)/2))
dist_dt$midPosDist_log10 <- log10(dist_dt$midPosDist)
dist_dt$signif_lab <- ifelse(dist_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")

totSignif <- table(dist_dt$signif_lab)

plot_var <- "midPosDist_log10"
plotTit <- paste0(plot_var, " dist.")
mySub <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
legTitle <- ""

p <- plot_density(ggdensity(dist_dt,
                            x = plot_var,
                            y = "..density..",
                            # combine = TRUE,                  # Combine the 3 plots
                            xlab = plot_var,
                            # add = "median",                  # Add median line.
                            rug = FALSE,                      # Add marginal rug
                            color = "signif_lab",
                            fill = "signif_lab",
                            palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density"))

outFile <- file.path(outFolder, paste0(paste0(plot_var, "_signif_notSignif_dist_density."), plotTypeGG))
ggsave(p, file=outFile, height=ggHeight, width=ggWidth*1.5)
cat(paste0("... written: ", outFile, "\n"))

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

################### DENSPLOT nConvergent

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



################### BARPLOT nConvergent

# plot(density(tad_conv_dt$nConvergent))
# range(tad_conv_dt$nConvergent)

tad_conv_dt$signif_lab <- ifelse(tad_conv_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")

tad_conv_dt$nConvergent_lab <- ifelse(tad_conv_dt$nConvergent == 0, "0",
                                          ifelse(tad_conv_dt$nConvergent > 0 & tad_conv_dt$nConvergent <= 3, ">0 & <=3",
                                                 ifelse(tad_conv_dt$nConvergent > 3 & tad_conv_dt$nConvergent <= 10, ">3 & <=6",
                                                        ifelse(tad_conv_dt$nConvergent > 6 & tad_conv_dt$nConvergent <= 15, ">6 & <=9",
                                                               ifelse(tad_conv_dt$nConvergent > 9, ">9", NA)))))

stopifnot(!is.na(tad_conv_dt$nConvergent_lab))

nConvergent_levels <- c("0", ">0 & <=3", ">3 & <=6", ">6 & <=9", ">9")

plot_dt <- aggregate(region ~ nConvergent_lab+signif_lab, data=tad_conv_dt, FUN=length)
plot_dt$nConvergent_lab <- factor(plot_dt$nConvergent_lab, levels=rev(nConvergent_levels))
stopifnot(!is.na(plot_dt$nConvergent_lab))

tmp <- aggregate(region~signif_lab, data=plot_dt, FUN=sum)

totSignif <- setNames(tmp$region, tmp$signif_lab)

stopifnot(sum(ds_final_dt$adjPvalComb <= signifThresh) == totSignif[ "signif."])
stopifnot(sum(ds_final_dt$adjPvalComb > signifThresh) == totSignif["not signif."])

plot_dt$signif_lab <- as.character(plot_dt$signif_lab)
plot_dt$region_ratio <- plot_dt$region/totSignif[plot_dt$signif_lab]

stopifnot(!is.na(plot_dt$region_ratio))


plotTit <- "Ratio domains by # of nConvergent"
subTit <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")

p <- ggbarplot(plot_dt, x="signif_lab", y="region_ratio", fill="nConvergent_lab", 
               xlab = "", ylab = "Ratio of domains")+
  scale_fill_nejm() + 
  labs(fill="nConvergent") + 
  ggtitle(plotTit, subtitle=subTit)+
  theme(
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  )



outFile <- file.path(outFolder, paste0("nConvergent_byTAD_bySignif_ratio_barplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))


################### DENSITY signif not signif ratioMaxClust

totSignif <- table(tad_conv_dt$signif_lab)

plot_var <- "nConvergent"
plotTit <- paste0(plot_var, " dist.")
mySub <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
legTitle <- ""

p <- plot_density(ggdensity(tad_conv_dt,
                            x = plot_var,
                            y = "..density..",
                            # combine = TRUE,                  # Combine the 3 plots
                            xlab = plot_var,
                            # add = "median",                  # Add median line.
                            rug = FALSE,                      # Add marginal rug
                            color = "signif_lab",
                            fill = "signif_lab",
                            palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density"))

outFile <- file.path(outFolder, paste0(paste0(plot_var, "_signif_notSignif_dist_density."), plotTypeGG))
ggsave(p, file=outFile, height=ggHeight, width=ggWidth*1.5)
cat(paste0("... written: ", outFile, "\n"))



tmp_dt <- clust_merged_dt
tmp_dt$orientation <- tmp_dt$nConvergent <- NULL
stopifnot(is.numeric(tmp_dt$nInClust))
tmp_dt$nInClust[is.na(tmp_dt$nInClust)] <- 0
maxClust_tad_dt <- aggregate(nInClust~ ., data = tmp_dt, FUN=max)
colnames(maxClust_tad_dt) [colnames(maxClust_tad_dt) == "nInClust"] <- "max_nInClust"

ntot <- nrow(maxClust_tad_dt)
nDS <- length(unique(file.path(maxClust_tad_dt$hicds, maxClust_tad_dt$exprds)))


maxClust_tad_dt$adjPvalComb_log10 <- -log10(maxClust_tad_dt$adjPvalComb)

stopifnot(!duplicated(file.path(maxClust_tad_dt$hicds, maxClust_tad_dt$exprds, maxClust_tad_dt$region)))





################### DENSPLOT max_nInClust

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

################### BARPLOT max_nInClust

plot(density(maxClust_tad_dt$max_nInClust))
range(maxClust_tad_dt$max_nInClust)

maxClust_tad_dt$signif_lab <- ifelse(maxClust_tad_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")

maxClust_tad_dt$max_nInClust_lab <- ifelse(maxClust_tad_dt$max_nInClust == 0, "0",
                                      ifelse(maxClust_tad_dt$max_nInClust > 0 & maxClust_tad_dt$max_nInClust <= 3, ">0 & <=3",
                                             ifelse(maxClust_tad_dt$max_nInClust > 3 & maxClust_tad_dt$max_nInClust <= 10, ">3 & <=6",
                                                    ifelse(maxClust_tad_dt$max_nInClust > 6 & maxClust_tad_dt$max_nInClust <= 15, ">6 & <=9",
                                                           ifelse(maxClust_tad_dt$max_nInClust > 9, ">9", NA)))))

stopifnot(!is.na(maxClust_tad_dt$max_nInClust_lab))

maxClust_levels <- c("0", ">0 & <=3", ">3 & <=6", ">6 & <=9", ">9")

plot_dt <- aggregate(region ~ max_nInClust_lab+signif_lab, data=maxClust_tad_dt, FUN=length)
plot_dt$max_nInClust_lab <- factor(plot_dt$max_nInClust_lab, levels=rev(nConvergent_levels))
stopifnot(!is.na(plot_dt$max_nInClust_lab))

tmp <- aggregate(region~signif_lab, data=plot_dt, FUN=sum)

totSignif <- setNames(tmp$region, tmp$signif_lab)

stopifnot(sum(ds_final_dt$adjPvalComb <= signifThresh) == totSignif[ "signif."])
stopifnot(sum(ds_final_dt$adjPvalComb > signifThresh) == totSignif["not signif."])

plot_dt$signif_lab <- as.character(plot_dt$signif_lab)
plot_dt$region_ratio <- plot_dt$region/totSignif[plot_dt$signif_lab]

stopifnot(!is.na(plot_dt$region_ratio))


plotTit <- "Ratio domains by # of max_nInClust"
subTit <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")

p <- ggbarplot(plot_dt, x="signif_lab", y="region_ratio", fill="max_nInClust_lab", 
               xlab = "", ylab = "Ratio of domains")+
  scale_fill_nejm() + 
  labs(fill="max_nInClust") + 
  ggtitle(plotTit, subtitle=subTit)+
  theme(
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  )



outFile <- file.path(outFolder, paste0("max_nInClust_byTAD_bySignif_ratio_barplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))




################### DENSITY signif not signif ratioMaxClust

totSignif <- table(maxClust_tad_dt$signif_lab)

plot_var <- "max_nInClust"
plotTit <- paste0(plot_var, " dist.")
mySub <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
legTitle <- ""

p <- plot_density(ggdensity(maxClust_tad_dt,
                            x = plot_var,
                            y = "..density..",
                            # combine = TRUE,                  # Combine the 3 plots
                            xlab = plot_var,
                            # add = "median",                  # Add median line.
                            rug = FALSE,                      # Add marginal rug
                            color = "signif_lab",
                            fill = "signif_lab",
                            palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density"))

outFile <- file.path(outFolder, paste0(paste0(plot_var, "_signif_notSignif_dist_density."), plotTypeGG))
ggsave(p, file=outFile, height=ggHeight, width=ggWidth*1.5)
cat(paste0("... written: ", outFile, "\n"))







# just for check
tmp0 <- clust_merged_dt
tmp0$nInClust[is.na(tmp0$nInClust)] <- 0
tmp1 <- aggregate(nInClust ~ region + hicds + exprds, data=tmp0, FUN=sum)
x1 <- setNames(tmp1$nInClust, tmp1$region)
x2 <- setNames(agg_merged_dt$CTCF_totCount, agg_merged_dt$region)
stopifnot(setequal(x1,x2))
stopifnot(setequal(names(x1),names(x2)))

ratioMax_dt <- merge(maxClust_tad_dt,agg_merged_dt, by=intersect(colnames(maxClust_tad_dt), colnames(agg_merged_dt)), all=T )
ratioMax_dt <- ratioMax_dt[ratioMax_dt$CTCF_totCount > 0,]   ### RETAIN ONLY DOMAINS THAT HAVE CTCF
ratioMax_dt$ratioMaxClust <- ratioMax_dt$max_nInClust/ ratioMax_dt$CTCF_totCount
stopifnot(ratioMax_dt$ratioMaxClust <= 1)

ratioMax_dt$signif_lab <- ifelse(ratioMax_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")


################### DENSPLOT ratioMaxClust

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

################### DENSITY signif not signif ratioMaxClust

totSignif <- table(ratioMax_dt$signif_lab)

plot_var <- "ratioMaxClust"
plotTit <- paste0(plot_var, " dist.")
mySub <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
legTitle <- ""

p <- plot_density(ggdensity(ratioMax_dt,
                x = plot_var,
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = plot_var,
                # add = "median",                  # Add median line.
                rug = FALSE,                      # Add marginal rug
                color = "signif_lab",
                fill = "signif_lab",
                palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density"))

outFile <- file.path(outFolder, paste0(paste0(plot_var, "_signif_notSignif_dist_density."), plotTypeGG))
ggsave(p, file=outFile, height=ggHeight, width=ggWidth*1.5)
cat(paste0("... written: ", outFile, "\n"))


