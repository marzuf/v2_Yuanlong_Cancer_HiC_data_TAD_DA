
# Rscript ctcf_and_da_allDS.R

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

signifThresh <- 0.01


runFolder <- "." 


source("ctcf_da_utils.R")

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
ds_final_dt <- final_dt
ds_final_dt <- final_dt[final_dt$hicds %in% unlist(all_obs_hicds) & final_dt$exprds %in% unlist(all_obs_exprds),  ]
# ds_final_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds, ]
# stopifnot(nrow(ds_final_dt) > 0)

buildTable <- TRUE

outFolder <- file.path("CTCF_AND_DA_ALLDS")
dir.create(outFolder, recursive = TRUE)

init_ctcf_dt <- read_excel("13059_2020_2108_MOESM2_ESM.xlsx", sheet="CTCFs")
init_ctcf_dt <- as.data.frame(init_ctcf_dt)
init_ctcf_dt <- init_ctcf_dt[, 1:7]
init_ctcf_dt$chr <- as.character(init_ctcf_dt$chr)

if(buildTable){
  ctcf2tad_dt <- foreach(hicds = all_obs_hicds, .combine='rbind') %do% {
    
    cat(paste0("... start ", hicds, " \n"))
    
    tad_dt <- read.delim(file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt"), stringsAsFactors = FALSE, 
                         header=F, col.names = c("chromo", "region", "start", "end"))
    
    ### KEEP ONLY TAD REGIONS
    tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
    stopifnot(!grepl("BOUND", tad_dt$region))
    
    # assign ctcf BS to tads
    ctcf_dt <- init_ctcf_dt[init_ctcf_dt$chr %in% tad_dt$chromo,]
    stopifnot(nrow(ctcf_dt) > 0)
    stopifnot(is.numeric(ctcf_dt$start))
    stopifnot(is.numeric(ctcf_dt$end))
    stopifnot(ctcf_dt$start <= ctcf_dt$end)
    ctcf_dt$region <- NA
    ctcf_dt$hicds <- hicds
    
    i=1
    i_ctcf2tad_dt <- foreach(i = 1:nrow(ctcf_dt), .combine='rbind') %dopar% {
      
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
    i_ctcf2tad_dt
  } 
  
  outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
  save(ctcf2tad_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,"\n"))
  
}else {
  outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
  ctcf2tad_dt <- get(load(outFile))
}

# stop("-ok\n")
# ctcf2tad_dt=get(load("CTCF_AND_DA_ALLDS/ctcf2tad_dt.Rdata"))

###################### todo: TO CAHNGE KEEP HICDS + EXPRDS + KEEP THE 0s !!!


# load("CTCF_AND_DA/ENCSR444WCZ_A549_40kb/TCGAluad_mutKRAS_mutEGFR/ctcf2tad_dt.Rdata")
cat(paste0("# of CTCF BS:\t",nrow(ctcf2tad_dt), "\n"))
cat(paste0("# of CTCF BS in TADs:\t",sum(!is.na(ctcf2tad_dt$region)), "\n"))
cat(paste0("# of CTCF BS out of TADs:\t",sum(is.na(ctcf2tad_dt$region)), "\n"))

colnames(ds_final_dt)[colnames(ds_final_dt) == "start"] <- "tad_start"
colnames(ds_final_dt)[colnames(ds_final_dt) == "end"] <- "tad_end"

tmpx <- unique(file.path(ds_final_dt$hicds))
tmpy <- unique(file.path(ctcf2tad_dt$hicds))

cat(tmpy[!tmpy%in%tmpx], "\n")
cat(tmpx[!tmpx%in%tmpy], "\n")

stopifnot(setequal(file.path(ds_final_dt$hicds),file.path(ctcf2tad_dt$hicds) ))


merged_dt <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanLogFC", "meanCorr", "adjPvalComb", "tad_start", "tad_end")],
                   ctcf2tad_dt,
                   by=c("region", "hicds"), all=FALSE)
tmp <- merged_dt$region
tmp <- gsub("(.+)_.+", "\\1", tmp)
stopifnot(tmp == merged_dt$chr)

merged_dt <- merged_dt[order(merged_dt$chr, merged_dt$start, merged_dt$end ),]
# stopifnot(diff(merged_dt$start) >= 0) # not true because multiple chromo
merged_dt$region <- as.character(merged_dt$region)

#############################################################################
# presence of motifs
#############################################################################

aggByOrientation_dt <- aggregate(chr ~ orientation + hicds + exprds + region + meanCorr + meanLogFC + adjPvalComb, FUN=length, data=merged_dt)
colnames(aggByOrientation_dt)[colnames(aggByOrientation_dt) == "chr"] <- "CTCF_count"

aggByOrientation_dt$orientation_lab <- ifelse(aggByOrientation_dt$orientation == ">", "forward", 
                                              ifelse(aggByOrientation_dt$orientation == "<", "reverse", NA))
stopifnot(!is.na(aggByOrientation_dt$orientation_lab))

wide_aggByOrientation_dt <- reshape(aggByOrientation_dt[,c("hicds", "exprds","region", "orientation_lab","CTCF_count")], 
                                    idvar=c("hicds", "exprds", "region"), direction="wide", timevar = "orientation_lab")



agg_dt <- aggregate(chr ~ hicds + exprds + region + meanCorr + meanLogFC + adjPvalComb, FUN=length, data=merged_dt)
colnames(agg_dt)[colnames(agg_dt) == "chr"] <- "CTCF_totCount"

stopifnot(setequal(file.path(wide_aggByOrientation_dt$hicds, wide_aggByOrientation_dt$exprds),
                   file.path(agg_dt$hicds, agg_dt$exprds) ))

agg_merged_dt <- merge(wide_aggByOrientation_dt, agg_dt, by=c("hicds", "exprds", "region"), all=TRUE)

agg_merged_dt$CTCF_count.forward[is.na(agg_merged_dt$CTCF_count.forward)] <- 0
agg_merged_dt$CTCF_count.reverse[is.na(agg_merged_dt$CTCF_count.reverse)] <- 0
stopifnot(!is.na(agg_merged_dt))
stopifnot(agg_merged_dt$CTCF_count.reverse + agg_merged_dt$CTCF_count.forward == agg_merged_dt$CTCF_totCount)


### RETRIEVE THE REGION WITH 0 CTCF
stopifnot(setequal(file.path(ds_final_dt$hicds, ds_final_dt$exprds),
                   file.path(agg_merged_dt$hicds, agg_merged_dt$exprds) ))

agg_merged_dt2 <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanCorr", "meanLogFC", "adjPvalComb")], 
                        agg_merged_dt, all.x=T, all.y=T, by=c("hicds", "exprds", "region", "meanCorr", "meanLogFC", "adjPvalComb"))
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


################### BARPLOT CTCF count

agg_merged_dt$signif_lab <- ifelse(agg_merged_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")

agg_merged_dt$CTCF_totCount_lab <- ifelse(agg_merged_dt$CTCF_totCount == 0, "0",
                                          ifelse(agg_merged_dt$CTCF_totCount > 0 & agg_merged_dt$CTCF_totCount <= 5, ">0 & <=5",
                                                 ifelse(agg_merged_dt$CTCF_totCount > 5 & agg_merged_dt$CTCF_totCount <= 10, ">5 & <=10",
                                                        ifelse(agg_merged_dt$CTCF_totCount > 10 & agg_merged_dt$CTCF_totCount <= 15, ">10 & <=15",
                                                               ifelse(agg_merged_dt$CTCF_totCount > 15 & agg_merged_dt$CTCF_totCount <= 20, ">15 & <=20",
                                                                      ifelse(agg_merged_dt$CTCF_totCount > 20, ">20", NA))))))

stopifnot(!is.na(agg_merged_dt$CTCF_totCount_lab))
stopifnot(sum(table(agg_merged_dt$CTCF_totCount_lab)) == nrow(agg_merged_dt))

CTCF_totCount_breaks <- c(0,5,10,15,20)
# cat(paste0("... running get_fract_lab2\n"))
check_labs <- get_fract_lab2(vect_values=agg_merged_dt$CTCF_totCount, range_levels = CTCF_totCount_breaks)
stopifnot(gsub("<=0", "0",check_labs) == agg_merged_dt$CTCF_totCount_lab)
# cat(paste0("... running get_fract_lab0\n"))
# check_labs <- get_fract_lab0(vect_values=agg_merged_dt$CTCF_totCount, range_levels = CTCF_totCount_breaks)
# stopifnot(gsub("<=0", "0",check_labs) == agg_merged_dt$CTCF_totCount_lab)
# # CTCFcount_levels <- c("0", ">0 & <=5", ">5 & <=10", ">10 & <=15", ">15 & <=20", ">20")
CTCFcount_levels <- gsub("<=0", "0",get_level_labs(CTCF_totCount_breaks))

plot_dt <- aggregate(region ~ CTCF_totCount_lab+signif_lab, data=agg_merged_dt, FUN=length)
plot_dt$CTCF_totCount_lab <- factor(plot_dt$CTCF_totCount_lab, levels=rev(CTCFcount_levels))
stopifnot(!is.na(plot_dt$CTCF_totCount_lab))

# stop("---ok\n")

tmp <- aggregate(region~signif_lab, data=plot_dt, FUN=sum)

totSignif <- setNames(tmp$region, tmp$signif_lab)

stopifnot(sum(ds_final_dt$adjPvalComb <= signifThresh) == totSignif[ "signif."])
stopifnot(sum(ds_final_dt$adjPvalComb > signifThresh) == totSignif["not signif."])

plot_dt$signif_lab <- as.character(plot_dt$signif_lab)
plot_dt$region_ratio <- plot_dt$region/totSignif[plot_dt$signif_lab]

stopifnot(!is.na(plot_dt$region_ratio))


plotTit <- "Ratio domains by # of CTCF_totCount"
subTit <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")

p <- ggbarplot(plot_dt, x="signif_lab", y="region_ratio", fill="CTCF_totCount_lab", 
                    xlab = "", ylab = "Ratio of domains")+
  scale_fill_nejm() + 
  labs(fill="CTCF_totCount") + 
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
# look by cluster                                       ##### use clustByTAD_dt
#############################################################################

merged_dt$region_id <- file.path(merged_dt$hicds, merged_dt$exprds, merged_dt$region)

if(buildTable) {
  
  
  
  clustByTAD_dt <- foreach(region = unique(merged_dt$region_id), .combine='rbind') %dopar% {
    
    tad_agg_dt <- merged_dt[merged_dt$region_id == region,]  
    stopifnot(diff(tad_agg_dt$start) >= 0) 
    
    nConvergent <- str_count(paste0(tad_agg_dt$orientation, collapse=""), "><")
    
    data.frame(
      hicds = dirname(dirname(region)),
      exprds=basename(dirname(region)),
      region  = basename(region),
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

stopifnot(setequal(file.path(ds_final_dt$hicds, ds_final_dt$exprds),
                   file.path(clustByTAD_dt$hicds, clustByTAD_dt$exprds) ))

clust_merged_dt <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanLogFC", "meanCorr", "adjPvalComb", "tad_start", "tad_end")],
                   clustByTAD_dt,
                   by=c("hicds", "exprds", "region"), all.x=TRUE, all.y=FALSE)
clust_merged_dt$chr <- gsub("(.+)_.+", "\\1", as.character(clust_merged_dt$region))
# stopifnot(clust_merged_dt$chr %in% tad_dt$chromo)
stopifnot(clust_merged_dt$chr %in% paste0("chr", 1:22))


#############************************************************************
#############*#############************************************************************ nConvergent
#############************************************************************



tad_conv_dt <- clust_merged_dt
tad_conv_dt$nInClust <- tad_conv_dt$orientation <- NULL
tad_conv_dt <- unique(tad_conv_dt)
stopifnot(!duplicated(file.path(tad_conv_dt$hicds, tad_conv_dt$exprds, tad_conv_dt$region)))
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
                                                 ifelse(tad_conv_dt$nConvergent > 3 & tad_conv_dt$nConvergent <= 6, ">3 & <=6",
                                                        ifelse(tad_conv_dt$nConvergent > 6 & tad_conv_dt$nConvergent <= 9, ">6 & <=9",
                                                               ifelse(tad_conv_dt$nConvergent > 9, ">9", NA)))))

stopifnot(!is.na(tad_conv_dt$nConvergent_lab))
stopifnot(sum(table(tad_conv_dt$nConvergent_lab)) == nrow(tad_conv_dt))

nConvergent_breaks <- c(0,3,6,9)
# cat(paste0("... running get_fract_lab2\n"))
check_labs <- get_fract_lab2(vect_values=tad_conv_dt$nConvergent, range_levels = nConvergent_breaks)
stopifnot(gsub("<=0", "0",check_labs) == tad_conv_dt$nConvergent_lab)
# cat(paste0("... running get_fract_lab0\n"))
# check_labs <- get_fract_lab0(vect_values=tad_conv_dt$nConvergent, range_levels = nConvergent_breaks)
# stopifnot(gsub("<=0", "0",check_labs) == tad_conv_dt$nConvergent_lab)
# # nConvergent_levels <- c("0", ">0 & <=3", ">3 & <=6", ">6 & <=9", ">9")
nConvergent_levels <- gsub("<=0", "0",get_level_labs(nConvergent_breaks))


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


################### DENSITY signif not signif nConvergent 

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





#############************************************************************
#############*#############************************************************************ Singletons
#############************************************************************
#############*
##### singleton: ratio

tad_ratioSing_dt <- aggregate(nInClust ~ hicds + exprds + region+meanLogFC+meanCorr+adjPvalComb, data =clust_merged_dt, FUN=function(x) mean(x == 1))
colnames(tad_ratioSing_dt)[colnames(tad_ratioSing_dt) == "nInClust"] <- "ratioSingleton"
stopifnot(!duplicated(file.path(tad_ratioSing_dt$hicds, tad_ratioSing_dt$exprds, tad_ratioSing_dt$region)))
stopifnot(!is.na(tad_ratioSing_dt$nSingleton))

tad_ratioSing_dt$adjPvalComb_log10 <- -log10(tad_ratioSing_dt$adjPvalComb)

################### DENSPLOT ratioSingletons

save(tad_ratioSing_dt, file="tad_ratioSing_dt.Rdata", version=2)

xvar <- "ratioSingleton"
all_yvars <-  c("meanCorr", "meanLogFC", "adjPvalComb_log10")
yvar <- "adjPvalComb_log10"
for(yvar in all_yvars){
  outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_densplot_withCorr(xvar, yvar,tad_ratioSing_dt )
  mtext(side=3, text = paste0("n=", ntot, "; nDS=", nDS))
  title(main = paste0(yvar, " vs. ", xvar))
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
}

################### DENSITY signif not signif ratioSingleton 
tad_ratioSing_dt$signif_lab <- ifelse(tad_ratioSing_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")
totSignif <- table(tad_ratioSing_dt$signif_lab)

plot_var <- "ratioSingleton"
plotTit <- paste0(plot_var, " dist.")
mySub <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
legTitle <- ""

p <- plot_density(ggdensity(tad_ratioSing_dt,
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






##### singleton: number 

# clust_merged_dt <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanLogFC", "meanCorr", "adjPvalComb", "tad_start", "tad_end")],
#                          clustByTAD_dt,
#                          by=c("hicds", "exprds", "region"), all.x=TRUE, all.y=FALSE)
# clust_merged_dt$chr <- gsub("(.+)_.+", "\\1", as.character(clust_merged_dt$region))
# # stopifnot(clust_merged_dt$chr %in% tad_dt$chromo)
# stopifnot(clust_merged_dt$chr %in% paste0("chr", 1:22))


tad_nbrSing_dt <- clust_merged_dt
tad_nbrSing_dt$nInClust[is.na(tad_nbrSing_dt$nInClust)] <- 0  # for the actual number -> replace number by 0 [not done for ratio !!!]
tad_nbrSing_dt <-  aggregate(nInClust ~ hicds + exprds + region+meanLogFC+meanCorr+adjPvalComb, data =tad_nbrSing_dt, FUN=function(x) sum(x == 1))
colnames(tad_nbrSing_dt)[colnames(tad_nbrSing_dt) == "nInClust"] <- "nSingleton"
stopifnot(!duplicated(file.path(tad_nbrSing_dt$hicds, tad_nbrSing_dt$exprds, tad_nbrSing_dt$region)))
stopifnot(!is.na(tad_nbrSing_dt$nSingleton))

tad_nbrSing_dt$adjPvalComb_log10 <- -log10(tad_nbrSing_dt$adjPvalComb)

stopifnot(setequal(file.path(tad_nbrSing_dt$hicds, tad_nbrSing_dt$exprds, tad_nbrSing_dt$region),
file.path(ds_final_dt$hicds, ds_final_dt$exprds, ds_final_dt$region)))

################### DENSPLOT nSingletons

xvar <- "nSingleton"
all_yvars <-  c("meanCorr", "meanLogFC", "adjPvalComb_log10")
yvar <- "adjPvalComb_log10"
for(yvar in all_yvars){
  outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_densplot_withCorr(xvar, yvar,tad_nbrSing_dt )
  mtext(side=3, text = paste0("n=", ntot, "; nDS=", nDS))
  title(main = paste0(yvar, " vs. ", xvar))
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
}



################### BARPLOT nSingleton

# plot(density(tad_nbrSing_dt$nSingleton))
# range(tad_nbrSing_dt$nSingleton) # 0-21

tad_nbrSing_dt$signif_lab <- ifelse(tad_nbrSing_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")

tad_nbrSing_dt$nSingleton_lab <- ifelse(tad_nbrSing_dt$nSingleton == 0, "0",
                                      ifelse(tad_nbrSing_dt$nSingleton > 0 & tad_nbrSing_dt$nSingleton <= 5, ">0 & <=5",
                                             ifelse(tad_nbrSing_dt$nSingleton > 5 & tad_nbrSing_dt$nSingleton <= 10, ">5 & <=10",
                                                    ifelse(tad_nbrSing_dt$nSingleton > 10 & tad_nbrSing_dt$nSingleton <= 15, ">10 & <=15",
                                                           ifelse(tad_nbrSing_dt$nSingleton > 10, ">15", NA)))))

stopifnot(!is.na(tad_nbrSing_dt$nSingleton_lab))
stopifnot(sum(table(tad_nbrSing_dt$nSingleton_lab)) == nrow(tad_nbrSing_dt))

nSingleton_breaks <- c(0,5,10,15)
# cat(paste0("... running get_fract_lab2\n"))
check_labs <- get_fract_lab2(vect_values=tad_nbrSing_dt$nSingleton, range_levels = nSingleton_breaks)
stopifnot(gsub("<=0", "0",check_labs) == tad_nbrSing_dt$nSingleton_lab)
# cat(paste0("... running get_fract_lab0\n"))
# check_labs <- get_fract_lab0(vect_values=tad_nbrSing_dt$nSingleton, range_levels = nSingleton_breaks)
# stopifnot(gsub("<=0", "0",check_labs) == tad_nbrSing_dt$nSingleton_lab)
# nSingleton_levels <- c("0", ">0 & <=5", ">5 & <=10", ">10 & <=15", ">15")
nSingleton_levels <- gsub("<=0", "0",get_level_labs(nSingleton_breaks))


plot_dt <- aggregate(region ~ nSingleton_lab+signif_lab, data=tad_nbrSing_dt, FUN=length)
plot_dt$nSingleton_lab <- factor(plot_dt$nSingleton_lab, levels=rev(nSingleton_levels))
stopifnot(!is.na(plot_dt$nSingleton_lab))

tmp <- aggregate(region~signif_lab, data=plot_dt, FUN=sum)

totSignif <- setNames(tmp$region, tmp$signif_lab)

save(tad_nbrSing_dt, file="tad_nbrSing_dt.Rdata", version=2)
stopifnot(sum(ds_final_dt$adjPvalComb <= signifThresh) == totSignif[ "signif."])  #### NEEEED TO CHECK !!!!!!!!!!!!!!!!!!!
stopifnot(sum(ds_final_dt$adjPvalComb > signifThresh) == totSignif["not signif."])

plot_dt$signif_lab <- as.character(plot_dt$signif_lab)
plot_dt$region_ratio <- plot_dt$region/totSignif[plot_dt$signif_lab]

stopifnot(!is.na(plot_dt$region_ratio))


plotTit <- "Ratio domains by # of nSingleton"
subTit <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")

p <- ggbarplot(plot_dt, x="signif_lab", y="region_ratio", fill="nSingleton_lab", 
               xlab = "", ylab = "Ratio of domains")+
  scale_fill_nejm() + 
  labs(fill="nSingleton") + 
  ggtitle(plotTit, subtitle=subTit)+
  theme(
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  )

outFile <- file.path(outFolder, paste0("nSingleton_byTAD_bySignif_ratio_barplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))

################### DENSITY signif not signif nbrSingleton 
tad_nbrSing_dt$signif_lab <- ifelse(tad_nbrSing_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")
totSignif <- table(tad_nbrSing_dt$signif_lab)

plot_var <- "nSingleton"
plotTit <- paste0(plot_var, " dist.")
mySub <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
legTitle <- ""

p <- plot_density(ggdensity(tad_nbrSing_dt,
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





















#############************************************************************
#############*#############************************************************************ max/mean InClust
#############************************************************************




tmp_dt <- clust_merged_dt
tmp_dt$orientation <- tmp_dt$nConvergent <- NULL
stopifnot(is.numeric(tmp_dt$nInClust))
tmp_dt$nInClust[is.na(tmp_dt$nInClust)] <- 0

all_agg_funs <- c("max", "mean")  # >>>> do for maxInClust and meanInClust

for(agg_fun in all_agg_funs) {
  
  funClust_tad_dt <- aggregate(nInClust~ ., data = tmp_dt, FUN=agg_fun)
  colnames(funClust_tad_dt) [colnames(funClust_tad_dt) == "nInClust"] <- paste0(agg_fun, "_nInClust")
  
  ntot <- nrow(funClust_tad_dt)
  nDS <- length(unique(file.path(funClust_tad_dt$hicds, funClust_tad_dt$exprds)))
  
  
  funClust_tad_dt$adjPvalComb_log10 <- -log10(funClust_tad_dt$adjPvalComb)
  
  stopifnot(!duplicated(file.path(funClust_tad_dt$hicds, funClust_tad_dt$exprds, funClust_tad_dt$region)))
  
  ################### DENSPLOT max_nInClust
  
  xvar <- paste0(agg_fun, "_nInClust")
  all_yvars <-  c("meanCorr", "meanLogFC", "adjPvalComb_log10")
  yvar <- "meanCorr"
  for(yvar in all_yvars){
    outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "densplot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    do_densplot_withCorr(xvar, yvar,funClust_tad_dt )
    mtext(side=3, text = paste0("n=", ntot, "; nDS=", nDS))
    title(main = paste0(yvar, " vs. ", xvar))
    foo <- dev.off()
    cat(paste0("... written: ", outFile,"\n"))
  }
  
  ################### BARPLOT max_nInClust
  
  plot(density(funClust_tad_dt[,paste0(agg_fun, "_nInClust")]))
  range(funClust_tad_dt[,paste0(agg_fun, "_nInClust")])
  
  cat(paste0("range(funClust_tad_dt[",agg_fun, "_nInClust)]", "\t=", range(funClust_tad_dt[,paste0(agg_fun, "_nInClust")]), "\n"))

  if(agg_fun=="mean") save(funClust_tad_dt, file="funClust_tad_dt.Rdata", version=2)
    
  funClust_tad_dt$signif_lab <- ifelse(funClust_tad_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")
  
  funClust_tad_dt[,paste0(agg_fun, "_nInClust_lab")] <- ifelse(funClust_tad_dt[,paste0(agg_fun, "_nInClust")] == 0, "0",
                                             ifelse(funClust_tad_dt[,paste0(agg_fun, "_nInClust")] > 0 & funClust_tad_dt[,paste0(agg_fun, "_nInClust")] <= 3, ">0 & <=3",
                                                    ifelse(funClust_tad_dt[,paste0(agg_fun, "_nInClust")] > 3 & funClust_tad_dt[,paste0(agg_fun, "_nInClust")] <= 6, ">3 & <=6",
                                                           ifelse(funClust_tad_dt[,paste0(agg_fun, "_nInClust")] > 6 & funClust_tad_dt[,paste0(agg_fun, "_nInClust")] <= 9, ">6 & <=9",
                                                                  ifelse(funClust_tad_dt[,paste0(agg_fun, "_nInClust")] > 9, ">9", NA)))))
  stopifnot(!is.na(funClust_tad_dt[,paste0(agg_fun, "_nInClust_lab")]))
  stopifnot(sum(table(funClust_tad_dt[,paste0(agg_fun, "_nInClust_lab")])) == nrow(funClust_tad_dt))
  
  
  maxClust_breaks <- c(0,3,6,9) # same for inClust
  cat(paste0("... running get_fract_lab2\n"))
  check_labs <- get_fract_lab2(vect_values=funClust_tad_dt[,paste0(agg_fun, "_nInClust")] , range_levels = maxClust_breaks)
  stopifnot(gsub("<=0", "0",check_labs) ==   funClust_tad_dt[,paste0(agg_fun, "_nInClust_lab")] )
  # cat(paste0("... running get_fract_lab0\n"))
  # check_labs <- get_fract_lab0(vect_values=agg_merged_dt$CTCF_totCount, range_levels = CTCF_totCount_breaks)
  # stopifnot(gsub("<=0", "0",check_labs) == agg_merged_dt$CTCF_totCount_lab)
  #   maxClust_levels <- c("0", ">0 & <=3", ">3 & <=6", ">6 & <=9", ">9")
  maxClust_levels <- gsub("<=0", "0",get_level_labs(maxClust_breaks))
  
  
  plot_dt <- aggregate(as.formula(paste0("region ~ ", agg_fun, "_nInClust_lab + signif_lab")), data=funClust_tad_dt, FUN=length)
  plot_dt[,paste0(agg_fun, "_nInClust_lab")] <- factor(plot_dt[,paste0(agg_fun, "_nInClust_lab")], levels=rev(maxClust_levels))
  stopifnot(!is.na(plot_dt[,paste0(agg_fun, "_nInClust_lab")]))
  
  tmp <- aggregate(region~signif_lab, data=plot_dt, FUN=sum)
  
  totSignif <- setNames(tmp$region, tmp$signif_lab)
  
  stopifnot(sum(ds_final_dt$adjPvalComb <= signifThresh) == totSignif[ "signif."])
  stopifnot(sum(ds_final_dt$adjPvalComb > signifThresh) == totSignif["not signif."])
  
  plot_dt$signif_lab <- as.character(plot_dt$signif_lab)
  plot_dt$region_ratio <- plot_dt$region/totSignif[plot_dt$signif_lab]
  
  stopifnot(!is.na(plot_dt$region_ratio))
  
  
  plotTit <- paste0("Ratio domains by # of ", agg_fun, "_nInClust")
  subTit <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
  
  p <- ggbarplot(plot_dt, x="signif_lab", y="region_ratio", fill=paste0(agg_fun, "_nInClust_lab"), 
                 xlab = "", ylab = "Ratio of domains")+
    scale_fill_nejm() + 
    labs(fill=paste0(agg_fun, "_nInClust")) + 
    ggtitle(plotTit, subtitle=subTit)+
    theme(
      plot.title = element_text(size=16, face = "bold", hjust=0.5),
      plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
    )
  
  
  
  outFile <- file.path(outFolder, paste0(agg_fun, "_nInClust_byTAD_bySignif_ratio_barplot.", plotTypeGG))
  ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  
  
  
  ################### DENSITY signif not signif ratioMaxClust
  
  totSignif <- table(funClust_tad_dt$signif_lab)
  
  plot_var <- paste0(agg_fun, "_nInClust")
  plotTit <- paste0(plot_var, " dist.")
  mySub <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
  legTitle <- ""
  
  p <- plot_density(ggdensity(funClust_tad_dt,
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
  
  if(agg_fun == "max") 
    maxClust_tad_dt <- funClust_tad_dt
  
}








# just for check
tmp0 <- clust_merged_dt
tmp0$nInClust[is.na(tmp0$nInClust)] <- 0
tmp1 <- aggregate(nInClust ~ region + hicds + exprds, data=tmp0, FUN=sum)
x1 <- setNames(tmp1$nInClust, tmp1$region)
x2 <- setNames(agg_merged_dt$CTCF_totCount, agg_merged_dt$region)
stopifnot(setequal(x1,x2))
stopifnot(setequal(names(x1),names(x2)))

stopifnot(setequal(file.path(maxClust_tad_dt$hicds, maxClust_tad_dt$exprds),
                   file.path(agg_merged_dt$hicds, agg_merged_dt$exprds) ))


#############************************************************************
#############*#############************************************************************ ratioMaxClust
#############************************************************************



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



#############************************************************************
#############*#############************************************************************ scores
#############************************************************************

current_score <- c("MotifScore")
agg_fun <- "sum"

all_scores <- c("MotifScore", "ChipSeqScore")
all_agg_funs <- c("sum", "mean", "max")


merged_dt <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanLogFC", "meanCorr", "adjPvalComb", "tad_start", "tad_end")],
                   ctcf2tad_dt,
                   by=c("region", "hicds"), all=FALSE)
tmp <- merged_dt$region
tmp <- gsub("(.+)_.+", "\\1", tmp)
stopifnot(tmp == merged_dt$chr)

merged_dt <- merged_dt[order(merged_dt$chr, merged_dt$start, merged_dt$end ),]
# stopifnot(diff(merged_dt$start) >= 0) # not true because multiple chromo
merged_dt$region <- as.character(merged_dt$region)


for(current_score in all_scores) {
  
  
  
  
  
  for(agg_fun in all_agg_funs) {
    newCol <- paste0(agg_fun, "_", current_score)
    
    agg_score_dt <- aggregate(as.formula(paste0(current_score, " ~ region + hicds+exprds+meanLogFC+meanCorr+adjPvalComb")), data=merged_dt, FUN=agg_fun)
    agg_score_dt$adjPvalComb_log10 <- -log10(agg_score_dt$adjPvalComb)
    colnames(agg_score_dt)[ colnames(agg_score_dt) == current_score] <- newCol
    
    xvar <- newCol
    all_yvars <-  c("meanCorr", "meanLogFC", "adjPvalComb_log10")
    yvar <- "adjPvalComb_log10"
    for(yvar in all_yvars){
      outFile <- file.path(outFolder, paste0(xvar, "_vs_", yvar, "densplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      do_densplot_withCorr(xvar, yvar,agg_score_dt )
      mtext(side=3, text = paste0("n=", ntot, "; nDS=", nDS))
      title(main = paste0(yvar, " vs. ", xvar))
      foo <- dev.off()
      cat(paste0("... written: ", outFile,"\n"))
    }
    
    
    ################### DENSITY signif not signif score aggfun
    
    agg_score_dt$signif_lab <- ifelse(agg_score_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")
    
    totSignif <- table(agg_score_dt$signif_lab)
    
    plot_var <- newCol
    plotTit <- paste0(plot_var, " dist.")
    mySub <- paste0("#", names(totSignif), "=", totSignif, collapse="; ")
    legTitle <- ""
    
    p <- plot_density(ggdensity(agg_score_dt,
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
    
    
  }
}

