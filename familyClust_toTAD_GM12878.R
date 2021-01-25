# Rscript familyClust_toTAD_GM12878.R

require(doMC)
require(foreach)
registerDoMC(40)
require(ggpubr)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


# famClust_dt <- get(load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV2/ctcf4000_family260000/all_fams_agg_cluster_dt.Rdata"))
# famClust_dt2 <- get(load("CTCF_AND_FAMILY_CLUSTERS/ctcf4000_family260000/all_fams_agg_cluster_dt.Rdata"))  
# all.equal(famClust_dt,famClust_dt2)
# # TRUE
inFolder <- "CTCF_AND_FAMILY_CLUSTERS_RANDOMV2"
d_famClust <- 260000
d_ctcfClust <- 10000

outFolder <- file.path("FAMILYCLUST_TOTAD_GM12878",  paste0("ctcf", d_ctcfClust, "_family", d_famClust))
dir.create(outFolder, recursive = TRUE)

gm_tad_dt <- read.delim(file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878/GM12878_40kb/genes2tad/all_assigned_regions.txt"),
                        col.names = c("chromo", "region", "start", "end"), stringsAsFactors = FALSE, header=F)

signifThresh <- 0.01
signif_gm_dt <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878/PSEUDO_FINAL_TABLE/all_result_dt.Rdata"))
signif_gm_dt$region <- as.character(signif_gm_dt$region)
signif_tads <- signif_gm_dt$region[signif_gm_dt$adjPvalComb <= signifThresh]
                    
ctcf_clust_dt <- get(load(file.path(inFolder, paste0("ctcf", d_ctcfClust, "_family", d_famClust), "agg_ctcf_cluster_dt.Rdata")))
stopifnot(!duplicated(ctcf_clust_dt$clusterID))
all_ctcf_clustSize <- setNames(ctcf_clust_dt$clusterSize, paste0("CTCF_", ctcf_clust_dt$clusterID))
all_ctcf_clustLength <- setNames(ctcf_clust_dt$maxPos-ctcf_clust_dt$minPos+1, paste0("CTCF_", ctcf_clust_dt$clusterID))
stopifnot( all_ctcf_clustLength >= 0)

famClust_dt <- get(load(file.path(inFolder, paste0("ctcf", d_ctcfClust, "_family", d_famClust), "all_fams_agg_cluster_dt.Rdata")))

obsAndPerm_dt <- get(load(file.path(inFolder, paste0("ctcf", d_ctcfClust, "_family", d_famClust), 
                                    "obsAndPermut_closest_clusters_dt.Rdata")))
all_obs_dt <- do.call('rbind', lapply(obsAndPerm_dt, function(x) x[["obs_final_fams_cluster_dt"]]))
  # this one as only cluster size > 1

stopifnot(is.numeric( gm_tad_dt$end))
stopifnot(is.numeric( gm_tad_dt$start))

gm_tad_dt <- gm_tad_dt[order(gm_tad_dt$chromo, gm_tad_dt$start, gm_tad_dt$end),]
gm_tad_dt$chromo <- as.character(gm_tad_dt$chromo)

famClust_dt$chromo <- as.character(famClust_dt$chromo)
famClust_dt <- famClust_dt[famClust_dt$chromo %in% gm_tad_dt$chromo,]
stopifnot(nrow(famClust_dt) > 0)

famClust_dt$nTADs <- foreach(i  = 1:nrow(famClust_dt), .combine='c') %dopar% {
  i_min <- which(gm_tad_dt$chromo == famClust_dt$chromo[i] & 
                   gm_tad_dt$start <= famClust_dt$minPos[i] & gm_tad_dt$end >= famClust_dt$minPos[i])
  stopifnot(length(i_min) == 1)
  i_max <- which(gm_tad_dt$chromo == famClust_dt$chromo[i] & 
                   gm_tad_dt$start <= famClust_dt$maxPos[i] & gm_tad_dt$end >= famClust_dt$maxPos[i])
  stopifnot(length(i_max) == 1)
  nTADs <- i_max - i_min + 1
  
  stopifnot(nTADs >= 1)
  nTADs
}
famClust_dt$nTADs_ratio <- famClust_dt$nTADs/famClust_dt$clusterSize

famClust_dt$nSignifTADs <- foreach(i  = 1:nrow(famClust_dt), .combine='c') %dopar% {
  i_min <- which(gm_tad_dt$chromo == famClust_dt$chromo[i] & 
                   gm_tad_dt$start <= famClust_dt$minPos[i] & gm_tad_dt$end >= famClust_dt$minPos[i])
  stopifnot(length(i_min) == 1)
  i_max <- which(gm_tad_dt$chromo == famClust_dt$chromo[i] & 
                   gm_tad_dt$start <= famClust_dt$maxPos[i] & gm_tad_dt$end >= famClust_dt$maxPos[i])
  stopifnot(length(i_max) == 1)
  nTADs <- i_max - i_min + 1
  stopifnot(nTADs >= 1)
  sum(gm_tad_dt$region[i_min:i_max] %in% signif_tads)
}


all_obs_dt$chromo <- as.character(all_obs_dt$chromo)
all_obs_dt <- all_obs_dt[all_obs_dt$chromo %in% gm_tad_dt$chromo,]
stopifnot(nrow(all_obs_dt) > 0)


all_obs_dt$nTADs <- foreach(i  = 1:nrow(all_obs_dt), .combine='c') %dopar% {
  i_min <- which(gm_tad_dt$chromo == all_obs_dt$chromo[i] & 
                   gm_tad_dt$start <= all_obs_dt$minPos[i] & gm_tad_dt$end >= all_obs_dt$minPos[i])
  stopifnot(length(i_min) == 1)
  i_max <- which(gm_tad_dt$chromo == all_obs_dt$chromo[i] & 
                   gm_tad_dt$start <= all_obs_dt$maxPos[i] & gm_tad_dt$end >= all_obs_dt$maxPos[i])
  stopifnot(length(i_max) == 1)
  nTADs <- i_max - i_min + 1
  stopifnot(nTADs >= 1)
  nTADs
}
all_obs_dt$nTADs_ratio <- all_obs_dt$nTADs/all_obs_dt$clusterSize

all_obs_dt$nSignifTADs <- foreach(i  = 1:nrow(all_obs_dt), .combine='c') %dopar% {
  i_min <- which(gm_tad_dt$chromo == all_obs_dt$chromo[i] & 
                   gm_tad_dt$start <= all_obs_dt$minPos[i] & gm_tad_dt$end >= all_obs_dt$minPos[i])
  stopifnot(length(i_min) == 1)
  i_max <- which(gm_tad_dt$chromo == all_obs_dt$chromo[i] & 
                   gm_tad_dt$start <= all_obs_dt$maxPos[i] & gm_tad_dt$end >= all_obs_dt$maxPos[i])
  stopifnot(length(i_max) == 1)
  nTADs <- i_max - i_min + 1
  stopifnot(nTADs >= 1)
  sum(gm_tad_dt$region[i_min:i_max] %in% signif_tads)
  
}

stopifnot(all_obs_dt$closest_CTCFclust %in% names(all_ctcf_clustLength))
stopifnot(all_obs_dt$closest_CTCFclust %in% names(all_ctcf_clustSize))

all_obs_dt$closest_CTCFclust_length <- all_ctcf_clustLength[as.character(all_obs_dt$closest_CTCFclust)]
stopifnot(!is.na(all_obs_dt$closest_CTCFclust_length))

all_obs_dt$closest_CTCFclust_size <- all_ctcf_clustSize[as.character(all_obs_dt$closest_CTCFclust)]
stopifnot(!is.na(all_obs_dt$closest_CTCFclust_size))


myx <- all_obs_dt$closest_CTCFclust_size
myy <- all_obs_dt$nTADs_ratio
outFile <- file.path(outFolder, "nTADsRatio_by_closestClustSize_densplot.png")
png(outFile, height=400, width=400)
densplot(x=myx,
         y=myy,
         xlab ="closestCTCFclust_size",
         ylab = "nTADs_ratio")
mtext(side=3, text=paste0("GM12878 TADs; -d = ", d_famClust))
foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))

myx <- all_obs_dt$closest_CTCFclust_size
myy <- all_obs_dt$nTADs
outFile <- file.path(outFolder, "nTADs_by_closestClustSize_densplot.png")
png(outFile, height=400, width=400)
densplot(x=myx,
         y=myy,
         xlab ="closestCTCFclust_size",
         ylab = "nTADs")
mtext(side=3, text=paste0("GM12878 TADs; -d = ", d_famClust))
foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))

myx <- log10(all_obs_dt$closest_CTCFclust_length)
myy <- all_obs_dt$nTADs_ratio
outFile <- file.path(outFolder, "nTADsRatio_by_closestClustLengthLog10_densplot.png")
png(outFile, height=400, width=400)
densplot(x=myx,
         y=myy,
         xlab ="closestCTCFclust_length [log10]",
         ylab = "nTADs_ratio")
mtext(side=3, text=paste0("GM12878 TADs; -d = ", d_famClust))
foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))

outFile <- file.path(outFolder, "nTADs_by_famClustSize_boxplot.png")
png(outFile, height=400, width=400)
boxplot(nTADs~clusterSize, data=all_obs_dt,
        main=paste0("# of TADs in family cluster"))
mtext(side=3, text=paste0("GM12878 TADs; -d = ", d_famClust))
foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))


myx <- log10(all_obs_dt$minDist_CTCFclust)
myy <- all_obs_dt$nTADs_ratio
outFile <- file.path(outFolder, "nTADsRatio_by_minDistLog10_densplot.png")
png(outFile, height=400, width=400)
densplot(x=myx,
  y=myy,
  xlab ="minDist_CTCFclust [log10]",
  ylab = "nTADs_ratio")
mtext(side=3, text=paste0("GM12878 TADs; -d = ", d_famClust))
foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))

outFile <- file.path(outFolder, "nTADs_by_famClustSize_boxplot.png")
png(outFile, height=400, width=400)
boxplot(nTADs~clusterSize, data=all_obs_dt,
        main=paste0("# of TADs in family cluster"))
mtext(side=3, text=paste0("-d = ", d_famClust))
foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))


p <- gghistogram(all_obs_dt, x = "clusterSize", fill = "lightgray", binwidth=1, boundary=0,
                 add = "mean", rug = F) +
  ggtitle("Distribution family cluster size", subtitle=paste0("GM12878 TADs; -d = ", d_famClust)) +
  labs(x="cluster size", y="count")+
  scale_y_continuous(breaks=scales::pretty_breaks(n = 10)) + 
  # scale_x_continuous(breaks=scales::pretty_breaks(n = 10)) 
  scale_x_continuous(breaks=sort(unique(all_obs_dt$clustSize)),
                     labels=sort(unique(all_obs_dt$clustSize)))+
  # scale_x_continuous(breaks=1:10,labels=sort(unique(aggSize_clust_dt$clustSize))) + 
  theme(
    plot.subtitle = element_text(size=14, face="italic", hjust=0.5),
    plot.title = element_text(size=16, face="bold", hjust=0.5),
    axis.text = element_text(size=12)
  )
#> Warning: Using `bins = 30` by default. Pick better value with the argument `bins`.

outFile <- file.path(outFolder, paste0("family_clustSize_dist_d-", d_famClust, ".svg"))
ggsave(p, filename = outFile, height=5, width=6)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("all_obs_dt.Rdata"))
save(all_obs_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

