# Rscript look_inside_clusters.R


ctcf_clustSize <- 10000
# ctcf_clustSize <- 10000
# ctcf_clustSize <- 20000
family_clustSize <- 260000
# family_clustSize <- 130000
outFolder <- file.path("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3", paste0("ctcf", ctcf_clustSize, "_family", family_clustSize))

# load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3/ctcf4000_family260000/obsAndPermut_closest_clusters_dt.Rdata")
obsAndPermut_closest_clusters_dt <- get(load(file.path(outFolder, "obsAndPermut_closest_clusters_dt.Rdata")))
all_obs_dt <- do.call('rbind', lapply(obsAndPermut_closest_clusters_dt, function(x) x[["obs_final_fams_cluster_dt"]]))

# load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3/ctcf4000_family260000/agg_ctcf_cluster_dt.Rdata")
agg_ctcf_cluster_dt <- get(load(file.path(outFolder, "agg_ctcf_cluster_dt.Rdata")))
ctcf_clustFile <- file.path(outFolder, paste0("ctcf_bs_nanni2020_cluster", ctcf_clustSize, "bp.bed"))
ctcf_cluster_dt <- read.delim(ctcf_clustFile, header=FALSE,
                              stringsAsFactors = FALSE, col.names = c("chromo", "start", "end", "orientation", "clusterID"))
ctcf_cluster_dt$clusterID <- paste0("cluster", ctcf_cluster_dt$clusterID)



# load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3/ctcf4000_family260000/family_cluster_dt.Rdata")
family_cluster_dt <- get(load(file.path(outFolder, "family_cluster_dt.Rdata")))

head(family_cluster_dt)

#*********************************************************************
 # ctcf4000_family260000
#*********************************************************************
# > family_cluster_dt[grepl("family684_4", family_cluster_dt$clusterID),]
# chromo    start      end entrezID symbol                   family   clusterID
# 15094  chr17 46606807 46608272     3211  HOXB1 HOXL subclass homeoboxes family684_4
# 15095  chr17 46620017 46622393     3212  HOXB2 HOXL subclass homeoboxes family684_4
# 15096  chr17 46626232 46667634     3213  HOXB3 HOXL subclass homeoboxes family684_4
# 15097  chr17 46652869 46655743     3214  HOXB4 HOXL subclass homeoboxes family684_4
# 15098  chr17 46668619 46671103     3215  HOXB5 HOXL subclass homeoboxes family684_4
# 15099  chr17 46673099 46682343     3216  HOXB6 HOXL subclass homeoboxes family684_4
# 15100  chr17 46684594 46688383     3217  HOXB7 HOXL subclass homeoboxes family684_4
# 15101  chr17 46689708 46692474     3218  HOXB8 HOXL subclass homeoboxes family684_4
# 15102  chr17 46698518 46703835     3219  HOXB9 HOXL subclass homeoboxes family684_4
# 15103  chr17 46802125 46806111    10481 HOXB13 HOXL subclass homeoboxes family684_4
# 
# > all_obs_dt[grepl("family684_4", all_obs_dt$clusterID),]
# clusterID chromo clusterSize   minPos   maxPos   midPos                   family
# 11317 family684_4  chr17          10 46606807 46806111 46706459 HOXL subclass homeoboxes
# minDist_CTCFclust closest_CTCFclust clusterLength
# 11317           15912.5 CTCF_cluster19602        199305

# > agg_ctcf_cluster_dt[agg_ctcf_cluster_dt$clusterID == "cluster19602", ]
# clusterID chromo clusterSize   minPos   maxPos   midPos
# 10672 cluster19602  chr17           2 46690202 46690891 46690546
# 
# > ctcf_cluster_dt[ctcf_cluster_dt$clusterID == "cluster19602", ]
# chromo    start      end orientation    clusterID
# 24016  chr17 46690202 46690562           > cluster19602
# 24017  chr17 46690641 46690891           < cluster19602

#*********************************************************************
# ctcf10000_family260000
#*********************************************************************

# >  all_obs_dt[grepl("family684_4", all_obs_dt$clusterID),]
# clusterID chromo clusterSize   minPos   maxPos   midPos                   family
# 11317 family684_4  chr17          10 46606807 46806111 46706459 HOXL subclass homeoboxes
# minDist_CTCFclust closest_CTCFclust clusterLength
# 11317           13250.5 CTCF_cluster15910        199305
# # > agg_ctcf_cluster_dt[agg_ctcf_cluster_dt$clusterID == "cluster15910", ]
# clusterID chromo clusterSize   minPos   maxPos   midPos
# 6570 cluster15910  chr17           3 46690202 46696215 46693208
# # >  ctcf_cluster_dt[ctcf_cluster_dt$clusterID == "cluster15910", ]
# chromo    start      end orientation    clusterID
# 24016  chr17 46690202 46690562           > cluster15910
# 24017  chr17 46690641 46690891           < cluster15910
# 24018  chr17 46696009 46696215           < cluster15910

#*********************************************************************
#*********************************************************************
#*********************************************************************

outFolder <- "LOOK_INSIDE_CLUSTERS"
dir.create(outFolder, recursive = TRUE)


all_ctcf_clustSize <- setNames(agg_ctcf_cluster_dt$clusterSize, paste0("CTCF_", agg_ctcf_cluster_dt$clusterID))
stopifnot(all_obs_dt$closest_CTCFclust %in% names(all_ctcf_clustSize))
all_obs_dt$closest_CTCFclustSize <- all_ctcf_clustSize[paste0(all_obs_dt$closest_CTCFclust)]

xqt <- 0.1
yqt <- 0.9

plotCex <- 1.2

myx <- log10(all_obs_dt$minDist_CTCFclust)
myy <- all_obs_dt$closest_CTCFclustSize

xqt_val <- quantile(myx, probs=xqt)
yqt_val <- quantile(myy, probs=yqt)

out_dt <- all_obs_dt[which(myx <= xqt_val & myy >= yqt_val),]


outFile <- file.path(outFolder, 
                     paste0("ctcf", ctcf_clustSize, "_family", family_clustSize, "_xqt", xqt, "_yqt",yqt, "_closestDist_vs_closestSize_withQt.png"))
do.call("png", list(outFile, height=400, width=400))
plot(x=myx,
     y=myy,
     xlab = "minDist_CTCFclust [log10]",
     ylab = "closest_CTCFclustSize",
     pch=16, 
     cex=0.7,
     cex.main=plotCex,
     cex.lab=plotCex,
     cex.axis=plotCex)
abline(h = yqt_val, col="red")
abline(v = xqt_val, col="red")
mtext(side=3, text=paste0("x-qt = ", xqt, " (", round(xqt_val, 2) , "); y-qt = ", yqt, " (",round(yqt_val,2), ")"))
legend("topleft", bty="n", legend=paste0(nrow(out_dt), "/", nrow(all_obs_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

out_dt <- out_dt[order(out_dt$clusterSize, out_dt$closest_CTCFclustSize, decreasing = TRUE),]
out_dt <- out_dt[,c("clusterID", "family", "chromo", "clusterSize", "closest_CTCFclust" , "closest_CTCFclustSize")]

outFile <- file.path(outFolder, 
                     paste0("ctcf", ctcf_clustSize, "_family", family_clustSize, "_xqt", xqt, "_yqt",yqt, "_dt.txt"))
write.table(out_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=FALSE)
cat(paste0("... written: ", outFile, "\n"))
