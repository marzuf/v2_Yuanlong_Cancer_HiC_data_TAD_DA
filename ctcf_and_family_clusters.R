# CTCF_CLUSTER/ctcf_bs_nanni2020.bed             CTCF_CLUSTER/gene2family.bed     
# CTCF_CLUSTER/meanSize.Rdata  CTCF_CLUSTER/test_ctcf_bs_nanni2020.bed

# Rscript ctcf_and_family_clusters.R

require(doMC)
require(foreach)
registerDoMC(40)
require(ggpubr)
require(ggsci)
require(ggplot2)


ggHeight <- 5
ggWidth <- 7
plotTypeGG <- "svg"

chromoSize_dt <- read.delim("CTCF_CLUSTER/hg19_chromo_sizes.txt", stringsAsFactors=FALSE, col.names=c("chromo", "size"), header=F)

nPermut <- 1000

ctcf_clustSize <- 4000
family_clustSize <- 260000


outFolder <- file.path("CTCF_AND_FAMILY_CLUSTERS", paste0("ctcf", ctcf_clustSize, "_family", family_clustSize))
dir.create(outFolder, recursive = TRUE)
                       

###########################
# prep CTCF data
###########################

# #../bedtools cluster -i ctcf_bs_nanni2020.bed -d 40000 > ctcf_bs_nanni2020_cluster4kb.bed


ctcf_clustFile <- file.path(outFolder, paste0("ctcf_bs_nanni2020_cluster", ctcf_clustSize, "bp.bed"))

mycmd <- paste0("./bedtools cluster -i CTCF_CLUSTER/ctcf_bs_nanni2020.bed -d ", ctcf_clustSize, 
" > ", ctcf_clustFile)
cat(paste0("> ", mycmd, "\n"))
system(mycmd)

ctcf_cluster_dt <- read.delim(ctcf_clustFile, header=FALSE,
                              stringsAsFactors = FALSE, col.names = c("chromo", "start", "end", "orientation", "clusterID"))
ctcf_cluster_dt$clusterID <- paste0("cluster", ctcf_cluster_dt$clusterID)

# for each ctcf cluster -> cluster size and midposition
aggSize_ctcf_cluster_dt <- aggregate(orientation~clusterID+chromo, FUN=length, data=ctcf_cluster_dt)
colnames(aggSize_ctcf_cluster_dt)[colnames(aggSize_ctcf_cluster_dt) == "orientation"] <- "clusterSize"

aggMinPos_ctcf_cluster_dt <- aggregate(start~clusterID+chromo, FUN=min, data=ctcf_cluster_dt)
colnames(aggMinPos_ctcf_cluster_dt)[colnames(aggMinPos_ctcf_cluster_dt) == "start"] <- "minPos"

aggMaxPos_ctcf_cluster_dt <- aggregate(end~clusterID+chromo, FUN=max, data=ctcf_cluster_dt)
colnames(aggMaxPos_ctcf_cluster_dt)[colnames(aggMaxPos_ctcf_cluster_dt) == "end"] <- "maxPos"

agg_ctcf_cluster_dt <- merge(merge(aggSize_ctcf_cluster_dt, aggMinPos_ctcf_cluster_dt, by=c("clusterID","chromo"), all=TRUE),
                             aggMaxPos_ctcf_cluster_dt, by=c("clusterID", "chromo"), all=TRUE)
cat(paste0("# clusters CTCF:\t", nrow(agg_ctcf_cluster_dt), "\n"))
cat(paste0("mean size clusters CTCF:\t", mean(agg_ctcf_cluster_dt$clusterSize), "\n"))
agg_ctcf_cluster_dt$midPos <- (agg_ctcf_cluster_dt$minPos+agg_ctcf_cluster_dt$maxPos)/2

outFile <- file.path(outFolder, "agg_ctcf_cluster_dt.Rdata")
save(agg_ctcf_cluster_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

###########################
# prep family data
###########################


family_clustFile <- file.path(outFolder, paste0("gene2family_cluster", family_clustSize, "bp.bed"))

mycmd <- paste0("./bedtools cluster -i CTCF_CLUSTER/gene2family.bed -d ", family_clustSize, 
                " > ", family_clustFile)
cat(paste0("> ", mycmd, "\n"))
system(mycmd)


family_cluster_dt <- read.delim(family_clustFile, header=FALSE,
                                stringsAsFactors = FALSE, 
                                col.names = c("chromo", "start", "end", "entrezID", "symbol", "family", "clusterID"))
family_cluster_dt$clusterID <- paste0("cluster", family_cluster_dt$clusterID)


# for each family, for each gene cluster -> cluster size and midposition
# fam = unique(family_cluster_dt$family)[1]
all_fams_agg_cluster_dt <- foreach(fam = unique(family_cluster_dt$family), .combine='rbind') %dopar% {
  
  subFamCluster_dt <- family_cluster_dt[family_cluster_dt$family == fam,]
  stopifnot(nrow(subFamCluster_dt) > 0)
  
  # for each gene cluster -> cluster size and midposition
  aggSize_subFamCluster_dt <- aggregate(family~clusterID+chromo, FUN=length, data=subFamCluster_dt)
  colnames(aggSize_subFamCluster_dt)[colnames(aggSize_subFamCluster_dt) == "family"] <- "clusterSize"
  
  aggMinPos_subFamCluster_dt <- aggregate(start~clusterID+chromo, FUN=min, data=subFamCluster_dt)
  colnames(aggMinPos_subFamCluster_dt)[colnames(aggMinPos_subFamCluster_dt) == "start"] <- "minPos"
  
  aggMaxPos_subFamCluster_dt <- aggregate(end~clusterID+chromo, FUN=max, data=subFamCluster_dt)
  colnames(aggMaxPos_subFamCluster_dt)[colnames(aggMaxPos_subFamCluster_dt) == "end"] <- "maxPos"
  
  agg_subFamCluster_dt <- merge(merge(aggSize_subFamCluster_dt, aggMinPos_subFamCluster_dt, by=c("clusterID","chromo"), all=TRUE),
                               aggMaxPos_subFamCluster_dt, by=c("clusterID","chromo"), all=TRUE)
  # cat(paste0("# clusters CTCF:\t", nrow(agg_subFamCluster_dt), "\n"))
  # cat(paste0("mean size clusters CTCF:\t", mean(agg_subFamCluster_dt$clusterSize), "\n"))
  agg_subFamCluster_dt$midPos <- (agg_subFamCluster_dt$minPos+agg_subFamCluster_dt$maxPos)/2
  agg_subFamCluster_dt$family <- fam
  agg_subFamCluster_dt
}
outFile <- file.path(outFolder, "all_fams_agg_cluster_dt.Rdata")
save(all_fams_agg_cluster_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

nrow(all_fams_agg_cluster_dt) # 11228
sum(all_fams_agg_cluster_dt$clusterSize > 1) # 1744

stopifnot(!duplicated(file.path(all_fams_agg_cluster_dt$family,all_fams_agg_cluster_dt$clusterID)))
stopifnot(!duplicated(file.path(agg_ctcf_cluster_dt$clusterID)))

final_fams_cluster_dt <- all_fams_agg_cluster_dt[all_fams_agg_cluster_dt$clusterSize > 1,]
final_ctcf_cluster_dt <- agg_ctcf_cluster_dt[agg_ctcf_cluster_dt$clusterSize > 1,]

i=1

final_fams_cluster_dt$minDist_CTCFclust <- NA
final_fams_cluster_dt$closest_CTCFclust <- NA

final_fams_cluster_dt$clusterLength <- final_fams_cluster_dt$maxPos - final_fams_cluster_dt$minPos + 1

obsAndPermut_closest_clusters_dt <- foreach(i = 1:nrow(final_fams_cluster_dt)) %dopar% {
  
  curr_chromo <- final_fams_cluster_dt$chromo[i]
  curr_midpos <- final_fams_cluster_dt$midPos[i]
  
  curr_length <- final_fams_cluster_dt$clusterLength[i]
  
  sub_ctcf_dt <- final_ctcf_cluster_dt[final_ctcf_cluster_dt$chromo == curr_chromo,]
  

  i_minDist <- which.min(abs(curr_midpos - sub_ctcf_dt$midPos))
  
  final_fams_cluster_dt$closest_CTCFclust[i] <- paste0("CTCF_", sub_ctcf_dt$clusterID[i_minDist])
  final_fams_cluster_dt$minDist_CTCFclust[i] <- min(abs(curr_midpos - sub_ctcf_dt$midPos))
  obs_final_fams_cluster_dt <- final_fams_cluster_dt[i,]
  
  maxSize <- chromoSize_dt$size[chromoSize_dt==curr_chromo]
  stopifnot(length(maxSize) == 1)
  stopifnot(is.numeric(maxSize))
  
  permut_minDist <- foreach(i_permut=1:nPermut, .combine='c') %do% {
    
    # for a chromo of size 11, and length of 5
    # last position can be 11-5+1=7 (7-8-9-10-11)
    randomStartPos <- sample(1:(maxSize-curr_length+1), size=1)
    # if randomStartPos is 5; randomEndPos is 7+5-1
    randomEndPos <- randomStartPos + curr_length - 1
    random_midpos <- (randomStartPos+randomEndPos)/2
    random_minDist <- min(abs(random_midpos - sub_ctcf_dt$midPos))
    random_minDist
  }
  list(
    obs_final_fams_cluster_dt=obs_final_fams_cluster_dt,
    permut_minDist = permut_minDist
  )
}


outFile <- file.path(outFolder, "obsAndPermut_closest_clusters_dt.Rdata")
save(obsAndPermut_closest_clusters_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

# load("CTCF_AND_FAMILY_CLUSTERS/obsAndPermut_closest_clusters_dt.Rdata")
# plot # clust per family; mean clust size

all_obs_dt <- do.call('rbind', lapply(obsAndPermut_closest_clusters_dt, function(x) x[["obs_final_fams_cluster_dt"]]))
all_permut_dt <- unlist(lapply(obsAndPermut_closest_clusters_dt, function(x) x[["permut_minDist"]]))

plot_dt <- data.frame(
  dataType=c(rep("observed", nrow(all_obs_dt)), rep("random", length(all_permut_dt))),
  minDist = c(all_obs_dt$minDist_CTCFclust, all_permut_dt),
  stringsAsFactors=FALSE
)
# 
plot_dt$minDist_log10 <- log10(plot_dt$minDist)

plotTit <- "min dist. to CTCF cluster"
subTit <- paste0("# permut = ", nPermut, "; # family clusters = ", nrow(all_obs_dt), " (# fam. = ", length(unique(all_obs_dt$family)), ")") 

p <- ggdensity(plot_dt,
               x = paste0("minDist_log10"),
               y = "..density..",
               # combine = TRUE,                  # Combine the 3 plots
               xlab = paste0("minDist to CTCF cluster [log10]"),
               # add = "median",                  # Add median line.
               rug = FALSE,                      # Add marginal rug
               color = "dataType",
               fill = "dataType",
               palette = "d3"
)+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
  theme(plot.title=element_text(face="bold"),
        plot.subtitle=element_text(face="italic")
        )
outFile <- file.path(outFolder, paste0("minDist_toCTCFcluster_obsPermut_densityplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))


plotTit <- "Family cluster size"
subTit <- paste0("# family clusters = ", nrow(all_obs_dt), " (# fam. = ", length(unique(all_obs_dt$family)), ")") 

p <- ggdensity(all_fams_agg_cluster_dt,
               x = paste0("clusterSize"),
               y = "..density..",
               # combine = TRUE,                  # Combine the 3 plots
               xlab = paste0("minDist to CTCF cluster [log10]"),
               # add = "median",                  # Add median line.
               rug = FALSE,                      # Add marginal rug
               # color = "dataType",
               # fill = "dataType",
               palette = "d3"
)+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
  theme(plot.title=element_text(face="bold"),
        plot.subtitle=element_text(face="italic")
  )
outFile <- file.path(outFolder, paste0("minDist_toCTCFcluster_obsPermut_densityplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))

p <- ggdensity(all_fams_agg_cluster_dt[all_fams_agg_cluster_dt$clusterSize>1,],
               x = paste0("clusterSize"),
               y = "..density..",
               # combine = TRUE,                  # Combine the 3 plots
               xlab = paste0("minDist to CTCF cluster [log10]"),
               # add = "median",                  # Add median line.
               rug = FALSE,                      # Add marginal rug
               # color = "dataType",
               # fill = "dataType",
               palette = "d3"
)+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
  theme(plot.title=element_text(face="bold"),
        plot.subtitle=element_text(face="italic")
  )
outFile <- file.path(outFolder, paste0("minDist_toCTCFcluster_obsPermut_densityplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))


# for each family cluster -> min dist to next ctcf cluster
# random genomic stretch of same size -> min dist to next ctcf cluster

#   CTCF_CLUSTER/gene2family_cluster260kb.bed  CTCF_CLUSTER/prep_data.R
# 
# chr1    14362   29370   653635  WASH7P  Wiskott-Aldrich Syndrome protein family 1
# chr1    30366   30503   100302278       MIR1302-2       MicroRNAs       1
# chr1    52453   53396   79504   OR4G4P  Olfactory receptors, family 4   1
# chr1    63016   63885   403263  OR4G11P Olfactory receptors, family 4   1
# chr1    69091   70008   79501   OR4F5   Olfactory receptors, family 4   1
# chr1    367659  368597  729759  OR4F29  Olfactory receptors, family 4   2
# chr1    621096  622034  81399   OR4F16  Olfactory receptors, family 4   2
# chr1    859993  879961  148398  SAMD11  Sterile alpha motif domain containing   2
# chr1    879583  894679  26155   NOC2L   Protein phosphatase 1 regulatory subunits       2
# 
# 
# chr1    237593  237953  >       1
# chr1    521337  521697  >       2
# chr1    714087  714447  >       3
# chr1    805232  805362  >       4
# chr1    839966  840326  >       4
# chr1    848092  848452  >       4
# chr1    856454  856704  <       4
# chr1    864292  864508  <       4
# chr1    873488  873848  <       4
# chr1    880699  880915  <       4


# for each family cluster -> should filter same family inside cluster !
# filter more than two

# should do something with random genes
# mid pos of the cluster

# smaller dist to midpos of ctcf clusters
