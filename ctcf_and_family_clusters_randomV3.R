# CTCF_CLUSTER/ctcf_bs_nanni2020.bed             CTCF_CLUSTER/gene2family.bed     
# CTCF_CLUSTER/meanSize.Rdata  CTCF_CLUSTER/test_ctcf_bs_nanni2020.bed

# Rscript ctcf_and_family_clusters_randomV3.R


require(doMC)
require(foreach)
registerDoMC(40)
require(ggpubr)
require(ggsci)
require(ggplot2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

ggHeight <- 5
ggWidth <- 7
plotTypeGG <- "svg"

chromoSize_dt <- read.delim("CTCF_CLUSTER/hg19_chromo_sizes.txt", stringsAsFactors=FALSE, col.names=c("chromo", "size"), header=F)

nPermut <- 100

# ctcf_clustSize <- 4000
ctcf_clustSize <- 10000
# ctcf_clustSize <- 20000
family_clustSize <- 260000
# family_clustSize <- 130000


outFolder <- file.path("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3", paste0("ctcf", ctcf_clustSize, "_family", family_clustSize))
dir.create(outFolder, recursive = TRUE)
                       
tmpFolder <- file.path(outFolder, "TMP")
dir.create(tmpFolder, recursive = TRUE)

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
stopifnot(!duplicated(aggSize_ctcf_cluster_dt$clusterID))
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
# load("CTCF_AND_FAMILY_CLUSTERS/ctcf4000_family260000/agg_ctcf_cluster_dt.Rdata")


###########################
# prep family data
###########################

all_family_bed <- read.delim("CTCF_CLUSTER/gene2family.bed", header=FALSE,
                             stringsAsFactors = FALSE, 
                             col.names = c("chromo", "start", "end", "entrezID", "symbol", "family"))
allfams <- unique(all_family_bed$family)

fam=allfams[28]
i_fam=185

family_cluster_dt <- foreach(i_fam=1:length(allfams), .combine='rbind') %dopar% {
  
  fam <- allfams[i_fam]
  
  sub_family_bed <- all_family_bed[all_family_bed$family == fam,]
  stopifnot(nrow(sub_family_bed) > 0)
  
  sub_family_file <- file.path(tmpFolder, paste0(gsub(" ", "_", gsub("/", "_", gsub("\\(", "", gsub("\\)", "", gsub("'", "", fam))))),
                                                 "_gene2family.bed"))
  write.table(sub_family_bed, file=sub_family_file, col.names=F, row.names=F, sep="\t")
  cat(paste0("... written: ", sub_family_file, "\n"))
  
  mycmd <- paste0("./bedtools cluster -i ", sub_family_file, " -d ", family_clustSize)
                  # " > ", family_clustFile)
  cat(paste0("> ", mycmd, "\n"))
  sub_family_cluster_dt <-  read.table(textConnection(system(mycmd, intern=TRUE)), header=FALSE,
                                       col.names = c("chromo", "start", "end", "entrezID", "symbol", "family", "clusterID"))
  sub_family_cluster_dt$clusterID <- paste0("family", i_fam, "_",sub_family_cluster_dt$clusterID )
  sub_family_cluster_dt
}

outFile <- file.path(outFolder, "family_cluster_dt.Rdata")
#load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3/ctcf4000_family260000/family_cluster_dt.Rdata")
save(family_cluster_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

xx <- table(family_cluster_dt$clusterID)
length(xx)
# [1] 12465
nrow(family_cluster_dt)
# [1] 17244
xx2 <- xx[xx > 1]
length(xx2)
# 1380
# range(xx2)
# [1]   2 171

# stop("-ok\n")

# WRONG VERSION
# family_clustFile <- file.path(outFolder, paste0("gene2family_cluster", family_clustSize, "bp.bed"))
# 
# mycmd <- paste0("./bedtools cluster -i CTCF_CLUSTER/gene2family.bed -d ", family_clustSize, 
#                 " > ", family_clustFile)
# cat(paste0("> ", mycmd, "\n"))
# system(mycmd)
# 
# 
# family_cluster_dt <- read.delim(family_clustFile, header=FALSE,
#                                 stringsAsFactors = FALSE, 
#                                 col.names = c("chromo", "start", "end", "entrezID", "symbol", "family", "clusterID"))


###########################
# prep all genes data
###########################
setDir="/media/electron"
setDir=""
allgenes_dt <- read.delim(file.path(setDir,
                                    "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                          header=T, stringsAsFactors = FALSE)
allgenes_dt$entrezID <- as.character(allgenes_dt$entrezID)
allgenes_dt$symbol <- as.character(allgenes_dt$symbol)
stopifnot(!duplicated(allgenes_dt$entrezID))
allgenes_dt <- allgenes_dt[,c("chromo", "start", "end", "entrezID")]
allgenes_dt <- allgenes_dt[order(allgenes_dt$chromo, allgenes_dt$start, allgenes_dt$end),]

tmp_gene_file <- file.path(tmpFolder, "all_genes.bed")
write.table(allgenes_dt, file=tmp_gene_file, col.names=F, row.names=F, quote=F, sep="\t")

allgenes_clustFile <- file.path(tmpFolder, paste0("all_genes_clust", family_clustSize, ".bed"))

mycmd <- paste0("./bedtools cluster -i ", tmp_gene_file, " -d ", family_clustSize, 
                " > ", allgenes_clustFile)
cat(paste0("> ", mycmd, "\n"))
system(mycmd)

allgenes_cluster_dt <- read.delim(allgenes_clustFile, header=FALSE,
                                  stringsAsFactors = FALSE, col.names = c("chromo", "start", "end", "entrezID", "clusterID"))
allgenes_cluster_dt$clusterID <- paste0("allgenes", allgenes_cluster_dt$clusterID)
# > allgenes_clustFile="CTCF_AND_FAMILY_CLUSTERS_RANDOMV2/ctcf4000_family260000/TMP/all_genes_clust260000.bed"
allgenes_cluster_dt$entrezID <- as.character(allgenes_cluster_dt$entrezID)
all_family_bed$entrezID <- as.character(all_family_bed$entrezID)
withFam_allgenes_cluster_dt <- merge(allgenes_cluster_dt,all_family_bed, by=c("entrezID", "chromo", "start", "end"),
                                     all.x=TRUE, all.y=FALSE)
save(withFam_allgenes_cluster_dt, file="withFam_allgenes_cluster_dt.Rdata", version=2)
stopifnot(!all(is.na(withFam_allgenes_cluster_dt$family)))
stopifnot(!duplicated(withFam_allgenes_cluster_dt$entrezID))

withFam_allgenes_cluster_dt$family[is.na(withFam_allgenes_cluster_dt$family)] <- "unknown"

agg_allGenesFamClust_dt <- aggregate(entrezID~clusterID, data = withFam_allgenes_cluster_dt, FUN=length)
colnames(agg_allGenesFamClust_dt)[colnames(agg_allGenesFamClust_dt) == "entrezID"] <- "clusterSize"

# retrieve the max ratio of a same family
maxSameFam_allGenesFamClust_dt <- aggregate(family~clusterID, data= withFam_allgenes_cluster_dt,
                                         FUN=function(x){
                                           ratioFams <- table(x)/length(x)
                                           stopifnot(abs(sum(ratioFams)-1) < 10^-5)
                                           max(ratioFams[names(ratioFams) != "unknown"]) # if only unknown will return -Inf
                                         })

sameFamThresh <- 0.5

toKeep <- maxSameFam_allGenesFamClust_dt$clusterID[maxSameFam_allGenesFamClust_dt$family <= sameFamThresh]
# > length(toKeep)
# [1] 1227
# > nrow(maxSameFam_allGenesFamClust_dt)
# [1] 1268
# > 
filtered_agg_allGenesFamClust_dt <- agg_allGenesFamClust_dt[agg_allGenesFamClust_dt$clusterID %in% toKeep,]
stopifnot(nrow(filtered_agg_allGenesFamClust_dt) == length(toKeep))

all_genes_minPos <- aggregate(start~clusterID , data=allgenes_cluster_dt, FUN=min)
colnames(all_genes_minPos)[colnames(all_genes_minPos) == "start"] <- "minPos"
all_genes_maxPos <- aggregate(end~clusterID , data=allgenes_cluster_dt, FUN=max)
colnames(all_genes_maxPos)[colnames(all_genes_maxPos) == "end"] <- "maxPos"

all_genes_clust_dt <- merge(filtered_agg_allGenesFamClust_dt, 
                            merge(all_genes_maxPos,all_genes_minPos, by=c("clusterID"),all=FALSE),
                                                                     by=c("clusterID"),all=FALSE)
stopifnot(!is.na(all_genes_clust_dt))
stopifnot(nrow(all_genes_clust_dt) == nrow(filtered_agg_allGenesFamClust_dt))

# stop("-ok\n")


# > rd=unique(agg_allGenesFamClust_dt$entrezID)
# [1]  142   65    1    3    9   18   35    2    8  213    4   75   19   14   15    7  115  147   82  134
# [21]   10    6   23   32  262  211   29    5   26   25  148  237  131   34  102  111   47   54   53   37
# [41]   12   48   42   31   67   27   44   40  149   20   49   17   21  129  174  120   16   95   13   22
# [61]  251   11   93  101  343   28   45   24   30  292  412  140   85  112   41   38   74  187   52   64
# [81]  177   55  117  620  184   98   60   50   33  375  429   78  542   51   77   56  130  159   46   92
# [101]  138  106   36   39   76  401   72  207  216  103  186  224   71  107  153  166   68   73  137  750
# [121]  104  295   63  121  291   91  192  179  351  176  497  190   58  271   81   90  915 1196  168   61
# [141]   83   80  154  135   99  194   43  110  150  214  128  125   62  339   79  124  444  290  272   87
# [161]  230   70  304  161  109  182  773
# > obs=all_fams_agg_cluster_dt$clusterSize
# > obs=unique(obs)
# > obs
# [1]   1   2   3   4   6  53   7   5  51   8  15  21  11  22  25   9  12  27  23  30  38  16  10  18  97
# [26]  47  33  20  17  41  81  13  43  14  37  26  44  57 117 171  36  88  63  19  60
# > all(obs %in% rd)
# [1] FALSE
# > obs[which(!obs %in% rd)]
# [1]  97  57 171  88




############################################################ RUN THE MATCH
commonChromo <- intersect(family_cluster_dt$chromo, agg_ctcf_cluster_dt$chromo)
stopifnot(commonChromo %in% paste0("chr", c(1:22, "X")))

family_cluster_dt <- family_cluster_dt[family_cluster_dt$chromo %in% commonChromo,]
stopifnot(nrow(family_cluster_dt) > 0)

agg_ctcf_cluster_dt <- agg_ctcf_cluster_dt[agg_ctcf_cluster_dt$chromo %in% commonChromo,]
stopifnot(nrow(agg_ctcf_cluster_dt) > 0)

# for each family, for each gene cluster -> cluster size and midposition
# fam = unique(family_cluster_dt$family)[1]
all_fams_agg_cluster_dt <- foreach(fam = unique(family_cluster_dt$family), .combine='rbind') %dopar% {
  
  subFamCluster_dt <- family_cluster_dt[family_cluster_dt$family == fam,]
  stopifnot(nrow(subFamCluster_dt) > 0)
  
  # for each gene cluster -> cluster size, minpos, maxpos and midposition
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
#  load("CTCF_AND_FAMILY_CLUSTERS/ctcf4000_family260000/all_fams_agg_cluster_dt.Rdata")

nrow(all_fams_agg_cluster_dt) # 11228
sum(all_fams_agg_cluster_dt$clusterSize > 1) # 1744

stopifnot(!duplicated(file.path(all_fams_agg_cluster_dt$family,all_fams_agg_cluster_dt$clusterID)))
stopifnot(!duplicated(file.path(agg_ctcf_cluster_dt$clusterID)))

#### take only "real clusters" i.e. > 1

final_fams_cluster_dt <- all_fams_agg_cluster_dt[all_fams_agg_cluster_dt$clusterSize > 1,]
final_ctcf_cluster_dt <- agg_ctcf_cluster_dt[agg_ctcf_cluster_dt$clusterSize > 1,]
final_ctcf_cluster_dt$clusterID_full <- paste0("CTCF_", final_ctcf_cluster_dt$clusterID)
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

                # final_fams_cluster_dt$closest_CTCFclust[i] <- paste0("CTCF_", sub_ctcf_dt$clusterID[i_minDist])
                final_fams_cluster_dt$closest_CTCFclust[i] <- sub_ctcf_dt$clusterID_full[i_minDist]
                final_fams_cluster_dt$minDist_CTCFclust[i] <- min(abs(curr_midpos - sub_ctcf_dt$midPos))
                obs_final_fams_cluster_dt <- final_fams_cluster_dt[i,]
                
                stopifnot(  final_fams_cluster_dt$closest_CTCFclust[i]  %in% final_ctcf_cluster_dt$clusterID_full)
                
                
                # find closest in size all genes cluster
                # v3 => sample among the min values instead of taking the first matching (which.min)
                curr_size <- final_fams_cluster_dt$clusterSize[i]
                
                i_sameSizeV0 <- which.min(curr_size-all_genes_clust_dt$clusterSize)
                allDiff <- curr_size-all_genes_clust_dt$clusterSize
                whichSmallestDiff <- which(allDiff == min(allDiff))
                stopifnot(i_sameSizeV0 %in% whichSmallestDiff)
                i_sameSize <- sample(whichSmallestDiff, size=1)
                
                
                  randomStartPos <- all_genes_clust_dt$minPos[i_sameSize]
                  randomEndPos <- all_genes_clust_dt$maxPos[i_sameSize]
                  random_midpos <- (randomStartPos+randomEndPos)/2
                  
                  i_random_minDist <- which.min(abs(random_midpos - sub_ctcf_dt$midPos))
                  random_closest  <- sub_ctcf_dt$clusterID_full[i_random_minDist]
                  
                  stopifnot(  random_closest  %in% final_ctcf_cluster_dt$clusterID_full)
                  
                  
                  random_minDist <- min(abs(random_midpos - sub_ctcf_dt$midPos))
  
                  
                list(
                  obs_final_fams_cluster_dt=obs_final_fams_cluster_dt,
                  random_closest=random_closest,
                  random_minDist=random_minDist
                )
              }


              outFile <- file.path(outFolder, "obsAndPermut_closest_clusters_dt.Rdata")
              save(obsAndPermut_closest_clusters_dt, file=outFile, version=2)
              cat(paste0("... written: ", outFile, "\n"))

# load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3/ctcf4000_family260000/obsAndPermut_closest_clusters_dt.Rdata")


all_obs_dt <- do.call('rbind', lapply(obsAndPermut_closest_clusters_dt, function(x) x[["obs_final_fams_cluster_dt"]]))
all_permutDist <- unlist(lapply(obsAndPermut_closest_clusters_dt, function(x) x[["random_minDist"]]))
all_permutDist <- as.numeric(all_permutDist)
stopifnot(!is.na(all_permutDist))
all_permutClosest <- unlist(lapply(obsAndPermut_closest_clusters_dt, function(x) x[["random_closest"]]))
all_permutClosest <- as.character(all_permutClosest)
stopifnot(!is.na(all_permutClosest))
stopifnot(length(all_permutDist) == length(all_permutClosest))

stopifnot(all_obs_dt$closest_CTCFclust %in% final_ctcf_cluster_dt$clusterID_full)

tmp_fam_dt <- all_obs_dt[,c("clusterID", "clusterSize", "closest_CTCFclust")]
colnames(tmp_fam_dt)<- c("family_clusterID", "family_clusterSize", "CTCF_clusterID")
tmp_ctcf_dt <- agg_ctcf_cluster_dt[,c("clusterID", "clusterSize")]
colnames(tmp_ctcf_dt) <- c("CTCF_clusterID", "CTCF_clusterSize")
tmp_ctcf_dt$CTCF_clusterID <- paste0("CTCF_", tmp_ctcf_dt$CTCF_clusterID)
tmp_merged_dt <- merge(tmp_fam_dt, tmp_ctcf_dt, all.x=TRUE, all.y=FALSE, by="CTCF_clusterID")
stopifnot(!is.na(tmp_merged_dt))

plotCex <- 1.2
outFile <- file.path(outFolder, paste0("famClusterSize_vs_CTCFclusterSize_densplot.png"))
do.call("png", list(outFile, height=400, width=400))
densplot(x=tmp_merged_dt$CTCF_clusterSize,
         y=tmp_merged_dt$family_clusterSize,
         xlab="CTCF_clusterSize",
         ylab="family_clusterSize",
         cex.main=plotCex,
         cex.lab=plotCex,
         cex.axis=plotCex)
mtext(side=3, text = paste0("n = ", nrow(tmp_merged_dt)))
addCorr(x=tmp_merged_dt$CTCF_clusterSize,
        y=tmp_merged_dt$family_clusterSize,bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))

# stop("-ok\n")

################################################################################################
######################### clust size and dist of the of the closest ctcf, + filter max dist
################################################################################################


agg_ctcf_cluster_dt$CTCFid <- paste0("CTCF_", agg_ctcf_cluster_dt$clusterID)

random_match_dt <- data.frame(
  CTCFid = all_permutClosest,
  minDist_CTCFclust = as.numeric(all_permutDist),
  stringsAsFactors=FALSE
)
stopifnot(!is.na(random_match_dt))
mean(random_match_dt$minDist_CTCFclust)
stopifnot(random_match_dt$CTCFid %in% final_ctcf_cluster_dt$clusterID_full)

merged_random_dt <- merge(random_match_dt, agg_ctcf_cluster_dt, by="CTCFid", all.x=TRUE, all.y=FALSE, suffixes = c("", "_CTCF"))
stopifnot(!is.na(merged_random_dt))
mean(merged_random_dt$clusterSize)

mean(all_obs_dt$clusterSize)
mean(all_obs_dt$minDist_CTCFclust)
all_obs_dt$CTCFid <- all_obs_dt$closest_CTCFclust
merged_obs_dt <- merge(all_obs_dt, agg_ctcf_cluster_dt, by="CTCFid", all.x=TRUE, all.y=FALSE, suffixes = c("", "_CTCF"))
stopifnot(!is.na(merged_obs_dt))


plot_dt <- data.frame(
  dataType=c(rep("observed", nrow(merged_obs_dt)), rep("random", nrow(merged_random_dt))),
  closestSize = c(merged_obs_dt$clusterSize_CTCF, merged_random_dt$clusterSize),
  closestDist = c(merged_obs_dt$minDist_CTCFclust, merged_random_dt$minDist_CTCFclust),
  stringsAsFactors = FALSE
)
plot_dt$closestDist_log10 <- log10(plot_dt$closestDist)
  
mean(plot_dt$closestSize[plot_dt$dataType=="observed"])
mean(merged_obs_dt$clusterSize_CTCF)

mean(plot_dt$closestSize[plot_dt$dataType=="random"])
mean(merged_random_dt$clusterSize)

plot_vars <- c("closestSize", "closestDist", "closestDist_log10")

distThresh <- 4000

for(plot_var in plot_vars) {
  
  plotTit <- paste0(plot_var, " - closest CTCF cluster")
  subTit <- paste0("# ", names(table(plot_dt$dataType)), "=", as.numeric(table(plot_dt$dataType)), collapse="; ") 
    #"family cluster's closest CTCF cluster size" #paste0("# permut = ", nPermut, "; # family clusters = ", nrow(all_obs_dt), " (# fam. = ", length(unique(all_obs_dt$family)), ")") 
  
  p <- ggdensity(plot_dt,
                 x = paste0(plot_var),
                 y = "..density..",
                 # combine = TRUE,                  # Combine the 3 plots
                 xlab = paste0(plot_var, " - closest CTCF cluster"),
                 # add = "median",                  # Add median line.
                 rug = FALSE,                      # Add marginal rug
                 color = "dataType",
                 fill = "dataType",
                 palette = "d3"
  )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
    theme(plot.title=element_text(face="bold"),
          plot.subtitle=element_text(face="italic")
    )
  outFile <- file.path(outFolder, paste0(plot_var, "_closestCTCF_obsPermut_densityplot.", plotTypeGG))
  ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  
  p <- ggboxplot(plot_dt,
                 x = paste0("dataType"),
                 y = paste0(plot_var),
                 # combine = TRUE,                  # Combine the 3 plots
                 ylab = paste0(plot_var, " - closest CTCF cluster"),
                 xlab = "",
                 add = "jitter",                  # Add median line.
                 rug = FALSE,                      # Add marginal rug
                 color = "dataType",
                 # fill = "dataType",
                 palette = "d3"
  )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
    theme(plot.title=element_text(face="bold"),
          plot.subtitle=element_text(face="italic")
    )
  outFile <- file.path(outFolder, paste0(plot_var, "_closestCTCF_obsPermut_boxplot.", plotTypeGG))
  ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  
  
  
  
  plotTit <- paste0(plot_var, " - closest CTCF  (max dist. ", distThresh, ")")
  subTit <- paste0("# ", names(table(plot_dt$dataType)), "=", as.numeric(table(plot_dt$dataType)), collapse="; ") 
  #"family cluster's closest CTCF cluster size" #paste0("# permut = ", nPermut, "; # family clusters = ", nrow(all_obs_dt), " (# fam. = ", length(unique(all_obs_dt$family)), ")") 
  
  p <- ggdensity(plot_dt[plot_dt$closestDist <= distThresh,],
                 x = paste0(plot_var),
                 y = "..density..",
                 # combine = TRUE,                  # Combine the 3 plots
                 xlab = paste0(plot_var, " - closest CTCF cluster"),
                 # add = "median",                  # Add median line.
                 rug = FALSE,                      # Add marginal rug
                 color = "dataType",
                 # fill = "dataType",
                 palette = "d3"
  )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
    theme(plot.title=element_text(face="bold"),
          plot.subtitle=element_text(face="italic")
    )
  outFile <- file.path(outFolder, paste0(plot_var, "_closestCTCF_obsPermut_maxDist", distThresh, "_densityplot.", plotTypeGG))
  ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  p <- ggboxplot(plot_dt[plot_dt$closestDist <= distThresh,],
                 y = paste0(plot_var),
                 x = "dataType",
                 # combine = TRUE,                  # Combine the 3 plots
                 xlab ="",
                 ylab = paste0(plot_var, " - closest CTCF cluster"),
                 add = "jitter",                  # Add median line.
                 rug = FALSE,                      # Add marginal rug
                 color = "dataType",
                 # fill = "dataType",
                 palette = "d3"
  )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
    theme(plot.title=element_text(face="bold"),
          plot.subtitle=element_text(face="italic")
    )
  outFile <- file.path(outFolder, paste0(plot_var, "_closestCTCF_obsPermut_maxDist", distThresh, "_boxplot.", plotTypeGG))
  ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
  cat(paste0("... written: ", outFile,  "\n"))
  
  
}

source("ctcf_da_utils.R")
sizeBreaks <- c(2,4,6)
# cat(paste0("... running get_fract_lab2\n"))

save(plot_dt, file="plot_dt.Rdata", version=2)

barplot_dt <- plot_dt
barplot_dt$closestSize_labs <- get_fract_lab2(vect_values=barplot_dt$closestSize, range_levels = sizeBreaks)

agg_barplot_dt <- aggregate(closestDist~dataType+closestSize, data=barplot_dt, FUN=length)


agg_barplot_dt <- aggregate(closestDist~dataType+closestSize_labs, data=barplot_dt, FUN=length)
colnames(agg_barplot_dt)[colnames(agg_barplot_dt)=="closestDist"] <- "totCount"
range(agg_barplot_dt$closestSize)
# 2-7
totType <- table(plot_dt$dataType)


agg_barplot_dt$totRatio <- agg_barplot_dt$totCount/as.numeric(totType[agg_barplot_dt$dataType])

plotTit <- "Ratio domains by closest CTCF cluster size"
subTit <- paste0("#", names(totType), "=", totType, collapse="; ")

p <- ggbarplot(agg_barplot_dt, x="dataType", y="totRatio", fill="closestSize_labs", 
               xlab = "", ylab = "Ratio of domains")+
  scale_fill_nejm() + 
  labs(fill="size closest CTCF cluster") + 
  ggtitle(plotTit, subtitle=subTit)+
  theme(
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  )


outFile <- file.path(outFolder, paste0("closestCTCFclustSize_byObsPermut_ratio_barplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))

#########################
######################### desc family and ctcf clusters
#########################

agg_ctcf_cluster_dt$clusterLength <- agg_ctcf_cluster_dt$maxPos-agg_ctcf_cluster_dt$midPos+1
agg_ctcf_cluster_dt$clusterLength_log10 <- log10(agg_ctcf_cluster_dt$clusterLength)
agg_ctcf_cluster_dt$clusterSize_log10 <- log10(agg_ctcf_cluster_dt$clusterSize)

all_fams_agg_cluster_dt$clusterLength <- all_fams_agg_cluster_dt$maxPos-all_fams_agg_cluster_dt$midPos+1
all_fams_agg_cluster_dt$clusterLength_log10 <- log10(all_fams_agg_cluster_dt$clusterLength)
all_fams_agg_cluster_dt$clusterSize_log10 <- log10(all_fams_agg_cluster_dt$clusterSize)


plot_vars <- c("clusterSize", "clusterSize_log10", "clusterLength", "clusterLength_log10")

for(feature in c("Family", "CTCF")) {
  
  if(feature == "Family") {
    plot_dt <- all_fams_agg_cluster_dt
    subTit <- paste0("# ", feature, " clusters = ", length(unique(file.path(plot_dt$clusterID, plot_dt$family))),
                     " (# fam. = ", length(unique(plot_dt$family)), ")") 
    plot_big1_dt <- plot_dt[plot_dt$clusterSize > 1,]
    stopifnot(nrow(plot_big1_dt) > 0)
    subTit2 <- paste0("# ", feature, " clusters = ", length(unique(file.path(plot_big1_dt$clusterID, plot_big1_dt$family))),
                      " (# fam. = ", length(unique(plot_big1_dt$family)), ")") 
    
    
  } else if(feature=="CTCF"){
    plot_dt <- agg_ctcf_cluster_dt
    subTit <- paste0("# ", feature, " clusters = ", length(unique(file.path(plot_dt$clusterID))))
    plot_big1_dt <- plot_dt[plot_dt$clusterSize > 1,]
    stopifnot(nrow(plot_big1_dt) > 0)
    subTit2 <- paste0("# ", feature, " clusters = ", length(unique(file.path(plot_big1_dt$clusterID))))
    
  }else{
    stop("invalid\n")
  }
  
  
  for(plot_var in plot_vars){
    
    plotTit <- paste0(feature, " cluster size")
    
    
    p <- ggdensity(plot_dt,
                   x = paste0(plot_var),
                   y = "..density..",
                   # combine = TRUE,                  # Combine the 3 plots
                   xlab = paste0(feature, " ", plot_var),
                   # add = "median",                  # Add median line.
                   rug = FALSE,                      # Add marginal rug
                   # color = "dataType",
                   # fill = "dataType",
                   palette = "d3"
    )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
      theme(plot.title=element_text(face="bold"),
            plot.subtitle=element_text(face="italic")
      )
    outFile <- file.path(outFolder, paste0(feature, "_", plot_var, "_densityplot.", plotTypeGG))
    ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
    cat(paste0("... written: ", outFile,  "\n"))
    
    plotTit <- paste0(feature, " cluster size")
    
    p <- ggdensity(plot_big1_dt,
                   x = paste0(plot_var),
                   y = "..density..",
                   # combine = TRUE,                  # Combine the 3 plots
                   xlab = paste0(feature, " ", plot_var),
                   # add = "median",                  # Add median line.
                   rug = FALSE,                      # Add marginal rug
                   # color = "dataType",
                   # fill = "dataType",
                   palette = "d3"
    )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit2)+
      theme(plot.title=element_text(face="bold"),
            plot.subtitle=element_text(face="italic")
      )
    outFile <- file.path(outFolder, paste0(feature, "_", plot_var, "_bigger1_densityplot.", plotTypeGG))
    ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
    cat(paste0("... written: ", outFile,  "\n"))
  }
  
  
  
  
}




# 
# plotTit <- "size of closest CTCF cluster"
# subTit <- "family cluster's closest CTCF cluster size" #paste0("# permut = ", nPermut, "; # family clusters = ", nrow(all_obs_dt), " (# fam. = ", length(unique(all_obs_dt$family)), ")") 
# 
# p <- ggdensity(plot_dt,
#                x = paste0("closestSize"),
#                y = "..density..",
#                # combine = TRUE,                  # Combine the 3 plots
#                xlab = paste0("closest CTCF cluster size"),
#                # add = "median",                  # Add median line.
#                rug = FALSE,                      # Add marginal rug
#                color = "dataType",
#                fill = "dataType",
#                palette = "d3"
# )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
#   theme(plot.title=element_text(face="bold"),
#         plot.subtitle=element_text(face="italic")
#   )
# outFile <- file.path(outFolder, paste0("sizeClosestCTCF_byFamilyClust_obsPermut_densityplot.", plotTypeGG))
# ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
# cat(paste0("... written: ", outFile,  "\n"))
# 
# 
# 
# distThresh <- 2000
# 
# plot_dt <- plot_dt[plot_dt$closestDist <= distThresh,]
# plotTit <- paste0("size of closest CTCF cluster (dist. <= ", distThresh, ")")
# subTit <- "family cluster's closest CTCF cluster size" #paste0("# permut = ", nPermut, "; # family clusters = ", nrow(all_obs_dt), " (# fam. = ", length(unique(all_obs_dt$family)), ")") 
# 
# p <- ggdensity(plot_dt,
#                x = paste0("closestSize"),
#                y = "..density..",
#                # combine = TRUE,                  # Combine the 3 plots
#                xlab = paste0("closest CTCF cluster size"),
#                # add = "median",                  # Add median line.
#                rug = FALSE,                      # Add marginal rug
#                color = "dataType",
#                fill = "dataType",
#                palette = "d3"
# )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
#   theme(plot.title=element_text(face="bold"),
#         plot.subtitle=element_text(face="italic")
#   )
# outFile <- file.path(outFolder, paste0("sizeClosestCTCF_byFamilyClust_obsPermut_densityplot.", plotTypeGG))
# ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
# cat(paste0("... written: ", outFile,  "\n"))
# 
# # plot # clust per family; mean clust size
# 
# plot_dt <- data.frame(
#   dataType=c(rep("observed", nrow(all_obs_dt)), rep("random", length(all_permutDist))),
#   minDist = c(all_obs_dt$minDist_CTCFclust, all_permutDist),
#   stringsAsFactors=FALSE
# )
# # 
# plot_dt$minDist_log10 <- log10(plot_dt$minDist)
# 
# plotTit <- "min dist. to CTCF cluster"
# subTit <- paste0("# permut = ", nPermut, "; # family clusters = ", nrow(all_obs_dt), " (# fam. = ", length(unique(all_obs_dt$family)), ")") 
# 
# p <- ggdensity(plot_dt,
#                x = paste0("minDist_log10"),
#                y = "..density..",
#                # combine = TRUE,                  # Combine the 3 plots
#                xlab = paste0("minDist to CTCF cluster [log10]"),
#                # add = "median",                  # Add median line.
#                rug = FALSE,                      # Add marginal rug
#                color = "dataType",
#                fill = "dataType",
#                palette = "d3"
# )+ labs(color="", fill="") + ggtitle(plotTit, subtitle=subTit)+
#   theme(plot.title=element_text(face="bold"),
#         plot.subtitle=element_text(face="italic")
#   )
# outFile <- file.path(outFolder, paste0("minDist_toCTCFcluster_obsPermut_densityplot.", plotTypeGG))
# ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
# cat(paste0("... written: ", outFile,  "\n"))
# 



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
