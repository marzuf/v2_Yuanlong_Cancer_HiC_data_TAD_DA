
# Rscript ctcf_and_family_clusters_v3_v2.R

# take gene family -> take all ctcf +/- 1000 kb
# obs vs random:
# - tot number of ctcf
# - mean cluster size of ctcf

# CTCF_CLUSTER/ctcf_bs_nanni2020.bed             CTCF_CLUSTER/gene2family.bed     
# CTCF_CLUSTER/meanSize.Rdata  CTCF_CLUSTER/test_ctcf_bs_nanni2020.bed


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

aroundFamWindow <- 2000


outFolder <- file.path("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3_V2", paste0("ctcf", ctcf_clustSize, "_family", family_clustSize))
dir.create(outFolder, recursive = TRUE)

tmpFolder <- file.path(outFolder, "TMP")
dir.create(tmpFolder, recursive = TRUE)

all_ctcf_dt <- read.delim("CTCF_CLUSTER/ctcf_bs_nanni2020.bed", stringsAsFactors = FALSE,
                          col.names=c("chromo", "start", "end", "orientation"), header=FALSE)

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

final_fams_cluster_dt$meanSize_matchedCTCFclust <- NA
final_fams_cluster_dt$nbr_matchedCTCFclust <- NA
final_fams_cluster_dt$nbr_matchedCTCF <- NA

final_fams_cluster_dt$clusterLength <- final_fams_cluster_dt$maxPos - final_fams_cluster_dt$minPos + 1

obsAndPermut_closest_clusters_dt <- foreach(i = 1:nrow(final_fams_cluster_dt)) %dopar% {
  
  curr_chromo <- final_fams_cluster_dt$chromo[i]
  # curr_midpos <- final_fams_cluster_dt$midPos[i]
  curr_minpos <- final_fams_cluster_dt$minPos[i]
  curr_maxpos <- final_fams_cluster_dt$maxPos[i]
  
  curr_length <- final_fams_cluster_dt$clusterLength[i]
  
  sub_ctcf_dt <- final_ctcf_cluster_dt[final_ctcf_cluster_dt$chromo == curr_chromo,]
  sub_ctcfBS_dt <- all_ctcf_dt[all_ctcf_dt$chromo == curr_chromo,]
  stopifnot(nrow(sub_ctcf_dt) > 0)
  stopifnot(nrow(sub_ctcfBS_dt) > 0)
  # changed here for v3_v2:
  
  matching_ctcfBS_dt <- sub_ctcfBS_dt[sub_ctcfBS_dt$start >= curr_minpos-aroundFamWindow &
                                        sub_ctcfBS_dt$end <= curr_maxpos+aroundFamWindow,]
  
  
  matching_ctcf_dt <- sub_ctcf_dt[sub_ctcf_dt$minPos >= curr_minpos-aroundFamWindow &
                sub_ctcf_dt$maxPos <= curr_maxpos+aroundFamWindow,]

  # final_fams_cluster_dt$closest_CTCFclust[i] <- paste0("CTCF_", sub_ctcf_dt$clusterID[i_minDist])
  final_fams_cluster_dt$meanSize_matchedCTCFclust[i] <- mean(matching_ctcf_dt$clusterSize)
  final_fams_cluster_dt$nbr_matchedCTCFclust[i] <- nrow(matching_ctcf_dt)
  final_fams_cluster_dt$nbr_matchedCTCF[i] <- nrow(matching_ctcfBS_dt)
  
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

    
  random_matching_ctcfBS_dt <- sub_ctcfBS_dt[sub_ctcfBS_dt$minPos >= randomStartPos-aroundFamWindow &
                                               sub_ctcfBS_dt$maxPos <= randomEndPos+aroundFamWindow,]
  
  
  random_matching_ctcf_dt <- sub_ctcf_dt[sub_ctcf_dt$minPos >= randomStartPos-aroundFamWindow &
                                    sub_ctcf_dt$maxPos <= randomEndPos+aroundFamWindow,]
  
  # final_fams_cluster_dt$closest_CTCFclust[i] <- paste0("CTCF_", sub_ctcf_dt$clusterID[i_minDist])
  random_meanSize_matchedCTCFclust <- mean(random_matching_ctcf_dt$clusterSize)
  random_nbr_matchedCTCFclust <- nrow(random_matching_ctcf_dt)
  random_nbr_matchedCTCF <- nrow(random_matching_ctcfBS_dt)
  
  
  
  list(
    obs_final_fams_cluster_dt=obs_final_fams_cluster_dt,
    random_meanSize_matchedCTCFclust=random_meanSize_matchedCTCFclust,
    random_nbr_matchedCTCFclust=random_nbr_matchedCTCFclust,
    random_nbr_matchedCTCF=random_nbr_matchedCTCF
  )
}

names(obsAndPermut_closest_clusters_dt) <- final_fams_cluster_dt$clusterID
                                                       
outFile <- file.path(outFolder, "obsAndPermut_closest_clusters_dt.Rdata")
save(obsAndPermut_closest_clusters_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

# load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3_V2/ctcf10000_family260000/obsAndPermut_closest_clusters_dt.Rdata")

all_obs_dt <- do.call('rbind', lapply(obsAndPermut_closest_clusters_dt, function(x) x[["obs_final_fams_cluster_dt"]]))

all_permutMeanSize <- unlist(lapply(obsAndPermut_closest_clusters_dt, function(x) x[["random_meanSize_matchedCTCFclust"]]))
all_permutMeanSize <- as.numeric(all_permutMeanSize)
# stopifnot(!is.na(all_permutMeanSize))

all_permutNclust <- unlist(lapply(obsAndPermut_closest_clusters_dt, function(x) x[["random_nbr_matchedCTCFclust"]]))
all_permutNclust <- as.numeric(all_permutNclust)
# stopifnot(!is.na(all_permutNclust))
stopifnot(length(all_permutNclust) == length(all_permutMeanSize))


all_permutNbs <- unlist(lapply(obsAndPermut_closest_clusters_dt, function(x) x[["random_nbr_matchedCTCF"]]))
all_permutNbs <- as.numeric(all_permutNbs)
# stopifnot(!is.na(all_permutNclust))
stopifnot(length(all_permutNclust) == length(all_permutNbs))


plot_dt <- data.frame(
  ctcfClust_meanSize=c(all_obs_dt$meanSize_matchedCTCFclust, all_permutMeanSize),
  dataType = c( rep("obs.", nrow(all_obs_dt)), rep("random", length(all_permutMeanSize))),
  stringsAsFactors = FALSE
)

plot2_dt <- data.frame(
  ctcfClust_nClust=c(all_obs_dt$nbr_matchedCTCFclust, all_permutNclust),
  dataType = c( rep("obs.", nrow(all_obs_dt)), rep("random", length(all_permutNclust))),
  stringsAsFactors = FALSE
)


plot3_dt <- data.frame(
  ctcfClust_nBS=c(all_obs_dt$nbr_matchedCTCF, all_permutNbs),
  dataType = c( rep("obs.", nrow(all_obs_dt)), rep("random", length(all_permutNbs))),
  stringsAsFactors = FALSE
)


plot_var <- "ctcfClust_meanSize"

plotTit <- paste0(plot_var, " - matched CTCF cluster")
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
outFile <- file.path(outFolder, paste0(plot_var, "_matchedCTCF_obsPermut_densityplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))


plot_var <- "ctcfClust_nClust"

plot2_dt <- plot2_dt[plot2_dt$ctcfClust_nClust<=quantile(plot2_dt$ctcfClust_nClust, probs = 0.95),]

plotTit <- paste0(plot_var, " - matched CTCF cluster")
subTit <- paste0("# ", names(table(plot2_dt$dataType)), "=", as.numeric(table(plot2_dt$dataType)), collapse="; ") 

p <- ggdensity(plot2_dt,
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
outFile <- file.path(outFolder, paste0(plot_var, "_matchedCTCF_obsPermut_densityplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))

plot_var <- "ctcfClust_nBS"

plotTit <- paste0(plot_var, " - matched CTCF BS")
subTit <- paste0("# ", names(table(plot3_dt$dataType)), "=", as.numeric(table(plot3_dt$dataType)), collapse="; ") 

plot3_dt <- plot3_dt[plot3_dt$ctcfClust_nBS<=quantile(plot3_dt$ctcfClust_nBS, probs = 0.95),]

p <- ggdensity(plot3_dt,
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
outFile <- file.path(outFolder, paste0(plot_var, "_matchedCTCFbs_obsPermut_densityplot.", plotTypeGG))
ggsave(p, filename=outFile, height=ggHeight, width=ggWidth)
cat(paste0("... written: ", outFile,  "\n"))


##########################################################################

signifThresh <- 0.01
signif_gm_dt <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878/PSEUDO_FINAL_TABLE/all_result_dt.Rdata"))
signif_gm_dt$region <- as.character(signif_gm_dt$region)
signif_tads <- signif_gm_dt$region[signif_gm_dt$adjPvalComb <= signifThresh]


gm_tad_dt <- read.delim(file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878/GM12878_40kb/genes2tad/all_assigned_regions.txt"),
                        col.names = c("chromo", "region", "start", "end"), stringsAsFactors = FALSE, header=F)
stopifnot(is.numeric( gm_tad_dt$end))
stopifnot(is.numeric( gm_tad_dt$start))
gm_tad_dt <- gm_tad_dt[order(gm_tad_dt$chromo, gm_tad_dt$start, gm_tad_dt$end),]
gm_tad_dt$chromo <- as.character(gm_tad_dt$chromo)

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

outFile <- file.path(outFolder, "all_obs_dt.Rdata")
save(all_obs_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

all_obs_dt$signifTADs_ratio <- all_obs_dt$nSignifTADs/all_obs_dt$nTADs 
stopifnot(all_obs_dt$signifTADs_ratio <= 1)
stopifnot(all_obs_dt$signifTADs_ratio >= 0)

plotType <- "png"
plotCex <- 1.4
plotHeight <- plotWidth <- 400

outFile <- file.path(outFolder, paste0("meanSize_matchedCTCFclust_vs_signifTADs_ratio_densityplot.", plotType))
do.call(plotType, list(outFile, height=plotHeight, width=plotWidth))
densplot(y=all_obs_dt$meanSize_matchedCTCFclust,x=all_obs_dt$signifTADs_ratio,
         cex.main=plotCex,cex.lab=plotCex,cex.axis=plotCex,
     xlab ="signifTADs_ratio", ylab="meanSize_matchedCTCFclust")
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))


outFile <- file.path(outFolder, paste0("nbr_matchedCTCFclust_vs_signifTADs_ratio_densityplot.", plotType))
do.call(plotType, list(outFile, height=plotHeight, width=plotWidth))
densplot(y=all_obs_dt$nbr_matchedCTCFclust,x=all_obs_dt$signifTADs_ratio,
         cex.main=plotCex,cex.lab=plotCex,cex.axis=plotCex,
     xlab="signifTADs_ratio", ylab="nbr_matchedCTCFclust")
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))

outFile <- file.path(outFolder, paste0("nbr_matchedCTCF_vs_signifTADs_ratio_densityplot.", plotType))
do.call(plotType, list(outFile, height=plotHeight, width=plotWidth))
densplot(y=all_obs_dt$nbr_matchedCTCF,x=all_obs_dt$signifTADs_ratio,
         cex.main=plotCex,cex.lab=plotCex,cex.axis=plotCex,
         xlab="signifTADs_ratio", ylab="nbr_matchedCTCF")
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))

