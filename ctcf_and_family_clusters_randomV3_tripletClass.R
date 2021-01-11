# CTCF_CLUSTER/ctcf_bs_nanni2020.bed             CTCF_CLUSTER/gene2family.bed     
# CTCF_CLUSTER/meanSize.Rdata  CTCF_CLUSTER/test_ctcf_bs_nanni2020.bed

# Rscript ctcf_and_family_clusters_randomV3_tripletClass.R


require(doMC)
require(foreach)
registerDoMC(40)
require(ggpubr)
require(ggsci)
require(ggplot2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

myHeight <- 400
myWidth <- 400
plotType <- "png"

plotCex <- 1.2

chromoSize_dt <- read.delim("CTCF_CLUSTER/hg19_chromo_sizes.txt", stringsAsFactors=FALSE, col.names=c("chromo", "size"), header=F)

nPermut <- 100

ctcf_clustSize <- 4000
# ctcf_clustSize <- 10000
# ctcf_clustSize <- 20000
family_clustSize <- 260000
# family_clustSize <- 130000


outFolder <- file.path("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3_TRIPLETCLASS", paste0("ctcf", ctcf_clustSize, "_family", family_clustSize))
dir.create(outFolder, recursive = TRUE)
                       
tmpFolder <- file.path(outFolder, "TMP")
dir.create(tmpFolder, recursive = TRUE)

###########################
# prep CTCF data
###########################

# #../bedtools cluster -i ctcf_bs_nanni2020.bed -d 40000 > ctcf_bs_nanni2020_cluster4kb.bed

# ctcf_clustFile="CTCF_AND_FAMILY_CLUSTERS_RANDOMV3_TRIPLETCLASS/ctcf4000_family260000//ctcf_bs_nanni2020_cluster4000bp.bed"
ctcf_clustFile <- file.path(outFolder, paste0("ctcf_bs_nanni2020_cluster", ctcf_clustSize, "bp.bed"))

mycmd <- paste0("./bedtools cluster -i CTCF_CLUSTER/ctcf_bs_nanni2020_withTC.bed -d ", ctcf_clustSize, 
" > ", ctcf_clustFile)
cat(paste0("> ", mycmd, "\n"))
system(mycmd)

ctcf_cluster_dt <- read.delim(ctcf_clustFile, header=FALSE,
                              stringsAsFactors = FALSE, 
                              col.names = c("chromo", "start", "end", "orientation", "Triplet_class", "clusterID"))
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
# load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3_TRIPLETCLASS/ctcf4000_family260000/agg_ctcf_cluster_dt.Rdata")


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
# > allgenes_clustFile="CTCF_AND_FAMILY_CLUSTERS_RANDOMV3_TRIPLETCLASS/ctcf4000_family260000/TMP/all_genes_clust260000.bed"
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
#  load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3_TRIPLETCLASS/ctcf4000_family260000/all_fams_agg_cluster_dt.Rdata")

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


# for each CTCF cluster 
# -> dist to closest Fam
# -> size of the closest Fam
# -> length of the closest Fam
# -> size of current clust
# -> length of current clust
# -> triplet class dist of current clust
# -> motif of current clust


# stop("-ok\n")

final_ctcf_cluster_dt$clusterLength <- final_ctcf_cluster_dt$maxPos - final_ctcf_cluster_dt$minPos + 1

final_fams_cluster_dt$clusterLength <- final_fams_cluster_dt$maxPos - final_fams_cluster_dt$minPos + 1

ctcf_cluster_dt$Triplet_class <- as.character(ctcf_cluster_dt$Triplet_class)

ctcf_cluster_dt$Triplet_class <- factor(ctcf_cluster_dt$Triplet_class,
                                        levels=sort(unique(ctcf_cluster_dt$Triplet_class)))

stopifnot(is.numeric(ctcf_cluster_dt$end))
stopifnot(is.numeric(ctcf_cluster_dt$start))

ctcf_cluster_dt <- ctcf_cluster_dt[order(ctcf_cluster_dt$chromo, ctcf_cluster_dt$start, ctcf_cluster_dt$end),]

ctcf2fam_cluster_dt <- foreach(i = 1:nrow(final_ctcf_cluster_dt), .combine='rbind') %dopar% {
  
  
  curr_cluster <- final_ctcf_cluster_dt$clusterID[i]
  
  curr_chromo <- final_ctcf_cluster_dt$chromo[i]
  curr_midpos <- final_ctcf_cluster_dt$midPos[i]
  
  curr_length <- final_ctcf_cluster_dt$clusterLength[i]
  curr_size <- final_ctcf_cluster_dt$clusterSize[i]
  
  
  member_ctcf_dt <- ctcf_cluster_dt[ctcf_cluster_dt$clusterID == curr_cluster,]
  stopifnot(nrow(member_ctcf_dt) > 0)
  stopifnot(nrow(member_ctcf_dt)  == curr_size)
  
  stopifnot(diff(member_ctcf_dt$start)>=0) # check ordered by start pos
  
  curr_motif <- paste0(member_ctcf_dt$orientation, collapse="")
  
  curr_tcDist <- table(member_ctcf_dt$Triplet_class)/nrow(member_ctcf_dt)
  
  
  sub_fam_dt <- final_fams_cluster_dt[final_fams_cluster_dt$chromo == curr_chromo,]
  stopifnot(nrow(sub_fam_dt) > 0)
  
  i_closestFamClust <- which.min(abs(curr_midpos - sub_fam_dt$midPos)) 
  closestFamClust_distTo <- min(abs(curr_midpos - sub_fam_dt$midPos)) 
  closestFamClust_clustID <- sub_fam_dt$clusterSize[i_closestFamClust]
  closestFamClust_family <- sub_fam_dt$family[i_closestFamClust]
  closestFamClust_clusterSize <- sub_fam_dt$clusterSize[i_closestFamClust]
  closestFamClust_clusterLength <- sub_fam_dt$clusterLength[i_closestFamClust]

  data.frame(
    CTCFclust_ID=curr_cluster,
    CTCFclust_chromo=curr_chromo,
    CTCFclust_midpos=curr_midpos,
    CTCFclust_length=curr_length,
    CTCFclust_size=curr_size,
    CTCFclust_ratioC=as.numeric(curr_tcDist["C"]),
    CTCFclust_ratioCD=as.numeric(curr_tcDist["CD"]),
    CTCFclust_ratioD=as.numeric(curr_tcDist["D"]),
    CTCFclust_ratioS=as.numeric(curr_tcDist["S"]),
    CTCFclust_ratioT=as.numeric(curr_tcDist["T"]),
    closestFamClust_distTo=closestFamClust_distTo,
    closestFamClust_clustID=closestFamClust_clustID,
    closestFamClust_family=closestFamClust_family,
    CTCFclust_length=closestFamClust_clusterSize,
    closestFamClust_size=closestFamClust_clusterSize,
  stringsAsFactors=FALSE
    )
  
  
}
outFile <- file.path(outFolder, "ctcf2fam_cluster_dt.Rdata")
save(ctcf2fam_cluster_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

stopifnot(nrow(ctcf2fam_cluster_dt) == nrow(final_ctcf_cluster_dt))

# load("CTCF_AND_FAMILY_CLUSTERS_RANDOMV3_TRIPLETCLASS/ctcf4000_family260000/ctcf2fam_cluster_dt.Rdata")

# stop("-ok\n")
  

all_y <- c("CTCFclust_ratioC", "CTCFclust_ratioCD","CTCFclust_ratioD","CTCFclust_ratioS", "CTCFclust_ratioT")
xvar1 <- "closestFamClust_distTo"
xvar2 <- "CTCFclust_size"

for(curr_y in all_y) {
  
  
  myy <- ctcf2fam_cluster_dt[,curr_y]
  myx <- log10(ctcf2fam_cluster_dt[,xvar1])
  
  outFile <- file.path(outFolder, paste0(curr_y, "_vs_", xvar1, "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x=myx,
    y=myy,
    xlab=paste0(xvar1, " [log10]"),
    ylab=curr_y,
    cex.axis=plotCex,
    cex.lab=plotCex,
    cex.main=plotCex
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
  
  
  myx <- ctcf2fam_cluster_dt[,xvar2]
  
  outFile <- file.path(outFolder, paste0(curr_y, "_vs_", xvar2, "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  
  densplot(
    x=myx,
    y=myy,
    xlab=xvar2,
    ylab=curr_y,
    cex.axis=plotCex,
    cex.lab=plotCex,
    cex.main=plotCex
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
  
  

    
  
  
}


yvar <- "CTCFclust_size"
xvar <- "closestFamClust_distTo"
myy <- ctcf2fam_cluster_dt[,yvar]
myx <- log10(ctcf2fam_cluster_dt[,xvar])

outFile <- file.path(outFolder, paste0(yvar, "_vs_", xvar, "densplot.", plotType))

do.call(plotType, list(outFile, height=myHeight, width=myWidth))

densplot(
  x=myx,
  y=myy,
  xlab=paste0(xvar, " [log10]"),
  ylab=yvar,
  cex.axis=plotCex,
  cex.lab=plotCex,
  cex.main=plotCex
)
foo <- dev.off()
cat(paste0("... written: ", outFile,"\n"))















