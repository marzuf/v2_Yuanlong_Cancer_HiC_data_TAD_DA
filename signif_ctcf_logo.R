library("Biostrings")
library("msa")
dnaSet = DNAStringSet(c("AACCTT","CCGGTTTT","AAAGGGTTT"))
res = msa(dnaSet,method="ClustalOmega")
print(res)

myctcf_seqs <- c(">>><<",">>>><","><><>>")
myctcf_seqs_hack <- gsub("<", "T", gsub(">", "A", myctcf_seqs))
dnaSet = DNAStringSet(myctcf_seqs_hack)
res = msa(dnaSet,method="ClustalOmega")
print(res)

aligned_ctcf_sequences <- as.character(res)

source("http://bioconductor.org/biocLite.R")
BiocManager::install("Logolas")
logomaker(aligned_ctcf_sequences, type = "Logo")

aligned_ctcf_sequences_hack <- gsub("T", "<", gsub("A", ">", aligned_ctcf_sequences))
logomaker(aligned_ctcf_sequences_hack, type = "Logo")


load("CTCF_AND_DA_ALLDS/ctcf2tad_dt.Rdata")
load("CREATE_FINAL_TABLE/all_result_dt.Rdata")

pthresh <- 0.01

signif_dt <- all_result_dt[all_result_dt$adjPvalComb <= pthresh,]

signif_ctcf2tad_dt <- merge(signif_dt,ctcf2tad_dt, by=c("hicds", "region"), all=FALSE , suffixes=c("_TAD", "_CTCF"))
signif_ctcf2tad_dt$chromo <- gsub("(chr.+)_.+", "\\1", signif_ctcf2tad_dt$region)
stopifnot(signif_ctcf2tad_dt$chromo %in% paste0("chr", 1:22))
signif_ctcf2tad_dt <- signif_ctcf2tad_dt[order(signif_ctcf2tad_dt$chromo, signif_ctcf2tad_dt$start_CTCF, signif_ctcf2tad_dt$end_CTCF),]

agg_dt <- aggregate(orientation~hicds + exprds + region, data=signif_ctcf2tad_dt, FUN=function(x) paste0(x, collapse=""))
agg_dt$regionID <- file.path(agg_dt$hicds, agg_dt$exprds, agg_dt$region)
nrow(agg_dt)
# 1442
length_agg_dt <- aggregate(orientation~hicds + exprds + region, data=signif_ctcf2tad_dt, FUN=length)
length_agg_dt$regionID <- file.path(length_agg_dt$hicds, length_agg_dt$exprds, length_agg_dt$region)
plot(density(length_agg_dt$orientation))
# take only those between 8 and 12
keepTADs <- length_agg_dt$regionID[length_agg_dt$orientation >= 8 & length_agg_dt$orientation <= 12]
length(keepTADs)
# 474

sub_agg_dt <- agg_dt[agg_dt$regionID %in% keepTADs,]
stopifnot(nrow(sub_agg_dt) == length(keepTADs))


myctcf_seqs <- sub_agg_dt$orientation
myctcf_seqs_hack <- gsub("<", "T", gsub(">", "A", myctcf_seqs))
dnaSet = DNAStringSet(myctcf_seqs_hack)
res = msa(dnaSet,method="ClustalOmega")
print(res)
aligned_ctcf_sequences <- as.character(res)
logomaker(aligned_ctcf_sequences, type = "Logo")

aligned_ctcf_sequences_hack <- gsub("T", "<", gsub("A", ">", aligned_ctcf_sequences))
logomaker(aligned_ctcf_sequences_hack, type = "Logo")


