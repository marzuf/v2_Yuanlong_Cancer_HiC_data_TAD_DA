

# Rscript check_familyCliques.R

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

require(doMC)
require(foreach)
registerDoMC(40)
require(reshape2)
require(igraph)

runFolder <- "."
pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
familyVar <- "hgnc_family_short"
minCliqueSize <- 3

nMaxSize <- 1

maxSameTAD <- 0.5
minGenes <- 3


plotType <- "svg"
myHeight <- 5
myWidth <- 7

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

inFolder <- file.path("PREP_FAMILYCLIQUES", nMaxSize)
stopifnot(dir.exists(inFolder))
outFolder <- file.path("CHECK_FAMILYCLIQUES", nMaxSize)
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

hicds = "Barutcu_MCF-10A_40kb"
exprds="TCGAbrca_lum_bas"
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]

buildData <- FALSE




maxCliqueDist_dt <- foreach(hicds = all_hicds) %do%{
  cat(paste0("... start: ", hicds, "\n"))
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
    cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
    
    # retrieve file 
    famMod_file <- file.path(inFolder, hicds, exprds, "all_fams_dt.Rdata")
    if(!file.exists(famMod_file)) cat(famMod_file, "\n")
    stopifnot(file.exists(famMod_file))
    fam_data <- get(load(famMod_file))
    
    all_maxCliqueDist <- as.numeric(unlist(lapply(fam_data, function(x) x[["maxCliqueDist"]])))
    all_maxPairDist <- unique(as.numeric(unlist(lapply(fam_data, function(x) x[["maxPairDist"]]))))
    
    
    list(
      hicds=hicds,
      exprds=exprds,
      maxCliqueDist=all_maxCliqueDist,
      maxPairDist=all_maxPairDist,
      stringsAsFactors = FALSE
    )
  }
  exprds_dt
}

# maxPairDist_dt <- foreach(hicds = all_hicds, .combine='rbind') %do%{
#   cat(paste0("... start: ", hicds, "\n"))
#   exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
#     cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
#     
#     # retrieve file 
#     famMod_file <- file.path(inFolder, hicds, exprds, "all_fams_dt.Rdata")
#     if(!file.exists(famMod_file)) cat(famMod_file, "\n")
#     stopifnot(file.exists(famMod_file))
#     fam_data <- get(load(famMod_file))
#     
#     
#     data.frame(
#       hicds=hicds,
#       exprds=exprds,
#       maxPairDist=all_maxPairDist,
#       stringsAsFactors = FALSE
#     )
#   }
# }
# 

nDS <- length(unique(file.path(maxCliqueDist_dt$hicds, maxCliqueDist_dt$exprds)))
stopifnot(nDS == length(unique(file.path(maxCliqueDist_dt$hicds, maxCliqueDist_dt$exprds))))

maxCliqueDist_dt <- unlist(maxCliqueDist_dt, recursive = FALSE)
nDS <- length(maxCliqueDist_dt)

outFile <- file.path(outFolder,  paste0("allDS_distCliquesGenes_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(maxPairDist = unlist(lapply(maxCliqueDist_dt, function(x)x[["maxPairDist"]])),
       maxCliqueDist =  unlist(lapply(maxCliqueDist_dt, function(x)x[["maxCliqueDist"]]))),
  plotTit = paste0("all datasets -  n=", nDS)
)
mtext(side=3, text = paste0("minCliqueSize=", minCliqueSize), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))


outFile <- file.path(outFolder,  paste0("allDS_distCliquesGenes_density_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(maxPairDist_log10 = log10( unlist(lapply(maxCliqueDist_dt, function(x)x[["maxPairDist"]]))),
       maxCliqueDist_log10 = log10( unlist(lapply(maxCliqueDist_dt, function(x)x[["maxCliqueDist"]])))),
  plotTit = paste0("all datasets -  n=", nDS)
)
mtext(side=3, text = paste0("minCliqueSize=", minCliqueSize), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))



inFolder4 <- file.path("FAMILYCLIQUES_RUNTADMEANCORRRATIODOWN/")
outFile <- file.path(inFolder4, nMaxSize, "all_meanCorr_ratioDown.Rdata")
all_meanCorr_rD <- get(load(outFile))


cliqueFeatures_dt <- foreach(hicds = all_hicds, .combine='rbind') %do%{
  cat(paste0("... start: ", hicds, "\n"))
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
    

    
    
    # retrieve file - corr data
    dsCorrRd_file <- file.path(inFolder4, nMaxSize, hicds, exprds, "all_meanCorr_ratioDown_famCls.Rdata")
    stopifnot(file.exists(dsCorrRd_file))
    dsCorrRd_data <- get(load(dsCorrRd_file))
    
    toKeep <- unlist(lapply(dsCorrRd_data, function(x) is.numeric(x [["ratioDown"]])))
    stopifnot(length(toKeep) == length(dsCorrRd_data))
    rD_data <- dsCorrRd_data[toKeep]
    nGenes_rD_all <- unlist(lapply(rD_data, function(x) length(x[["fc_keptGenes"]])))
    stopifnot(nGenes_rD_all >= minGenes)
    nDiscSameTAD_rd <- sum(unlist(lapply(dsCorrRd_data, function(x)x[["ratioDown"]] == paste0("sameTAD>", maxSameTAD)) ))
    nDiscMinGenes_rd <- sum(unlist(lapply(dsCorrRd_data, function(x)x[["ratioDown"]] ==  paste0("<", minGenes, "genes") )))
    nCliquesWithValues_rd <- sum(toKeep)
    stopifnot(nCliquesWithValues_rd + nDiscMinGenes_rd + nDiscSameTAD_rd == length(dsCorrRd_data))
    meanGenesByClique_rd <- mean(unlist(nGenes_rD_all)) 
    
    
    toKeep <- unlist(lapply(dsCorrRd_data, function(x) is.numeric(x [["meanCorr"]])))
    stopifnot(length(toKeep) == length(dsCorrRd_data))
    corr_data <- dsCorrRd_data[toKeep]
    nGenes_corr_all <- unlist(lapply(corr_data, function(x) length(x[["corr_keptGenes"]])))
    stopifnot(nGenes_corr_all >= minGenes)
    nDiscSameTAD_corr <- sum(unlist(lapply(dsCorrRd_data, function(x)x[["meanCorr"]] == paste0("sameTAD>", maxSameTAD)) ))
    nDiscMinGenes_corr <- sum(unlist(lapply(dsCorrRd_data, function(x)x[["meanCorr"]] ==  paste0("<", minGenes, "genes") )))
    nCliquesWithValues_corr <- sum(toKeep)
    stopifnot(nCliquesWithValues_corr + nDiscMinGenes_corr + nDiscSameTAD_corr == length(dsCorrRd_data))
    meanGenesByClique_corr <- mean(unlist(nGenes_corr_all)) 
    
    
    

    
    data.frame(
      hicds=hicds,
      exprds=exprds,
      
      nDiscSameTAD_rd=  nDiscSameTAD_rd,
      nDiscMinGenes_rd=nDiscMinGenes_rd,
      nCliquesWithValues_rd=nCliquesWithValues_rd,
      meanGenesByClique_rd=meanGenesByClique_rd,
      
      nDiscSameTAD_corr=  nDiscSameTAD_corr,
      nDiscMinGenes_corr=nDiscMinGenes_corr,
      nCliquesWithValues_corr=nCliquesWithValues_corr,
      meanGenesByClique_corr=meanGenesByClique_corr,
      
      stringsAsFactors = FALSE
    )
  }
}
outFile <- file.path(outFolder, "cliqueFeatures_dt.Rdata")
save(cliqueFeatures_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile,"\n"))
# all(cliqueFeatures_dt$nDiscSameTAD_corr == cliqueFeatures_dt$nDiscSameTAD_fc)
# all(cliqueFeatures_dt$nDiscMinGenes_fc == cliqueFeatures_dt$nDiscMinGenes_corr)
load("CHECK_FAMILYCLIQUES/1/cliqueFeatures_dt.Rdata")
# load("cliqueFeatures_dt.Rdata")
require(reshape2)
plot_dt <- melt(cliqueFeatures_dt, id=c("hicds", "exprds"))
require(ggpubr)

plot_dt$variable <- factor(as.character(plot_dt$variable),
                           levels = c("nDiscSameTAD_rd", "nDiscSameTAD_corr",
                                      "nDiscMinGenes_rd", "nDiscMinGenes_corr",
                                      "nCliquesWithValues_rd", "nCliquesWithValues_corr",
                                      "meanGenesByClique_rd", "meanGenesByClique_corr"
                                      ))
stopifnot(!is.na(plot_dt$variable))

nfeatures_p <- ggboxplot(data=plot_dt[!grepl("mean", plot_dt$variable),], x="variable", y="value") +
  labs(x="", y="# of features")+
  theme(
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)
  )

outFile <- file.path(outFolder, "nFeatures_noMean.svg")
ggsave(nfeatures_p, file=outFile, height=7, width=8)
cat(paste0("... written: ", outFile,"\n"))


nfeatures_p <- ggboxplot(data=plot_dt[grepl("mean", plot_dt$variable),], x="variable", y="value") +
  labs(x="", y="# of features")+
  theme(
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)
  )

outFile <- file.path(outFolder, "nFeatures_onlyMean.svg")
ggsave(nfeatures_p, file=outFile, height=7, width=8)
cat(paste0("... written: ", outFile,"\n"))


cat(paste0("*** DONE: ", script_name, "\n"))
stop("-ok")

dsCorr_data <- get(load("FAMILYCLIQUES_RUNMEANTADCORR//1/Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/all_meanCorr_famCls.Rdata"))

# > dsCorr_data[[5]]
# $meanCorr
# [1] 0.0758462
# 
# $keptGenes
#   54583   56950   64754 
# "54583" "56950" "64754" 
# 
# $keptTADs
# [1] "chr1_TAD821" "chr1_TAD894" "chr1_TAD964"

dsCli_data <- get(load("PREP_FAMILYCLIQUES/1/Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/all_fams_dt.Rdata"))
# List of 3
# $ maxPairDist  : num 80998502
# $ maxCliqueDist: num 31809122
# $ fam_cl_dt    :'data.frame':	3 obs. of  2 variables:
#   ..$ entrezID: chr [1:3] "54583" "56950" "64754"
# ..$ clique  : chr [1:3] "Zinc fingers MYND-type_cl1" "Zinc fingers MYND-type_cl1" "Zinc fingers MYND-type_cl1"


######## TRASH

# inFolder2 <- file.path("WRONG_FAMILYCLIQUES_RUNMEANTADCORR/")
# outFile <- file.path(inFolder2, nMaxSize, "all_corr.Rdata")
# all_corr <- get(load(outFile))
# 
# inFolder3 <- file.path("WRONG_FAMILYCLIQUES_RUNTADRATIODOWN/")
# outFile <- file.path(inFolder3, nMaxSize, "all_ratioDown.Rdata")
# all_fc <- get(load(outFile))

# 
# # retrieve file - fc data
# dsFC_file <- file.path(inFolder3, nMaxSize, hicds, exprds, "all_ratioDown_famCls.Rdata")
# stopifnot(file.exists(dsFC_file))
# dsFC_data <- get(load(dsFC_file))
# fcOnly_data <- Filter(function(x) length(x) == 3, dsFC_data)
# nGenes_fc_all <- lapply(fcOnly_data, function(x) length(x[["keptGenes"]]))
# nDiscSameTAD_fc <- sum(dsFC_data == paste0("sameTAD>", maxSameTAD))
# nDiscMinGenes_fc <- sum(dsFC_data == paste0("<", minGenes, "genes"))
# nCliquesWithValues_fc <- sum(lengths(dsFC_data) == 3)
# stopifnot(nCliquesWithValues_fc + nDiscMinGenes_fc + nDiscSameTAD_fc == length(dsFC_data))
# meanGenesByClique_fc <- mean(unlist(nGenes_fc_all)) 
# 
# # retrieve file - corr data
# dsCorr_file <- file.path(inFolder2, nMaxSize, hicds, exprds, "all_meanCorr_famCls.Rdata")
# stopifnot(file.exists(dsCorr_file))
# dsCorr_data <- get(load(dsCorr_file))
# corrOnly_data <- Filter(function(x) length(x) == 3, dsCorr_data)
# nGenes_corr_all <- lapply(corrOnly_data, function(x) length(x[["keptGenes"]]))
# nDiscSameTAD_corr <- sum(dsCorr_data == paste0("sameTAD>", maxSameTAD))
# nDiscMinGenes_corr <- sum(dsCorr_data == paste0("<", minGenes, "genes"))
# nCliquesWithValues_corr <- sum(lengths(dsCorr_data) == 3)
# stopifnot(nCliquesWithValues_corr + nDiscMinGenes_corr + nDiscSameTAD_corr == length(dsCorr_data))
# meanGenesByClique_corr <- mean(unlist(nGenes_corr_all)) 
