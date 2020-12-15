
# Rscript ctcf_and_families.R


################# AJOUTER LE # TOT DANS LA LEGENDE !!!

source("ctcf_da_utils.R")

library(ggplot2)
library(ggsci)
library(ggpubr)
library(foreach)
library(doMC)
registerDoMC(40)


pthresh <- 0.01
signifCol <- "red"
notSignifCol <- "grey"
plotCex <- 1.2

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_vars_2 <- c("maxInClust", "meanInClust", "nConvergent", "totCTCF")
all_vars_1 <- c("nDiffFam")

plotType <- "png"
myHeight <- 400
myWidth <- 400


outFolder <- file.path("CTCF_AND_FAMILIES")
dir.create(outFolder, recursive = TRUE)

setDir <- "/media/electron"
setDir <- ""
hgnc_geneFamilyFile <- file.path(setDir, "/mnt/ed4/marie/family_data_2/hgnc_entrez_family.txt")
hgnc_geneFamilyDT <- read.delim(hgnc_geneFamilyFile, col.names=c("entrezID", "family"), header = F, stringsAsFactors = F)
hgnc_geneFamilyDT$entrezID <- as.character(hgnc_geneFamilyDT$entrezID)
hgnc_geneFamilyDT$family_short <- unlist(sapply(hgnc_geneFamilyDT$family, function(x) strsplit(x, "\\|")[[1]][1] ))
family_dt <- hgnc_geneFamilyDT[,c("entrezID", "family_short")]
family_dt <- na.omit(family_dt)
family_dt$entrezID <- as.character(family_dt$entrezID)

runFolder <- "."

buildTable <- F

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
all_obs_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_obs_exprds <- sapply(all_obs_hicds, function(x) list.files(file.path(pipFolder, x)))


inFolder <- "CREATE_FINAL_TABLE/"
inFile <- file.path(inFolder, "all_result_dt.Rdata")
result_dt <- get(load(inFile))

inFolder <- "CTCF_AND_DA_ALLDS"
inFile <- file.path(inFolder, "ctcf2tad_dt.Rdata")
ctcf2tad_dt <- get(load(inFile))
stopifnot(ctcf2tad_dt$hicds %in% all_obs_hicds)
totCTCF_dt <- aggregate(chr ~ region + hicds, data=ctcf2tad_dt, FUN=length)
totCTCF_dt <- na.omit(totCTCF_dt)
colnames(totCTCF_dt)[colnames(totCTCF_dt) == "chr"] <- "totCTCF"
stopifnot(!is.na(totCTCF_dt))

inFolder <- "CTCF_AND_DA_ALLDS"
inFile <- file.path(inFolder, "clustByTAD_dt.Rdata")
clustByTAD_dt <- get(load(inFile))
stopifnot(clustByTAD_dt$hicds %in% all_obs_hicds)
stopifnot(clustByTAD_dt$exprds %in% unlist(all_obs_exprds))


maxInClust_dt <- aggregate(nInClust ~ hicds+exprds+region, data = clustByTAD_dt, FUN=max)
colnames(maxInClust_dt)[colnames(maxInClust_dt) == "nInClust"] <- "maxInClust"
hicds_maxInClust_dt <- unique(maxInClust_dt[,c("hicds", "region", "maxInClust")])
stopifnot(!duplicated(file.path(hicds_maxInClust_dt$hicds, hicds_maxInClust_dt$region)))

meanInClust_dt <- aggregate(nInClust ~ hicds+exprds+region, data = clustByTAD_dt, FUN=mean)
colnames(meanInClust_dt)[colnames(meanInClust_dt) == "nInClust"] <- "meanInClust"
hicds_meanInClust_dt <- unique(meanInClust_dt[,c("hicds", "region", "meanInClust")])
stopifnot(!duplicated(file.path(hicds_meanInClust_dt$hicds, hicds_meanInClust_dt$region)))


hicds_nConvergent_dt <- unique(clustByTAD_dt[,c("hicds",  "region", "nConvergent")])
stopifnot(!duplicated(file.path(hicds_nConvergent_dt$hicds, hicds_nConvergent_dt$region)))

clust_dt <- merge(hicds_maxInClust_dt, merge(hicds_meanInClust_dt, hicds_nConvergent_dt, by=c("hicds", "region"),all=F),
                  by=c("hicds", "region"),all=F)

stopifnot(nrow(clust_dt) == nrow(hicds_meanInClust_dt))
stopifnot(nrow(hicds_maxInClust_dt) == nrow(hicds_meanInClust_dt))
stopifnot(nrow(hicds_maxInClust_dt) == nrow(hicds_nConvergent_dt))
stopifnot(!is.na(clust_dt))
stopifnot(!duplicated(file.path(clust_dt$hicds, clust_dt$region)))

outFile <- file.path(outFolder, "clust_dt.Rdata")
save(clust_dt, file=outFile, version=2)


if(buildTable){
  fam_dt <- foreach(hicds = all_obs_hicds, .combine='rbind') %dopar% {
    
    # hicds_dt <- foreach(exprds = all_obs_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat(paste0("... start ", hicds, " \n"))
      
      gene2tad_dt <- read.delim(file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt"), stringsAsFactors = FALSE, 
                                header=F, col.names = c("entrezID", "chromo", "start", "end", "region"))
      gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)
      
      # geneList <- get(load(file.path(runFolder, pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
      # 
      # stopifnot(geneList %in% gene2tad_dt$entrezID)
      # 
      # gene2tad_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% geneList,]
      # stopifnot(nrow(gene2tad_dt) == length(geneList))$
      
      # do not subset -> because otherwise would also need to subset genes with family annotation
      
      tadSize <- setNames(as.numeric(table(gene2tad_dt$region)), names(table(gene2tad_dt$region)))
      
      g_tad_fam_dt <- merge(gene2tad_dt,family_dt,by="entrezID", all=FALSE)
      
      nAnnot <- setNames(as.numeric(table(g_tad_fam_dt$region)), names(table(g_tad_fam_dt$region)))
      
      stopifnot(!duplicated(g_tad_fam_dt$entrezID))
 
      minAnnot <- 3
      keepTADs <- names(nAnnot)[nAnnot >= minAnnot]     
      
      sub_tad_fam_dt <- g_tad_fam_dt[g_tad_fam_dt$region %in% keepTADs,]
      sub_tad_fam_dt$hicds <- hicds
      nDiffFam_byTAD_dt <- aggregate(family_short~hicds+region, data=sub_tad_fam_dt, FUN=function(x) length(unique(x)))
      colnames(nDiffFam_byTAD_dt)[colnames(nDiffFam_byTAD_dt)=="family_short"] <- "nDiffFam"
      nDiffFam_byTAD_dt$region <- as.character(nDiffFam_byTAD_dt$region)
      nDiffFam_byTAD_dt$totSize <- tadSize[nDiffFam_byTAD_dt$region]
      nDiffFam_byTAD_dt$totAnnot <- nAnnot[nDiffFam_byTAD_dt$region]
      nDiffFam_byTAD_dt
  }
  outFile <- file.path(outFolder, "fam_dt.Rdata")
  save(fam_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,"\n"))
  
}else {
  outFile <- file.path(outFolder, "fam_dt.Rdata")
  fam_dt <- get(load(outFile))
}

all_dt <- merge(totCTCF_dt, merge(fam_dt, clust_dt, by=c("hicds", "region"), all=F), by=c("hicds", "region"), all=F)
outFile <- file.path(outFolder, "all_dt.Rdata")
save(all_dt, file=outFile, version=2)


##########################################
########################################## densplot
##########################################

for(plot_var1 in all_vars_1) {
  for(plot_var2 in all_vars_2) {
    
    my_x <- all_dt[,plot_var1]
    my_y <- all_dt[,plot_var2]
    
    outFile <- file.path(outFolder, paste0(plot_var2, "_vs_", plot_var1, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x=my_x,
      y=my_y,
      xlab = plot_var1,
      ylab = plot_var2,
      cex.main=plotCex,
      cex.lab=plotCex,
      cex.axis=plotCex
    )
    addCorr(x=my_x, y=my_y,legPos="topright", bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}

##########################################
########################################## color by signif
##########################################


signif_all_dt <- merge(all_dt, result_dt, by=c("region", "hicds"), all=FALSE)

signif_all_dt$signif_lab <- ifelse(signif_all_dt$adjPvalComb <= pthresh, "signif.", "not signif.")


for(plot_var1 in all_vars_1) {
  for(plot_var2 in all_vars_2) {
    
    my_x <- signif_all_dt[,plot_var1]
    my_y <- signif_all_dt[,plot_var2]
    
    outFile <- file.path(outFolder, paste0(plot_var2, "_vs_", plot_var1, "_bySignif.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(
      x=my_x,
      y=my_y,
      xlab = plot_var1,
      ylab = plot_var2,
      pch=16,
      cex=0.7,
      col=notSignifCol,
      cex.main=plotCex, 
      cex.lab=plotCex,
      cex.axis=plotCex
    )
    points(
      x=my_x[signif_all_dt$signif_lab == "signif."],
      y=my_y[signif_all_dt$signif_lab == "signif."],
      col = signifCol,
      pch=16,
      cex=0.7
    )
    legend("topright", bty="n", col = c(signifCol, notSignifCol), legend=c("signif.", "not signif."))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}






# all_obs_exprds <- "TCGAluad_mutKRAS_mutEGFR"
# all_obs_hicds <- "ENCSR489OCU_NCI-H460_40kb"
# inFolder <- "CTCF_AND_DA_ALLDS_CHECK1"


# nBreaks <- 100
# step_breaks <- 1/nBreaks
# 
# plotTypeGG <- "svg"
# ggHeight <- 5
# ggWidth <- 6

