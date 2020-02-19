require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

startTime <- Sys.time()

script_name <- "obs_permut_enhancerTADs_minDist.R"

# Rscript obs_permut_enhancerTADs_minDist.R

plotType <- "svg"
myHeight <- 7
myWidth <- 9
myHeightGG <- 7
myWidthGG <- 9

plotCex <- 1.4

outFolder <- "OBS_PERMUT_ENHANCERTADS_MINDIST"
dir.create(outFolder, recursive = TRUE)


final_DT <- get(load(file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")))
final_DT_permut <- get(load(file.path("CREATE_FINAL_TABLE_RANDOM", "all_result_dt.Rdata")))
result_dt <- rbind(final_DT, final_DT_permut)

exprds <- "TCGAluad_norm_luad"
famType <- "hgnc_family_short"

tad_signif_thresh <- 0.01

myhicds <- "ENCSR489OCU_NCI-H460"

all_hicds <- list.files(file.path("PIPELINE/OUTPUT_FOLDER"))

myhicds <- "ENCSR489OCU_NCI-H460"

exprds <- "TCGAluad_norm_luad"

buildData <- TRUE

all_hicds <- all_hicds[grep(myhicds, all_hicds)]

enh2entrez_dt <- get(load("enhanceratlas_data/PREP_TARGETGENES/enhancerEntrezTarget_dt.Rdata"))

if(buildData) {
  
  all_data_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    
    cat(paste0("... start ", hicds, "\n"))
    
    enh2entrez_g2t_dt <- enh2entrez_dt
    colnames(enh2entrez_g2t_dt)[colnames(enh2entrez_g2t_dt) == "target_entrezID"] <- "entrezID"
    enh2entrez_g2t_dt$entrezID <- as.character(enh2entrez_g2t_dt$entrezID)
    
    
    g2t_dt_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "gene_chromo", "gene_start", "gene_end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
    sum(enh2entrez_g2t_dt$entrezID %in% g2t_dt$entrezID)/nrow(enh2entrez_g2t_dt)
    
    
    enh2entrez_g2t_dt <- merge(enh2entrez_g2t_dt, g2t_dt[, c("entrezID", "gene_chromo", "gene_start", "gene_end", "region")], all.x=FALSE, all.y=FALSE, by="entrezID")
    colnames(enh2entrez_g2t_dt)[colnames(enh2entrez_g2t_dt) == "region"] <- "target_region"
    stopifnot(!is.na(enh2entrez_g2t_dt))
    cat(paste0(nrow(enh2entrez_dt) , "->", nrow(enh2entrez_g2t_dt), "\n"))
    
    enh2tad_dt <- get(load(file.path("enhanceratlas_data/PREP_ENHANCER2TAD_EP", paste0(hicds, "_enhancer2tad_dt.Rdata"))))
    colnames(enh2tad_dt)[colnames(enh2tad_dt) == "region"] <- "enhancer_region"
    
    
    enh2tad_dt$enhancer_id <- file.path(enh2tad_dt$enhancer_chromo, enh2tad_dt$enhancer_start,enh2tad_dt$enhancer_end)
    enh2entrez_g2t_dt$enhancer_id <- file.path(enh2entrez_g2t_dt$enhancer_chromo, enh2entrez_g2t_dt$enhancer_start,enh2entrez_g2t_dt$enhancer_end)
    stopifnot(enh2entrez_g2t_dt$enhancer_id %in% enh2tad_dt$enhancer_id)
    
    enh2tad_g2t_dt <- merge(enh2entrez_g2t_dt, unique(na.omit(enh2tad_dt[,c("enhancer_id", "enhancer_region")])), by=c("enhancer_id"))
    
    stopifnot(nrow(enh2tad_g2t_dt) == nrow(enh2entrez_g2t_dt))
    
    stopifnot(!is.na(enh2tad_g2t_dt))
    
    ### ADDED HERE: TAKE ONLY THE PAIRS THAT ARE 200 kb dist apart
    enh2tad_g2t_dt$enhancer_gene_dist <- abs(enh2tad_g2t_dt$gene_start - enh2tad_g2t_dt$enhancer_start)
    enh2tad_g2t_dt$enhancer_gene_dist[enh2tad_g2t_dt$enhancer_chromo != enh2tad_g2t_dt$gene_chromo ] <- NA
    
    
    enh2tad_g2t_dt$hicds <- hicds
    
    enh2tad_g2t_dt$exprds <- exprds
    
    enh2tad_g2t_dt
    
    # plot(density(na.omit(enh2tad_g2t_dt$enhancer_gene_dist)))
    # 
    # 
    # enh2tad_g2t_dt$target_entrez_sameTAD <- enh2tad_g2t_dt$enhancer_region == enh2tad_g2t_dt$target_region
    # 
    # tad_nEPsameTAD_dt <- aggregate(target_entrez_sameTAD ~ target_region, data=enh2tad_g2t_dt, FUN=function(x)sum(x))
    # colnames(tad_nEPsameTAD_dt) [colnames(tad_nEPsameTAD_dt) == "target_entrez_sameTAD"] <- "nEPsameTAD"
    # tad_nEPdiffTAD_dt <- aggregate(target_entrez_sameTAD ~ target_region, data=enh2tad_g2t_dt, FUN=function(x)sum(!x))
    # colnames(tad_nEPdiffTAD_dt) [colnames(tad_nEPdiffTAD_dt) == "target_entrez_sameTAD"] <- "nEPdiffTAD"
    # tad_nEP_dt <- aggregate(target_entrez_sameTAD ~ target_region, data=enh2tad_g2t_dt, FUN=length)
    # colnames(tad_nEP_dt) [colnames(tad_nEP_dt) == "target_entrez_sameTAD"] <- "nEP"
    # 
    # tad_nEP_all_dt <- merge(tad_nEP_dt, merge(tad_nEPdiffTAD_dt, tad_nEPsameTAD_dt, by="target_region"), by="target_region")
    # stopifnot(nrow(tad_nEP_all_dt) == length(unique(enh2tad_g2t_dt$target_region)))
    # stopifnot(!is.na(tad_nEP_all_dt))
    # 
    # colnames(tad_nEP_all_dt)[colnames(tad_nEP_all_dt)=="target_region"] <- "region"
    # 
    # if(hicds %in% final_DT$hicds) {
    #   result_dt <- final_DT[final_DT$hicds == hicds & final_DT$exprds == exprds,]
    # } else if(hicds %in% final_DT_permut$hicds) {
    #   result_dt <- final_DT_permut[final_DT_permut$hicds == hicds & final_DT_permut$exprds == exprds,]
    # } else {
    #   stop("error")
    # }
    # 
    # tad_nEP_all_dt_final <- merge(tad_nEP_all_dt, result_dt[,c("hicds", "exprds", "region", "meanCorr", "adjPvalComb")],by="region", all.x=FALSE, all.y=FALSE)
    # stopifnot(!is.na(tad_nEP_all_dt_final))
    # 
    # cat(paste0(nrow(result_dt) , " -> ", nrow(tad_nEP_all_dt_final), "\n"))
    # tad_nEP_all_dt_final
    
    
  }
  outFile <- file.path(outFolder, "all_data_dt.Rdata")
  save(all_data_dt, file =outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
} else{
  
  outFile <- file.path(outFolder, "all_data_dt.Rdata")
  all_data_dt <- get(load(outFile))
}


all_data_dt$hicds_lab <- gsub(myhicds, "", all_data_dt$hicds)

outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_EP_dist.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(log10(na.omit(all_data_dt)$enhancer_gene_dist), na.omit(all_data_dt)$hicds_lab),
               plotTit="Dist enhancer_start-gene_start [log10]", legPos = "topleft")
mtext(side=3, text = paste0(myhicds, " -", exprds, "; NA diff. chromo all DS: ", sum(is.na(all_data_dt$enhancer_gene_dist))))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



all_data_dt$target_entrez_sameTAD <- all_data_dt$enhancer_region == all_data_dt$target_region
all_data_dt <- all_data_dt[grepl("_TAD", all_data_dt$target_region),]


#################################################################################################################################### all

tad_nEPsameTAD_dt <- aggregate(target_entrez_sameTAD ~ target_region+hicds+hicds_lab, data=all_data_dt, FUN=function(x)sum(x))
colnames(tad_nEPsameTAD_dt) [colnames(tad_nEPsameTAD_dt) == "target_entrez_sameTAD"] <- "nEPsameTAD"
tad_nEPdiffTAD_dt <- aggregate(target_entrez_sameTAD ~ target_region+hicds+hicds_lab, data=all_data_dt, FUN=function(x)sum(!x))
colnames(tad_nEPdiffTAD_dt) [colnames(tad_nEPdiffTAD_dt) == "target_entrez_sameTAD"] <- "nEPdiffTAD"
tad_nEP_dt <- aggregate(target_entrez_sameTAD ~ target_region+hicds+hicds_lab, data=all_data_dt, FUN=length)
colnames(tad_nEP_dt) [colnames(tad_nEP_dt) == "target_entrez_sameTAD"] <- "nEP"

outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nEP_byTAD.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(log10(tad_nEP_dt$nEP), tad_nEP_dt$hicds_lab),
               plotTit="# EP by TAD [log10]", legPos = "topright")
mtext(side=3, text = paste0(myhicds, " -", exprds))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nEPsameTAD_byTAD.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(log10(tad_nEPsameTAD_dt$nEPsameTAD), tad_nEPsameTAD_dt$hicds_lab),
               plotTit="# EP same TAD by TAD [log10]", legPos = "topright")
mtext(side=3, text = paste0(myhicds, " -", exprds))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nEPdiffTAD_byTAD.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(log10(tad_nEPdiffTAD_dt$nEPdiffTAD), tad_nEPdiffTAD_dt$hicds_lab),
               plotTit="# EP diff TAD by TAD [log10]", legPos = "topright")
mtext(side=3, text = paste0(myhicds, " -", exprds))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

ratio_dt <- merge(tad_nEPsameTAD_dt, tad_nEPdiffTAD_dt, by = c("target_region", "hicds", "hicds_lab"))
ratio_dt$nEP_sameTAD_over_diffTAD <- (ratio_dt$nEPsameTAD+1)/(ratio_dt$nEPdiffTAD+1)                  

outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_ratio_nEPsamediffTAD_byTAD.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(log10(ratio_dt$nEP_sameTAD_over_diffTAD), ratio_dt$hicds_lab),
               plotTit="ratio # EP same over diff TAD by TAD [log10]", legPos = "topright")
mtext(side=3, text = paste0(myhicds, " -", exprds))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
stop("ok")


#################################################################################################################################### 

# tad_signif_thresh <- 0.01
# colnames(all_data_dt)[colnames(all_data_dt) == "target_region"] <- "region"
# data_result_dt <- merge(all_data_dt, result_dt, by=c("hicds", "exprds", "region" ), all.x=FALSE, all.y=FALSE)
# data_result_dt$signif <- ifelse(data_result_dt$adjPvalComb <= tad_signif_thresh, "signif.", "not signif.")


#################################################################################################################################### signif > NOT ENOUGH DATA POINTS

# signif_dt <- data_result_dt[data_result_dt$signif == "signif.",]
# 
# 
# tad_nEPsameTAD_dt <- aggregate(target_entrez_sameTAD ~ region+hicds+hicds_lab, data=signif_dt, FUN=function(x)sum(x))
# colnames(tad_nEPsameTAD_dt) [colnames(tad_nEPsameTAD_dt) == "target_entrez_sameTAD"] <- "nEPsameTAD"
# tad_nEPdiffTAD_dt <- aggregate(target_entrez_sameTAD ~ region+hicds+hicds_lab, data=signif_dt, FUN=function(x)sum(!x))
# colnames(tad_nEPdiffTAD_dt) [colnames(tad_nEPdiffTAD_dt) == "target_entrez_sameTAD"] <- "nEPdiffTAD"
# tad_nEP_dt <- aggregate(target_entrez_sameTAD ~ region+hicds+hicds_lab, data=signif_dt, FUN=length)
# colnames(tad_nEP_dt) [colnames(tad_nEP_dt) == "target_entrez_sameTAD"] <- "nEP"
# 
# outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nEP_byTAD_signifOnly.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot_multiDens(split(log10(tad_nEP_dt$nEP), tad_nEP_dt$hicds_lab),
#                plotTit="# EP by TAD [log10]", legPos = "topright")
# mtext(side=3, text = paste0(myhicds, " -", exprds, "; only signif. TADs p-val <= ", tad_signif_thresh))
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nEPsameTAD_byTAD_minDist.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot_multiDens(split(log10(tad_nEPsameTAD_dt$nEPsameTAD[!grepl("PERMUTG2T", tad_nEPsameTAD_dt$hicds_lab)]), tad_nEPsameTAD_dt$hicds_lab[!grepl("PERMUTG2T", tad_nEPsameTAD_dt$hicds_lab)]),
#                plotTit="# EP same TAD by TAD [log10]", legPos = "topright")
# mtext(side=3, text = paste0(myhicds, " -", exprds, "; E-G dist >= ", minDist/1000, " kb"))
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nEPdiffTAD_byTAD_minDist.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot_multiDens(split(log10(tad_nEPdiffTAD_dt$nEPdiffTAD[!grepl("PERMUTG2T", tad_nEPdiffTAD_dt$hicds_lab)]), tad_nEPdiffTAD_dt$hicds_lab[!grepl("PERMUTG2T", tad_nEPdiffTAD_dt$hicds_lab)]),
#                plotTit="# EP diff TAD by TAD [log10]", legPos = "topright")
# mtext(side=3, text = paste0(myhicds, " -", exprds, "; E-G dist >= ", minDist/1000, " kb"))
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))




################################################################################################# SAME WITH MIN DIST

#meanTADsize*0.5 = 128642.2
minDist <- 130*1000

minDist_dt <- all_data_dt[!is.na(all_data_dt$enhancer_gene_dist),]
minDist_dt <- minDist_dt[minDist_dt$enhancer_gene_dist >= minDist,]




tad_nEPsameTAD_dt <- aggregate(target_entrez_sameTAD ~ target_region+hicds+hicds_lab, data=minDist_dt, FUN=function(x)sum(x))
colnames(tad_nEPsameTAD_dt) [colnames(tad_nEPsameTAD_dt) == "target_entrez_sameTAD"] <- "nEPsameTAD"
tad_nEPdiffTAD_dt <- aggregate(target_entrez_sameTAD ~ target_region+hicds+hicds_lab, data=minDist_dt, FUN=function(x)sum(!x))
colnames(tad_nEPdiffTAD_dt) [colnames(tad_nEPdiffTAD_dt) == "target_entrez_sameTAD"] <- "nEPdiffTAD"
tad_nEP_dt <- aggregate(target_entrez_sameTAD ~ target_region+hicds+hicds_lab, data=minDist_dt, FUN=length)
colnames(tad_nEP_dt) [colnames(tad_nEP_dt) == "target_entrez_sameTAD"] <- "nEP"

outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nEP_byTAD_minDist.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(log10(tad_nEP_dt$nEP), tad_nEP_dt$hicds_lab),
               plotTit="# EP by TAD [log10]", legPos = "topright")
mtext(side=3, text = paste0(myhicds, " -", exprds, "; E-G dist >= ", minDist/1000, " kb"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nEPsameTAD_byTAD_minDist.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(log10(tad_nEPsameTAD_dt$nEPsameTAD[!grepl("PERMUTG2T", tad_nEPsameTAD_dt$hicds_lab)]), tad_nEPsameTAD_dt$hicds_lab[!grepl("PERMUTG2T", tad_nEPsameTAD_dt$hicds_lab)]),
               plotTit="# EP same TAD by TAD [log10]", legPos = "topright")
mtext(side=3, text = paste0(myhicds, " -", exprds, "; E-G dist >= ", minDist/1000, " kb"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nEPdiffTAD_byTAD_minDist.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(log10(tad_nEPdiffTAD_dt$nEPdiffTAD[!grepl("PERMUTG2T", tad_nEPdiffTAD_dt$hicds_lab)]), tad_nEPdiffTAD_dt$hicds_lab[!grepl("PERMUTG2T", tad_nEPdiffTAD_dt$hicds_lab)]),
               plotTit="# EP diff TAD by TAD [log10]", legPos = "topright")
mtext(side=3, text = paste0(myhicds, " -", exprds, "; E-G dist >= ", minDist/1000, " kb"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))










