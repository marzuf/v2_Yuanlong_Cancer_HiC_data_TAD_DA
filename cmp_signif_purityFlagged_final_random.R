options(scipen=100)

SSHFS=F

# Rscript cmp_signif_purityFlagged_final_random.R 

# _final:  discussion 04.08.2020 Giovanni - keep aran CPE data, first vial only > version for the RANDOM datasets

purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- paste0("Aran - ", pm)

script0_name <- "0_prepGeneData"

### HARD-CODED - MAIN SETTINGS

corMet <- "pearson"
transfExpr <- "log10"
signifThresh <- 0.01
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifcol <- paste0(signif_column, "_", signifThresh)


script_name <- "cmp_signif_purityFlagged_final_random.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(ggpubr)
require(patchwork)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

setDir <- ""
plotCex <- 1.2

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("plot_lolliTAD_funct.R")
# source("my_heatmap.2.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
# source("../MANUSCRIPT_FIGURES//settings.R")
# source("../MANUSCRIPT_FIGURES/full_dataset_names.R")


buildTable <- TRUE

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 5.5)
myWidth <- 5.5
axisCex <- 1.4


outFolder <- file.path("CMP_SIGNIF_PURITYFLAGGED_FINAL_RANDOM", purity_ds, pm, transfExpr)
dir.create(outFolder, recursive = TRUE)

# > setdiff(all_result_dt$dataset, obs_dt$dataset)
# [1] "GSE118588_Panc_beta_40kb/TCGApaad_wt_mutKRAS" "K562_40kb/TCGAlaml_wt_mutFLT3"               
# [3] "PA2_40kb/TCGApaad_wt_mutKRAS"                 "PA3_40kb/TCGApaad_wt_mutKRAS"                
# [5] "Panc1_rep12_40kb/TCGApaad_wt_mutKRAS"        

obs_dt <- get(load("SIGNIF_PURITYFLAGGED_FINAL/aran/CPE/log10/all_dt.Rdata"))
random_dt1 <- get(load("SIGNIF_PURITYFLAGGED_FINAL_RANDOMMIDPOS//aran/CPE/log10/all_dt.Rdata"))
random_dt2 <- get(load("SIGNIF_PURITYFLAGGED_FINAL_PERMUT//aran/CPE/log10/all_dt.Rdata"))
random_dt <- rbind(random_dt1, random_dt2)

obs_dt$ratioFlagged <- obs_dt$nPurityFlagged/obs_dt$nTot
obs_dt$nSignif_notFlagged <- obs_dt$nSignif - obs_dt$nSignifAndFlagged

obs_dt$ratioSignif_notFlagged <- 1 - obs_dt$ratioSignifFlagged

obs_dt$dotcols <- all_cols[all_cmps[basename(obs_dt$dataset)]]
obs_dt$dotpch <- ifelse(grepl("TCGAskcm_lowInf_highInf", obs_dt$dataset), 8, 16)

my_x <- obs_dt$ratioSignif
my_y <- obs_dt$ratioSignifFlagged

outfile <- file.path(outFolder, "ratioSignifFlagged_ratioSignif_withSkcm.svg")
svg(outfile, height=7,width=7)
plot(
  x= my_x,
  y= my_y,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex,
  main="",
  xlab="ratioSignif",
  ylab="ratioSignifFlagged",
  pch=obs_dt$dotpch,
  col=obs_dt$dotcols
)
text(y=my_y[grepl("TCGAskcm_lowInf_highInf", obs_dt$dataset)],
     x=my_x[grepl("TCGAskcm_lowInf_highInf", obs_dt$dataset)],
     labels=gsub("/", "\n", obs_dt$dataset[grepl("TCGAskcm_lowInf_highInf", obs_dt$dataset)]), cex=0.6
     )
addCorr(x=my_x, y=my_y, bty="n")
abline(lm(my_y~my_x), lty=2, col="darkgrey")
legend("topleft", col=all_cols, legend=names(all_cols), bty="n", pch=16, cex=0.8)

foo <- dev.off()

stop("-ok\n")


random_dt$ratioFlagged <- random_dt$nPurityFlagged/random_dt$nTot
random_dt$nSignif_notFlagged <- random_dt$nSignif - random_dt$nSignifAndFlagged

random_dt$ratioSignif_notFlagged <- 1 - random_dt$ratioSignifFlagged

all_patts <- c("RANDOMMIDPOS", "RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT", "PERMUTG2T")
pat=all_patts[1]
pat=all_patts[2]
pat="PERMUTG2T"
all_plotVars <- c("nSignif", "nSignif_notFlagged", "ratioSignif", "ratioSignif_notFlagged", "nSignifAndFlagged",  "ratioSignifFlagged")
plotVar <- "nSignif"
plotVar <- "ratioSignif_notFlagged" 

for(plotVar in all_plotVars){
  
  plotList <- list()
  
  
  cat(paste0("> START ", plotVar, "\n"))
  
  for(pat in all_patts) {
    
    cat(paste0("......", plotVar, " - ", pat, "\n"))
    
    ylab <- ifelse(plotVar == "nSignif", "# signif. TADs", 
                   ifelse(plotVar == "nSignif_notFlagged", "# signif. and not purity-flagged TADs",
                          ifelse(plotVar == "nSignifAndFlagged", "# signif. and purity-flagged TADs",
                          ifelse(plotVar == "ratioSignif", "ratio signif. TADs",
                                 ifelse(plotVar == "ratioSignifFlagged", "ratio signif. and purity-flagged TADs",
                                 ifelse(plotVar == "ratioSignif_notFlagged", "ratio signif. and not purity-flagged TADs", NA))))))
    
    stopifnot(!is.na(ylab))
    
    tmp_dt <- random_dt[grepl(paste0(pat,"_40kb" ), random_dt$dataset),]
    stopifnot(nrow(tmp_dt) > 0)
    
    tmp_dt$dataset <- gsub("_RANDOM.+_40kb|_PERMUT.+_40kb", "_40kb", tmp_dt$dataset)
    stopifnot(setequal(tmp_dt$dataset, obs_dt$dataset))
    
    plot_dt <- merge(tmp_dt, obs_dt, by="dataset", suffixes=c("_rd", "_obs"), all=TRUE)
    # stopifnot(!is.na(plot_dt))  # not TRUE because I have some division by 0
    
    
    colnames(plot_dt)[colnames(plot_dt) == paste0(plotVar, "_obs")] <- "observed"
    colnames(plot_dt)[colnames(plot_dt) == paste0(plotVar, "_rd")] <- paste0(pat)
    
    plot_dt <- plot_dt[,c("observed", pat)]
    if(nrow(na.omit(plot_dt)) > 0)      plot_dt <- na.omit(plot_dt)
    
    p <- ggpaired(data=plot_dt, cond1=paste0("observed"),  cond2=paste0(pat),
                  # = "supp", y = "len",
                  # color = "supp", 
                  # labels=c("observed", "random"),
                  title= paste0(plotVar, " - obs. vs. ", pat),
                  color="condition",
                  line.color = "gray", line.size = 0.4,
                  palette = "jco")+
      labs(x="", y=ylab, color="", subtitle = paste0("# datasets = ", nrow(plot_dt))) + 
      guides(color=FALSE)+
      stat_compare_means(paired = TRUE, method="wilcox.test", na.rm=TRUE) +
      theme(
        plot.title = element_text(size=16, hjust=0.5, face="bold"),
        plot.subtitle = element_text(size=14, hjust=0.5, face="italic"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        legend.text=element_text(size=12)
      )
    
    
    plotList[[pat]] <- p
    
    
  }
  save(plotList, file="plotList.Rdata", version=2)
  require(patchwork)
  stopifnot(length(plotList) == 4)
  x <- plotList[[1]] + plotList[[2]] + plotList[[3]] + plotList[[4]] + plot_layout(ncol=2)
  outFile <- file.path(outFolder, paste0(plotVar,"_pairedCmp_boxplot_allRandom.", plotType))
  ggsave(x, file=outFile, height=myHeight*2.5, width=myWidth*2.5)
  cat(paste0("... written: ", outFile, "\n"))
}


for(pat in all_patts) {
  
  tmp_dt <- random_dt[grepl(paste0(pat,"_40kb" ), random_dt$dataset),]
  stopifnot(nrow(tmp_dt) > 0)
  
  tmp_dt$dataset <- gsub("_RANDOM.+_40kb|_PERMUT.+_40kb", "_40kb", tmp_dt$dataset)
  stopifnot(setequal(tmp_dt$dataset, obs_dt$dataset))
  
  plot_dt <- merge(tmp_dt, obs_dt, by="dataset", suffixes=c("_rd", "_obs"), all=TRUE)
  # stopifnot(!is.na(plot_dt))  # not TRUE because I have some division by 0
  
  plot_dt$dotcols <- all_cols[all_cmps[basename(as.character(plot_dt$dataset))]]
  
  my_x <- plot_dt$nSignif_rd/plot_dt$nSignif_obs
  my_y <- plot_dt$nSignifAndFlagged_rd/plot_dt$nSignifAndFlagged_obs
  
  my_xlab <- "nSignif_rd/nSignif_obs"
  my_ylab <- "nSignifAndFlagged_rd/nSignifAndFlagged_obs"
  
  outFile <- file.path(outFolder, paste0("nSignif_rd_obs_vs_nSignifAndFlagged_rd_obs_", pat, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  plot(x=my_x,
       y=my_y,
       main=paste0(pat, "/ observed"),
       pch=16,
       col=plot_dt$dotcols,
       xlab=my_xlab,
       ylab=my_ylab)
  curve(1*x, add=T, col="darkgrey")
  # addCorr(x=my_x,y=my_y,bty="n", legPos="topleft")
  legend("bottomright", legend=names(all_cols), col=all_cols, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  my_x <- plot_dt$ratioSignif_rd/plot_dt$ratioSignif_obs
  my_y <- plot_dt$ratioSignifFlagged_rd/plot_dt$ratioSignifFlagged_obs
  
  my_xlab <- "ratioSignif_rd/ratioSignif_obs"
  my_ylab <- "ratioSignifFlagged_rd/nSignifAndFlagged_obs"
  
  outFile <- file.path(outFolder, paste0("ratioSignif_rd_obs_vs_ratioSignifFlagged_rd_obs_", pat, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  plot(x=my_x,
       y=my_y,
       main=paste0(pat, "/ observed"),
       pch=16,
       col=plot_dt$dotcols,
       xlab=my_xlab,
       ylab=my_ylab)
  curve(1*x, add=T, col="darkgrey")
  # addCorr(x=my_x,y=my_y,bty="n", legPos="topleft")
  legend("bottomright", legend=names(all_cols), col=all_cols, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  
}




stop("-ok")



















### output final table
setDir <- "/media/electron"
setDir <- ""
mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")
obs_dt <- get(load("SIGNIF_PURITYFLAGGED_FINAL/aran/CPE/log10/all_dt.Rdata"))
random_dt1 <- get(load("SIGNIF_PURITYFLAGGED_FINAL_RANDOMMIDPOS//aran/CPE/log10/all_dt.Rdata"))
random_dt1$data_type <- gsub(".+_(RANDOM.+)_40kb", "\\1", dirname(random_dt1$dataset))
random_dt1$dataset <- gsub("_RANDOM.+_40kb", "_40kb", random_dt1$dataset)
stopifnot(setequal(random_dt1$dataset, obs_dt$dataset))
random_dt2 <- get(load("SIGNIF_PURITYFLAGGED_FINAL_PERMUT//aran/CPE/log10/all_dt.Rdata"))
random_dt2$data_type <- gsub(".+_(PERMUT.+)_40kb", "\\1", dirname(random_dt2$dataset))
random_dt2$dataset <- gsub("_PERMUT.+_40kb", "_40kb", random_dt2$dataset)
stopifnot(setequal(random_dt2$dataset, obs_dt$dataset))
random_dt <- rbind(random_dt1, random_dt2)
obs_dt$ratioFlagged <- obs_dt$nPurityFlagged/obs_dt$nTot
random_dt$ratioFlagged <- random_dt$nPurityFlagged/random_dt$nTot
tmp_out_dt <- merge(obs_dt, random_dt, by="dataset", suffixes=c("_obs", "_rd"), all=TRUE)
all_ds <- unique(tmp_out_dt$dataset)
nSamp_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  cat(paste0("... start: ", ds, "\n"))
  hicds <- dirname(ds)
  exprds <- basename(ds)
  settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  samp1 <- get(load(file.path(setDir, sample1_file)))
  samp2 <- get(load(file.path(setDir, sample2_file)))
  data.frame(
    dataset = ds,
    cond1=cond1,
    cond2=cond2,
    nSamp1=length(samp1),
    nSamp2=length(samp2),
    stringsAsFactors = FALSE    
  )
}
out_dt <- merge(tmp_out_dt, nSamp_dt, by="dataset")

outFile <- file.path(outFolder, paste0("nSignif_purityTagged_obsRandom_withSamp.txt"))
write.table(out_dt, file=outFile, sep="\t", quote=F, append=F)
cat(paste0("... written: ", outFile, "\n"))

nSignif_purityTagged_obsRandom_withSamp <- out_dt
outFile <- file.path(outFolder, paste0("nSignif_purityTagged_obsRandom_withSamp.Rdata"))
save(nSignif_purityTagged_obsRandom_withSamp, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))



# all_result_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
# all_result_dt$dataset <- file.path(all_result_dt$hicds, all_result_dt$exprds)




# 
# 
# save(all_dt,file=file.path(outFolder, "all_dt.Rdata"), version=2)
# 
# all_dt$hicds_lab <- gsub(".+_(.+?)_40kb","\\1",  all_dt$hicds)
# all_dt$hicds_lab <- ifelse(grepl("RANDOM",all_dt$hicds_lab) | grepl("PERMUT", all_dt$hicds_lab),all_dt$hicds_lab,  
#                            "OBSERVED" )
# 
# all_dt$hicds_lab <- factor(all_dt$hicds_lab, 
#                            levels=c("OBSERVED", "RANDOMMIDPOS", "RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT", "PERMUTG2T"))
# 
# subTit <- paste0(paste0("# ", names(table(all_dt$hicds_lab)), "=", as.numeric(table(all_dt$hicds_lab))), collapse="; ")
# legText <- paste0("# ", names(table(all_dt$hicds_lab)), "=", as.numeric(table(all_dt$hicds_lab)))
# 
# 
# stopifnot(diff(table(all_dt$hicds_lab)) == 0)  # TO UNCOMMENT LATER!
# 
# outFile <- file.path(outFolder, paste0("nbrSignifTADs_densityplot.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# plot_multiDens(
#   split(all_dt$nSignif, all_dt$hicds_lab),
#   plotTit = paste0("# signif. TADs")
# )
# legend("topright", legend=legText, bty="n", cex=0.9)
# 
# mtext(side=3, text = paste0("adj. comb. p-val <= ", tad_signifThresh), font = 3)
# foo <- dev.off()
# 
# outFile <- file.path(outFolder, paste0("nbrSignifTADs_boxplot.", plotType))
# do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
# par(mar = par()$mar + c(10,3,0,0))
# boxplot(nSignif~hicds_lab, data = all_dt, outline=FALSE,
#         main = "# signif. TADs", xlab="", ylab="", 
#         cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
# 
# stripchart(nSignif~hicds_lab, vertical = TRUE, data = all_dt,
#            method = "jitter", add = TRUE, pch = 20, col = jitterCol)
# 
# legend("topright", legend=legText, bty="n", cex=0.9)
# 
# mtext(side=2, text="# signif. TADs", cex=plotCex, line=5)
# 
# mtext(side=3, text = paste0("adj. comb. p-val <= ", tad_signifThresh), font = 3)
# foo <- dev.off()
# 
# all_dt$ratioSignif <- all_dt$nSignif/all_dt$nTot
# outFile <- file.path(outFolder, paste0("ratioSignifTADs_boxplot.", plotType))
# do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
# par(mar = par()$mar + c(10,3,0,0))
# boxplot(ratioSignif~hicds_lab, outline=FALSE,
#         data = all_dt, main = "Ratio signif. TADs", 
#         xlab="", ylab="", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
# 
# stripchart(ratioSignif~hicds_lab, vertical = TRUE, data = all_dt,
#            method = "jitter", add = TRUE, pch = 20, col = jitterCol)
# legend("topright", legend=legText, bty="n", cex=0.9)
# 
# mtext(side=2, text="Ratio signif. TADs", cex=plotCex, line=5)
# mtext(side=3, text = paste0("adj. comb. p-val <= ", tad_signifThresh), font = 3)
# foo <- dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 





# 
# 
# 
# 
# 
# 
# 
# purity_file <- file.path("ALLTADS_AND_PURITY_FINAL_RANDOMMIDPOS", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata")  # here _final INPUT
# purityData <- get(load(purity_file))
# agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)
# 
# agg_purity$regID <- file.path(agg_purity$dataset, agg_purity$region)
# 
# result_file <- file.path("CREATE_FINAL_TABLE_RANDOM", "all_result_dt.Rdata")
# resultData <- get(load(result_file))
# resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
# 
# merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))
# merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
# purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))
# merge_dt$purityFlagged <- merge_dt$purityCorr <= purityCorrThresh
# merge_dt$signifFlagged <- merge_dt$signif & merge_dt$purityFlagged
# 
# aggSignif_merge_dt <- aggregate(signif~dataset, FUN=sum, data=merge_dt)
# colnames(aggSignif_merge_dt)[2] <- "nSignif"
# stopifnot(sum(aggSignif_merge_dt$nSignif) ==sum(merge_dt$signif))
# aggFlagged_merge_dt <- aggregate(purityFlagged~dataset, FUN=sum, data=merge_dt)
# colnames(aggFlagged_merge_dt)[2] <- "nPurityFlagged"
# stopifnot(sum(aggFlagged_merge_dt$nPurityFlagged) ==sum(merge_dt$purityFlagged))
# aggSignifFlagged_merge_dt <- aggregate(signifFlagged~dataset, FUN=sum, data=merge_dt)
# colnames(aggSignifFlagged_merge_dt)[2] <- "nSignifAndFlagged"
# stopifnot(sum(aggSignifFlagged_merge_dt$nSignifAndFlagged) ==sum(merge_dt$signif & merge_dt$purityFlagged ))
# 
# all_dt <- merge(merge(aggSignif_merge_dt, aggFlagged_merge_dt, by="dataset", all=TRUE ),aggSignifFlagged_merge_dt,by="dataset", all=TRUE)
# all_dt$ratioSignifFlagged <- all_dt$nSignifAndFlagged/all_dt$nSignif
# stopifnot(na.omit(all_dt)$ratioSignifFlagged >= 0 & na.omit(all_dt)$ratioSignifFlagged <= 1)
# all_dt <- all_dt[order(all_dt$ratioSignifFlagged, decreasing = TRUE),]      
# 
# 
# outFile <- file.path(outFolder, "all_dt.Rdata")
# save(all_dt, file=outFile, version=2)
# cat(paste0("... written: ", outFile, "\n"))
#     
# 
# all_dt$ratioSignifFlagged <- round(all_dt$ratioSignifFlagged,4)
# 
# outFile <- file.path(outFolder, "all_dt_signif_flagged.txt")
# write.table(all_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# all_patts <- c("RANDOMMIDPOS", "RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT")
# 
# for(pat in all_patts) {
# 
# tmp_dt <- all_dt[grepl(paste0(pat,"_40kb" ), all_dt$dataset),]
# stopifnot(nrow(tmp_dt) > 0)
# outFile <- file.path(outFolder, paste0(pat, "_all_dt_signif_flagged.txt"))
# write.table(tmp_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
# cat(paste0("... written: ", outFile, "\n"))
#  
# 
# }
# 
# 
# 
# resultData$regID <- file.path(resultData$hicds, resultData$exprds, resultData$region)
# flagged_signif_dt <- resultData[resultData$regID %in% merge_dt$regID[merge_dt$signifFlagged],c("regID", "region_genes")]
# flagged_signif_dt <- flagged_signif_dt[order(flagged_signif_dt$regID),]
# 
# flagged_signif_genes_dt <- do.call(rbind, apply(flagged_signif_dt, 1, function(x) data.frame(conserved_region=unique(x["regID"]), 
#                                                                                                                            symbol=unlist(strsplit(x["region_genes"], ",")),
#                                                                                                                            stringsAsFactors = FALSE)))
# rownames(flagged_signif_genes_dt) <- NULL
# 
# nFlagged_genes_dt <- data.frame(
#   symbol = names(table(flagged_signif_genes_dt$symbol)),
#   nSignifFlagged=as.numeric(table(flagged_signif_genes_dt$symbol)),
#   stringsAsFactors = FALSE
# )
# nFlagged_genes_dt <- nFlagged_genes_dt[order(nFlagged_genes_dt$nSignifFlagged, decreasing = TRUE),]
# 
# outFile <- file.path(outFolder, "nFlagged_genes.txt")
# write.table(nFlagged_genes_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F, append=F)
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
