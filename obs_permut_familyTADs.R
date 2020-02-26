require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

startTime <- Sys.time()

script_name <- "obs_permut_familyTADs.R"

# Rscript obs_permut_familyTADs.R

plotType <- "svg"
myHeight <- myWidth <- 7
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "OBS_PERMUT_FAMILY_TADS"
dir.create(outFolder, recursive = TRUE)


final_DT <- get(load(file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")))
final_DT_permut <- get(load(file.path("CREATE_FINAL_TABLE_RANDOM", "all_result_dt.Rdata")))

result_dt_all <- rbind(final_DT, final_DT_permut)

exprds <- "TCGAluad_norm_luad"
famType <- "hgnc_family_short"

tad_signif_thresh <- 0.05

myhicds <- "ENCSR489OCU_NCI-H460"

plotCex <- 1.4

buildData <- TRUE

all_hicds <- c(
  "ENCSR489OCU_NCI-H460_40kb",
  "ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb",
  "ENCSR489OCU_NCI-H460_RANDOMNBRGENES_40kb",
  "ENCSR489OCU_NCI-H460_RANDOMSHIFT_40kb",
  "ENCSR489OCU_NCI-H460_PERMUTG2T_40kb"
)

if(buildData) {
  all_data <- foreach(hicds = all_hicds) %dopar% {
    
    cat("... start ", hicds, "\n")
    
    fam_dt <- get(load(file.path("PREP_GENE_FAMILIES_TAD_DATA", hicds, "hgnc_entrezID_family_TAD_DT.Rdata")))
    fam_dt$entrezID <- as.character(fam_dt$entrezID)
    fam_dt <- fam_dt[order(fam_dt$entrezID),]
    stopifnot(!duplicated(fam_dt$entrezID))
    entrez2fam <- setNames(fam_dt[,paste0(famType)], fam_dt$entrezID)
    entrez2reg_fam <- setNames(fam_dt$region, fam_dt$entrezID)
    
    
    geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
    sum(geneList %in% fam_dt$entrezID)/length(geneList)
    
    g2t_dt_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
    stopifnot(geneList %in% g2t_dt$entrezID)
    g2t_pip_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
    
    fam_pip_dt <- fam_dt[fam_dt$entrezID %in% geneList,]
    
    tmp <- g2t_dt[g2t_dt$entrezID %in% fam_dt$entrezID,]
    tmp <- tmp[order(tmp$entrezID),]
    entrez2reg_tmp <- setNames(tmp$region, tmp$entrezID)
    stopifnot(entrez2reg_tmp == entrez2reg_fam)
    stopifnot(names(entrez2reg_tmp) == names(entrez2reg_fam))
    
    ###################################################################### AGGREG BY TAD
    # FOR THE TADs: # of families; # of families/# genes; meanCorr and pval
    ###################################################################### 
    
    nFamByTAD_dt <- aggregate(as.formula(paste0(famType, "~region")), data=fam_pip_dt, FUN=function(x) length(unique(x)))
    colnames(nFamByTAD_dt)[colnames(nFamByTAD_dt) == paste0(famType)] <- "nFams"
    nGenesByTAD_dt <- aggregate(entrezID ~  region, data=g2t_pip_dt, FUN=length)
    colnames(nGenesByTAD_dt)[colnames(nGenesByTAD_dt) == "entrezID"] <- "nGenes"
    # OR TO TAKE ONLY ANNOTATED GENES:
    # nFamByTAD_dt <- aggregate(as.formula(paste0(famType, "~region")), data=fam_dt, FUN=function(x) length(unique(x)))
    
    nGenesFamsByTAD_dt <- merge(nFamByTAD_dt, nGenesByTAD_dt, all.x=TRUE, all.y=FALSE, by="region")
    stopifnot(!is.na(nGenesFamsByTAD_dt$nGenes))
    
    
    result_dt <- result_dt_all[result_dt_all$hicds == hicds & result_dt_all$exprds == exprds,]
    stopifnot(!duplicated(result_dt$region))
    result_dt$region_rank <- rank(result_dt$adjPvalComb, ties="min")
    
    
    stopifnot(nGenesFamsByTAD_dt$region %in% result_dt$region)
    
    nGenesFamsByTAD_final_dt <- merge(nGenesFamsByTAD_dt, result_dt[,c("hicds", "exprds", "region", "meanCorr", "adjPvalComb", "region_rank")],
                                      all.x=TRUE, all.y=FALSE, by="region")
    stopifnot(!is.na(nGenesFamsByTAD_final_dt))
    stopifnot(nGenesFamsByTAD_final_dt$hicds == hicds)
    stopifnot(nGenesFamsByTAD_final_dt$exprds == exprds)
    
    ###################################################################### AGGREG BY FAM
    # FOR THE FAMs:  #  of TADs
    ###################################################################### 
    
    nTADsByFam_dt <- aggregate(as.formula(paste0("region ~ ", famType)), data=fam_pip_dt, FUN=function(x) length(unique(x)))
    colnames(nTADsByFam_dt)[colnames(nTADsByFam_dt) == "region"] <- "nTADs"
    stopifnot(fam_pip_dt$region %in% result_dt$region)
    
    fam_pip_result_dt <- merge(fam_pip_dt, result_dt[,c("hicds", "exprds", "region", "adjPvalComb", "region_rank")], all.x=TRUE, all.y=FALSE, by="region")
    stopifnot(!is.na(fam_pip_result_dt))
    stopifnot(fam_pip_result_dt$hicds == hicds)
    stopifnot(fam_pip_result_dt$exprds == exprds)
    
    tad_fam_pip_result_dt <- fam_pip_result_dt[,c("region", "adjPvalComb", "region_rank", paste0(famType))]
    tad_fam_pip_result_dt <- unique(tad_fam_pip_result_dt)
    
    nSignifTADsByFam_dt <- aggregate(as.formula(paste0("adjPvalComb ~ ", famType)), data=tad_fam_pip_result_dt, FUN=function(x) sum(x <= tad_signif_thresh))
    colnames(nSignifTADsByFam_dt)[ colnames(nSignifTADsByFam_dt) == "adjPvalComb"] <- "nSignifTADs"
    nNotSignifTADsByFam_dt <- aggregate(as.formula(paste0("adjPvalComb ~ ", famType)), data=tad_fam_pip_result_dt, FUN=function(x) sum(x > tad_signif_thresh))
    colnames(nNotSignifTADsByFam_dt)[ colnames(nNotSignifTADsByFam_dt) == "adjPvalComb"] <- "nNotSignifTADs"
    
    all_nTADsByFam_dt <- merge(nTADsByFam_dt, merge(  nSignifTADsByFam_dt, nNotSignifTADsByFam_dt, all=TRUE, by=paste0(famType)), 
                               all=TRUE, by=paste0(famType))
    
    stopifnot(!is.na(all_nTADsByFam_dt))
    stopifnot(all_nTADsByFam_dt$nTADs == all_nTADsByFam_dt$nSignifTADs + all_nTADsByFam_dt$nNotSignifTADs)
    
    all_nTADsByFam_dt$hicds <- hicds
    all_nTADsByFam_dt$exprds <- exprds
    
    list(
      nGenesFamsByTAD_final_dt = nGenesFamsByTAD_final_dt,
      all_nTADsByFam_dt=all_nTADsByFam_dt
    )
    
  }
  
  outFile <- file.path(outFolder, "all_data.Rdata")
  save(all_data, file =outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {
  outFile <- file.path(outFolder, "all_data.Rdata")
  all_data <- get(load(outFile))
}

allDS_nTADsByFam_dt <- do.call(rbind,
                               lapply(all_data, function(x)x[[paste0("all_nTADsByFam_dt")]]))

allDS_nGenesFamsByTAD_dt <- do.call(rbind,
                                    lapply(all_data, function(x)x[[paste0("nGenesFamsByTAD_final_dt")]]))

### => BY FAM
allDS_nTADsByFam_dt$ratioSignifTADs <- allDS_nTADsByFam_dt$nSignifTADs/allDS_nTADsByFam_dt$nTADs



allDS_nTADsByFam_dt$hicds_lab <- gsub(paste0(myhicds, "_(.+)"),"\\1",  allDS_nTADsByFam_dt$hicds)

for(plot_var in c("ratioSignifTADs", "nTADs")) {
  
  if(plot_var == "ratioSignifTADs"){
    mysub <- paste0(myhicds, " -", exprds, "; TAD signif. thresh <= ", tad_signif_thresh)
  } else {
    mysub <- paste0(myhicds, " -", exprds)
  }
  
  
  
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byFam_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(split(allDS_nTADsByFam_dt[,paste(plot_var)], allDS_nTADsByFam_dt$hicds_lab),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  if(plot_var != "ratioSignifTADs"){
  
    outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byFam_density_log10.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot_multiDens(split(log10(allDS_nTADsByFam_dt[,paste(plot_var)]), allDS_nTADsByFam_dt$hicds_lab),
                   plotTit=paste0(plot_var, " [log10]"), legPos = "topright")
    mtext(side=3, text = paste0(mysub))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  
  }
  
  p_box <- ggboxplot(data=allDS_nTADsByFam_dt, x="hicds_lab", y=paste0(plot_var))+
    ggtitle(plot_var, subtitle=mysub)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    theme( # Increase size of axis lines
      strip.text = element_text(size = 12),
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"))
  
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byFam_boxplot.", plotType))
  ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byFam_density_noG2t.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(split(allDS_nTADsByFam_dt[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds),paste(plot_var)], allDS_nTADsByFam_dt$hicds_lab[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds)]),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  if(plot_var != "ratioSignifTADs"){
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byFam_density_noG2t_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(split(
    log10(allDS_nTADsByFam_dt[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds),paste(plot_var)]),
    allDS_nTADsByFam_dt$hicds_lab[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds)]),
                 plotTit=paste0(plot_var, " [log10]"), legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  }
  
  p_box <- ggboxplot(data=allDS_nTADsByFam_dt[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds),], x="hicds_lab", y=paste0(plot_var))+
    ggtitle(plot_var, subtitle=mysub)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    theme( # Increase size of axis lines
      strip.text = element_text(size = 12),
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"))
  
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byFam_boxplot_noG2t.", plotType))
  ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  
}

### => BY TAD

allDS_nGenesFamsByTAD_dt$adjPvalComb_log10 <- -log10(allDS_nGenesFamsByTAD_dt$adjPvalComb )
allDS_nGenesFamsByTAD_dt$hicds_lab <- gsub("ENCSR489OCU_NCI-H460_(.+)","\\1",  allDS_nGenesFamsByTAD_dt$hicds)


for(plot_var in c("adjPvalComb_log10", "meanCorr")) {
  for(hicds in all_hicds){
    plot_dt <- allDS_nGenesFamsByTAD_dt[allDS_nGenesFamsByTAD_dt$hicds == hicds,]
    
    
      
      my_x <- plot_dt[,c("nFams")]
      my_y <- plot_dt[,c(plot_var)]
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", plot_var, "_vs_nFams_byTAD_densplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      
      
      densplot(
        x = my_x,
        y=my_y,
        xlab = "# fam. in TAD",
        ylab = paste0("TAD ", plot_var),  
        cex=0.7,
        main = paste0(hicds),
        cex.lab=plotCex,
        cex.axis = plotCex,
        cex.main = plotCex
      )  
    mtext(side = 3, text = exprds)
    addCorr(x=my_x, y=my_y, bty="n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    }
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_allTADs_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(split(allDS_nGenesFamsByTAD_dt[,paste0(plot_var)], allDS_nGenesFamsByTAD_dt$hicds_lab),
                 plotTit=plot_var, legPos = "topleft")
  mtext(side=3, text = paste0(myhicds, " -", exprds))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_allTADs_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(allDS_nGenesFamsByTAD_dt[!grepl("_PERMUTG2T",allDS_nGenesFamsByTAD_dt$hicds ),paste0(plot_var)], allDS_nGenesFamsByTAD_dt[!grepl("_PERMUTG2T",allDS_nGenesFamsByTAD_dt$hicds ),]$hicds_lab),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = paste0(myhicds, " -", exprds))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_allTADs_boxplot.", plotType))
  
  p_box <- ggboxplot(data=allDS_nGenesFamsByTAD_dt, x="hicds_lab", y=paste0(plot_var)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    theme( # Increase size of axis lines
      strip.text = element_text(size = 12),
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"))
  
  ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
}


for(plot_var in c("nFams")) {
  
  if(plot_var == "ratioSignifTADs"){
    mysub <- paste0(myhicds, " -", exprds, "; TAD signif. thresh <= ", tad_signif_thresh)
  } else {
    mysub <- paste0(myhicds, " -", exprds)
  }
  
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byTAD_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(split(allDS_nGenesFamsByTAD_dt[,paste(plot_var)], allDS_nGenesFamsByTAD_dt$hicds_lab),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  

  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byTAD_density_noG2t.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(split(allDS_nGenesFamsByTAD_dt[!grepl("PERMUTG2T", allDS_nGenesFamsByTAD_dt$hicds),paste(plot_var)], allDS_nGenesFamsByTAD_dt$hicds_lab[!grepl("PERMUTG2T", allDS_nGenesFamsByTAD_dt$hicds)]),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))

  
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byTAD_density_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(split(
    log10(allDS_nGenesFamsByTAD_dt[,paste(plot_var)]), allDS_nGenesFamsByTAD_dt$hicds_lab),
                 plotTit=paste0(plot_var, " [log10]"), legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_", plot_var, "_byTAD_density_noG2t_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(split(log10(allDS_nGenesFamsByTAD_dt[!grepl("PERMUTG2T", allDS_nGenesFamsByTAD_dt$hicds),paste(plot_var)]), allDS_nGenesFamsByTAD_dt$hicds_lab[!grepl("PERMUTG2T", allDS_nGenesFamsByTAD_dt$hicds)]),
                 plotTit=paste0(plot_var, " [log10]"), legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}

#############################################################################################################################

# nFAMs of topTADs

nTopRank <- 50

plot_dt <- allDS_nGenesFamsByTAD_dt
plot_dt <- plot_dt[plot_dt$region_rank <= nTopRank,]

sub <- paste0(nTopRank, " top-ranking TADs only")

outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nFams_", nTopRank, "nTopTADs_boxplot.", plotType))

p_box <- ggboxplot(data=plot_dt, x="hicds_lab", y=paste0("nFams")) +
  ggtitle("# Fams by TAD", subtitle = paste0(mysub, "; ", sub))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"))

ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0(myhicds, "_", exprds, "_nFams_", nTopRank, "nTopTADs_density.", plotType))

do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(plot_dt$nFams, plot_dt$hicds_lab),
               plotTit=paste0("# Fams by TAD"), legPos = "topright")
mtext(side=3, text = paste0(mysub, "; ", sub))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))












