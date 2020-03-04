require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

startTime <- Sys.time()

script_name <- "obs_permut_familyTADs_allDS_min3.R"

# Rscript obs_permut_familyTADs_allDS_min3.R

plotType <- "png"
myHeight <- myWidth <- 400
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "OBS_PERMUT_FAMILY_TADS_ALLDS_MIN3"
dir.create(outFolder, recursive = TRUE)

rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T" , "RANDOMMIDPOSDISC")

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))


famType <- "hgnc_family_short"

plotCex <- 1.4

buildData <- TRUE


minFamGenes <- 3

if(buildData) {
  all_data <- foreach(hicds = all_hicds) %dopar% {
    
    cat("... start ", hicds, "\n")
    
    fam_dt <- get(load(file.path("PREP_GENE_FAMILIES_TAD_DATA", hicds, "hgnc_entrezID_family_TAD_DT.Rdata")))
    fam_dt$entrezID <- as.character(fam_dt$entrezID)
    fam_dt <- fam_dt[order(fam_dt$entrezID),]
    
    stopifnot(!duplicated(fam_dt$entrezID))
    nByFams <- setNames(as.numeric(table(fam_dt[, paste0(famType)])), names(table(fam_dt[, paste0(famType)])))
    
    keepFams <- names(nByFams)[nByFams >= minFamGenes ]
    stopifnot(keepFams %in% fam_dt[,paste0(famType)])
    
    fam_dt <-  fam_dt[fam_dt[,paste0(famType)] %in% keepFams,]
    stopifnot(setequal(keepFams, fam_dt[,paste0(famType)]))
    
    stopifnot(!duplicated(fam_dt$entrezID))
    entrez2fam <- setNames(fam_dt[,paste0(famType)], fam_dt$entrezID)
    entrez2reg_fam <- setNames(fam_dt$region, fam_dt$entrezID)
    
    
    
    exprds_data <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      mean_corr <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata")))
      
      meanCorr_dt <- data.frame(
        region = names(mean_corr),
        meanCorr = as.numeric(mean_corr),
        stringsAsFactors = FALSE
      )
    
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
      
      # save(meanCorr_dt, file="meanCorr_dt.Rdata", version=2)
      # save(nGenesFamsByTAD_dt, file="nGenesFamsByTAD_dt.Rdata", version=2)
      
      nGenesFamsCorrByTAD_dt <- merge(meanCorr_dt, nGenesFamsByTAD_dt, by ="region", all.x=FALSE, all.y=FALSE ) # not all X -> where nFam=0, not retained
      # stopifnot(!is.na(nGenesFamsCorrByTAD_dt))
      
      nGenesFamsCorrByTAD_dt$hicds <- hicds
      nGenesFamsCorrByTAD_dt$exprds <- exprds
      
  
      ###################################################################### AGGREG BY FAM
      # FOR THE FAMs:  #  of TADs
      ###################################################################### 
      
      nTADsByFam_dt <- aggregate(as.formula(paste0("region ~ ", famType)), data=fam_pip_dt, FUN=function(x) length(unique(x)))
      colnames(nTADsByFam_dt)[colnames(nTADsByFam_dt) == "region"] <- "nTADs"
      
      nTADsByFam_dt$hicds <- hicds
      nTADsByFam_dt$exprds <- exprds
      
      
      
      
      list(
        nGenesFamsCorrByTAD_dt = nGenesFamsCorrByTAD_dt,
        nTADsByFam_dt=nTADsByFam_dt
      )
    } 
    names(exprds_data) <- all_exprds[[paste0(hicds)]]
    exprds_data
  }
  names(all_data) <- all_hicds
  
  outFile <- file.path(outFolder, "all_data.Rdata")
  save(all_data, file =outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {
  outFile <- file.path(outFolder, "all_data.Rdata")
  all_data <- get(load(outFile))
}


allDS_nTADsByFam_dt <- do.call(rbind,
                               lapply(all_data, function(subdata) do.call(rbind, lapply(subdata, function(x)x[[paste0("nTADsByFam_dt")]]))))

allDS_nGenesFamsByTAD_dt <- do.call(rbind,
                               lapply(all_data, function(subdata) do.call(rbind, lapply(subdata, function(x)x[[paste0("nGenesFamsCorrByTAD_dt")]]))))


### => BY FAM



allDS_nTADsByFam_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                                           gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                                gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                                     gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC", 
                                                     gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", allDS_nTADsByFam_dt$hicds)))))
allDS_nTADsByFam_dt$hicds_lab[! allDS_nTADsByFam_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"



mysub <- paste0( "all DS - n=", length(unique(file.path(allDS_nTADsByFam_dt$hicds, allDS_nTADsByFam_dt$exprds))))


  for(plot_var in c("nTADs")) {

  
  outFile <- file.path(outFolder, paste0( "allDS_", plot_var, "_byFam_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(allDS_nTADsByFam_dt[,paste(plot_var)], allDS_nTADsByFam_dt$hicds_lab),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = mysub)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0( "allDS_", plot_var, "_byFam_density_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(log10(allDS_nTADsByFam_dt[,paste(plot_var)]), allDS_nTADsByFam_dt$hicds_lab),
                 plotTit=paste0(plot_var, " [log10]"), legPos = "topright")
  mtext(side=3, text = mysub)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  p_box <- ggboxplot(data=allDS_nTADsByFam_dt, x="hicds_lab", y=paste0(plot_var))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    ggtitle(plot_var, subtitle = mysub)+
    theme( # Increase size of axis lines
      strip.text = element_text(size = 12),
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"))
  
  
  outFile <- file.path(outFolder, paste0("allDS_", plot_var, "_byFam_boxplot.", plotType))
  ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0( "allDS_", plot_var, "_byFam_density_noG2t.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(allDS_nTADsByFam_dt[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds),paste(plot_var)], allDS_nTADsByFam_dt$hicds_lab[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds)]),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0( "allDS_", plot_var, "_byFam_density_noG2t_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(log10(allDS_nTADsByFam_dt[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds),paste(plot_var)]), allDS_nTADsByFam_dt$hicds_lab[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds)]),
                 plotTit=paste0(plot_var, " [log10]"), legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  p_box <- ggboxplot(data=allDS_nTADsByFam_dt[!grepl("PERMUTG2T", allDS_nTADsByFam_dt$hicds),], x="hicds_lab", y=paste0(plot_var))+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    ggtitle(plot_var, subtitle = mysub)+
    theme( # Increase size of axis lines
      strip.text = element_text(size = 12),
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"))
  
  
  outFile <- file.path(outFolder, paste0( "allDS_", plot_var, "_byFam_boxplot_noG2t.", plotType))
  ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  
}

### => BY TAD


allDS_nGenesFamsByTAD_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                          gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                               gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                    gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC", 
                                    gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", allDS_nGenesFamsByTAD_dt$hicds)))))
allDS_nGenesFamsByTAD_dt$hicds_lab[! allDS_nGenesFamsByTAD_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"

all_hicds_lab <- unique(allDS_nGenesFamsByTAD_dt$hicds_lab)


for(plot_var in c( "meanCorr")) {
  for(hicds in all_hicds_lab){
    plot_dt <- allDS_nGenesFamsByTAD_dt[allDS_nGenesFamsByTAD_dt$hicds_lab == hicds,]
    
      
      my_x <- plot_dt[,c("nFams")]
      my_y <- plot_dt[,c(plot_var)]
      
      outFile <- file.path(outFolder, paste0(hicds, "_", "allDS_", plot_var, "_vs_nFams_byTAD_densplot.", plotType))
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
    mtext(side = 3, text = mysub)
    addCorr(x=my_x, y=my_y, bty="n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    }
  
  outFile <- file.path(outFolder, paste0(hicds, "_allDS_", plot_var, "_allTADs_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(allDS_nGenesFamsByTAD_dt[,paste0(plot_var)], allDS_nGenesFamsByTAD_dt$hicds_lab),
                 plotTit=plot_var, legPos = "topleft")
  mtext(side=3, text = paste0(hicds, " - all DS"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0(hicds, "_allDS_", plot_var, "_allTADs_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(allDS_nGenesFamsByTAD_dt[!grepl("_PERMUTG2T",allDS_nGenesFamsByTAD_dt$hicds ),paste0(plot_var)],
                       allDS_nGenesFamsByTAD_dt[!grepl("_PERMUTG2T",allDS_nGenesFamsByTAD_dt$hicds ),]$hicds_lab),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = paste0(hicds, " - all DS"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0(hicds, "_allDS_", plot_var, "_allTADs_boxplot.", plotType))
  
  p_box <- ggboxplot(data=allDS_nGenesFamsByTAD_dt, x="hicds_lab", y=paste0(plot_var)) +
    ggtitle(plot_var, subtitle=paste0(hicds, " - all DS"))+
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
  
  
  outFile <- file.path(outFolder, paste0( "allDS_", plot_var, "_byTAD_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(allDS_nGenesFamsByTAD_dt[,paste(plot_var)], allDS_nGenesFamsByTAD_dt$hicds_lab),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = mysub)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  

  outFile <- file.path(outFolder, paste0( "allDS_", plot_var, "_byTAD_density_noG2t.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(allDS_nGenesFamsByTAD_dt[!grepl("PERMUTG2T", allDS_nGenesFamsByTAD_dt$hicds),paste(plot_var)], allDS_nGenesFamsByTAD_dt$hicds_lab[!grepl("PERMUTG2T", allDS_nGenesFamsByTAD_dt$hicds)]),
                 plotTit=plot_var, legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0( "allDS_", plot_var, "_byTAD_density_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(log10(allDS_nGenesFamsByTAD_dt[,paste(plot_var)]), allDS_nGenesFamsByTAD_dt$hicds_lab),
                 plotTit=paste0(plot_var, " [log10]"), legPos = "topright")
  mtext(side=3, text = mysub)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0( "allDS_", plot_var, "_byTAD_density_noG2t_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
  plot_multiDens(split(log10(allDS_nGenesFamsByTAD_dt[!grepl("PERMUTG2T", allDS_nGenesFamsByTAD_dt$hicds),paste(plot_var)]), allDS_nGenesFamsByTAD_dt$hicds_lab[!grepl("PERMUTG2T", allDS_nGenesFamsByTAD_dt$hicds)]),
                 plotTit=paste0(plot_var, " [log10]"), legPos = "topright")
  mtext(side=3, text = paste0(mysub))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
}



#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))












