##### !!! NEED TO ADD:
# -VARIANCE
# - CORR WITH PURITY


require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
require(ggsci)

# Rscript check_skmel.R

outFolder <- file.path("CHECK_SMEL")
dir.create(outFolder,recursive = TRUE)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

step0_folder <- "0_prepGeneData"
step3_folder <- "3_runMeanTADLogFC"
step4_folder <- "4_runMeanTADCorr"
step11_folder <- "11sameNbr_runEmpPvalCombined"
script8_name <- "8cOnlyFCC_runAllDown"

plotType <-  "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4

runFolder <- "."

pipOutFolder <- "PIPELINE/OUTPUT_FOLDER"

check_exprds <- "TCGAskcm_lowInf_highInf"
my_cols <- setNames(pal_jama()(5)[c(3, 2)], c(check_exprds, "other"))

buildData <- TRUE

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[! (grepl("PERMUTG2T", all_hicds) | grepl("RANDOM", all_hicds))]
hicds = all_hicds[1]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

if(buildData) {
  
  auc_ratio_file <- file.path("../MANUSCRIPT_FIGURES/FIG_1/FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata")
  all_auc_dt <- get(load(auc_ratio_file))


var_file <- "../OLDER_v2_Yuanlong_Cancer_HiC_data_TAD_DA/EXPR_VARIANCE_BYTAD/LOG2FPKM/all_ds_geneVarDT.Rdata"
var_values <- get(load(var_file))  

corr_file <- "ALLTADS_AND_PURITY_FINAL/aran/CPE/log10/all_ds_corrPurity_dt.Rdata"
corr_dt <- get(load(corr_file))  

  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
sub_corr_dt <- corr_dt[corr_dt$dataset == file.path(hicds, exprds),]

      stopifnot(paste0(hicds,"_",exprds) %in% names(var_values)) 
var_data <- var_values[[paste0(hicds, "_", exprds)]][["tadMeanVar"]]

      genelist_file <- file.path(pipOutFolder, hicds, exprds, 
                                 step0_folder, "pipeline_geneList.Rdata")
      stopifnot(file.exists(genelist_file))
      pipeline_genes <- get(load(genelist_file))
      
      regionlist_file <- file.path(pipOutFolder, hicds, exprds, 
                                   step0_folder, "pipeline_regionList.Rdata")
      stopifnot(file.exists(regionlist_file))
      pipeline_regions <- get(load(regionlist_file))
      

if(nrow(sub_corr_dt) == 0){
puritycorr <- NA
} else {
puritycorr <- setNames(sub_corr_dt$purityCorr, sub_corr_dt$region)
stopifnot(setequal(pipeline_regions, names(puritycorr)))
}



stopifnot(pipeline_regions %in% names(var_data))
stopifnot(setequal(pipeline_regions, names(var_data)))
      
      g2t_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
      g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
      g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
      stopifnot(pipeline_genes %in% g2t_dt$entrezID)
      stopifnot(pipeline_regions %in% g2t_dt$region)
      g2t_dt <- g2t_dt[g2t_dt$region %in% pipeline_regions & g2t_dt$entrezID %in% pipeline_genes,]
      region_ngenes <- setNames(as.numeric(table(g2t_dt$region)), as.character(names(table(g2t_dt$region))))
      
      stopifnot(setequal(pipeline_regions, names(region_ngenes)))
      
      combPval_file <- file.path(pipOutFolder, hicds, exprds, 
                                 step11_folder, "emp_pval_combined.Rdata")
      stopifnot(file.exists(combPval_file))
      pval <- get(load(combPval_file))
      adjPval <- p.adjust(pval, method="BH")
      
      stopifnot(setequal(pipeline_regions, names(adjPval)))
      
      logFC_file <- file.path(pipOutFolder, hicds, exprds, 
                              step3_folder, "all_meanLogFC_TAD.Rdata")
      stopifnot(file.exists(logFC_file))
      logFC <- get(load(logFC_file))
      
      stopifnot(setequal(pipeline_regions, names(logFC)))    
      
      meanCorr_file <- file.path(pipOutFolder, hicds, exprds, 
                                 step4_folder, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(meanCorr_file))
      meanCorr <- get(load(meanCorr_file))
      
      stopifnot(setequal(pipeline_regions, names(meanCorr)))    
      
      fcc_file <- file.path(pipOutFolder, hicds, exprds,
                            script8_name, "all_obs_prodSignedRatio.Rdata")
      fcc <- get(load(fcc_file))
      stopifnot(setequal(pipeline_regions, names(fcc)))    
      
      data.frame(
        hicds=hicds,
        exprds=exprds,
        region = pipeline_regions,
        nGenes = region_ngenes[pipeline_regions],
        adjPval = adjPval[pipeline_regions],
        logFC = logFC[pipeline_regions],
        meanCorr = meanCorr[pipeline_regions],
        FCC = fcc[pipeline_regions],
        meanVar = var_data[pipeline_regions],
        stringsAsFactors = FALSE
      )
    }
    exprds_dt
  }
  rownames(all_dt) <- NULL
  outFile <- file.path(outFolder, "all_dt.Rdata")
  save(all_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  # all_dt <- get(load("CHECK_SMEL/all_dt.Rdata"))
  outFile <- file.path(outFolder, "all_dt.Rdata")
  all_dt <- get(load(outFile))
}


legTitle <- "RNAseq data:"

fontFamily <- "Hershey"

all_vars <- colnames(all_dt)[!colnames(all_dt) %in% c("hicds", "exprds", "region")]
all_dt$exprds_lab <- ifelse(all_dt$exprds == check_exprds, check_exprds, "other")

plot_var=all_vars[1]
for(plot_var in all_vars) {
  
  plotTit <- paste0(plot_var)
  mySub <- paste0(gsub("_", " ", gsub("TCGA", "", check_exprds)), " vs. other")
  legTitle <- ""
  
  p3 <- ggdensity(all_dt,
                  x = paste0(plot_var),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0(plot_var),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "exprds_lab",
                  fill = "exprds_lab",
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = mySub)+
    scale_color_manual(values=my_cols)+
    scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    ) 
  
  outFile <- file.path(outFolder, paste0(plot_var, "_", check_exprds, "_vs_all_density.", plotType))
  ggsave(p3, file=outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
}

pvalThresh=0.01
nTADs_dt <- aggregate(region~hicds + exprds + exprds_lab, data = all_dt, FUN=length)
colnames(nTADs_dt)[4] <- "nbrTADs"
nSignif_dt <-  aggregate(adjPval~hicds + exprds + exprds_lab, data = all_dt, 
                         FUN=function(x) sum(x<=pvalThresh))
colnames(nSignif_dt)[4] <- "nbrSignif"
meanGenes_dt <-  aggregate(nGenes~hicds + exprds + exprds_lab, data = all_dt, FUN=mean)
colnames(meanGenes_dt)[4] <- "meanNbrGenes"
totGenes_dt <-  aggregate(nGenes~hicds + exprds + exprds_lab, data = all_dt, FUN=sum)
colnames(totGenes_dt)[4] <- "totNbrGenes"

nbrGenes3_dt <-  aggregate(nGenes~hicds + exprds + exprds_lab, data = all_dt, FUN=function(x) sum(x==3))
colnames(nbrGenes3_dt)[4] <- "nbrGenes3"

nbrFCC1_dt <-  aggregate(FCC~hicds + exprds + exprds_lab, data = all_dt, FUN=function(x)sum(round(x,4) == 1))
colnames(nbrFCC1_dt)[4] <- "nbrFCC1"


all_dt2 <- merge(nbrFCC1_dt, merge(totGenes_dt, merge(meanGenes_dt, merge(nSignif_dt, merge(nTADs_dt, nbrGenes3_dt, 
                                  by=c("hicds", "exprds", "exprds_lab"), all=T),
                                   by=c("hicds", "exprds", "exprds_lab"), all=T),
                                   by=c("hicds", "exprds", "exprds_lab"), all=T),
                                   by=c("hicds", "exprds", "exprds_lab"), all=T),
                 by=c("hicds", "exprds", "exprds_lab"), all=T)



all_vars <- colnames(all_dt2)[!colnames(all_dt2) %in% c("hicds", "exprds", "exprds_lab")]

plot_var=all_vars[1]
for(plot_var in all_vars) {
  
  plotTit <- paste0(plot_var)
  mySub <- paste0(gsub("_", " ", gsub("TCGA", "", check_exprds)), " vs. other")
  legTitle <- ""
  
  all_dt2 <- all_dt2[order(all_dt2[,plot_var], decreasing = TRUE),]
  all_dt2$dataset <- file.path(as.character(all_dt2$hicds), as.character(all_dt2$exprds))
  all_dt2$dataset <- factor(as.character(all_dt2$dataset), levels=as.character(all_dt2$dataset))
  
  p3 <- ggbarplot(all_dt2,
                  x = "dataset",
                  y = paste0(plot_var),
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab ="datasets",
                  ylab = paste0(plot_var),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "exprds_lab",
                  fill = "exprds_lab",
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = mySub)+
    scale_color_manual(values=my_cols)+
    scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle)) +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    ) 
  
  outFile <- file.path(outFolder, paste0(plot_var, "_", check_exprds, "_vs_all_barplot.", plotType))
  ggsave(p3, file=outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
}

                                   
