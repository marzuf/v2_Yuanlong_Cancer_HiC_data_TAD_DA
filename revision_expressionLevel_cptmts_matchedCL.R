require(ggsci)
require(ggpubr)
require(ggplot2)
require(ggrepel)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("revision_settings.R")
require(stringr)

plotType <- "png"
myWidthGG <- 7
myHeightGG <- 5
myHeight <- 400
myWidth <- 400

plotCex <- 1.2

tadSignifThresh <- 0.01

geneWithExprThresh <- 3

# Rscript revision_expressionLevel_cptmts_matchedCL.R

tad_aggFun <- "mean"
plotCols <- c(paste0(tad_aggFun, "_geneExpr"), paste0(tad_aggFun, "_geneExprLog10"))

outFolder <- file.path("REVISION_EXPRESSIONLEVEL_CPTMTS_MATCHEDCL")
dir.create(outFolder, recursive=TRUE)

inFolder <- "REVISION_EXPRESSIONDATA_ALL_MATCHEDCL"
outFile <- file.path(inFolder, "all_matched_mRNA_dt.Rdata")
all_matched_mRNA_dt <- get(load(outFile))
nrow(all_matched_mRNA_dt)
# 215257
all_matched_mRNA_dt$hicds_region <- file.path(all_matched_mRNA_dt$hicds, all_matched_mRNA_dt$region)
notfiltered_dt <- all_matched_mRNA_dt
all_matched_mRNA_dt <- all_matched_mRNA_dt[all_matched_mRNA_dt$nGenes_withExpr >= geneWithExprThresh,]
nrow(all_matched_mRNA_dt)
# 104427

inFolder <- "REVISION_CHANGES_CPTMTLABELS_ALLDS_ALLTADS"
outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
tad2cptmt_dt <- get(load(outFile))

stopifnot(tad2cptmt_dt$hicds_region %in% notfiltered_dt$hicds_region)

stopifnot(tad2cptmt_dt$tad_binaryCptmtLab == tad2cptmt_dt$start_binaryCptmtLab)
stopifnot(tad2cptmt_dt$tad_eightCptmtLab == tad2cptmt_dt$start_eightCptmtLab)

tad2cptmt_dt <- tad2cptmt_dt[,c("hicds", "region", "tadCptmtNormRank", "tad_binaryCptmtLab", "tad_eightCptmtLab" )]
tad2cptmt_dt <- unique(tad2cptmt_dt)
tad2cptmt_dt$hicds_region <- file.path(tad2cptmt_dt$hicds, tad2cptmt_dt$region )
stopifnot(!duplicated(tad2cptmt_dt$hicds_region))

stopifnot(all_matched_mRNA_dt$hicds_region %in% tad2cptmt_dt$hicds_region)

meanExprValues <- setNames(all_matched_mRNA_dt$mean_geneExpr, all_matched_mRNA_dt$hicds_region)
# stopifnot(tad2cptmt_dt$hicds_region %in% names(meanExprValues)) # not true because I remove those that have too few genes !!!
tad2cptmt_dt$mean_geneExpr <- meanExprValues[paste0(tad2cptmt_dt$hicds_region)]
  
meanExprValues_log10 <- setNames(all_matched_mRNA_dt$mean_geneExprLog10, all_matched_mRNA_dt$hicds_region)
# stopifnot(tad2cptmt_dt$hicds_region %in% names(meanExprValues_log10)) # not true because I remove those that have too few genes !!!
tad2cptmt_dt$mean_geneExprLog10 <- meanExprValues_log10[paste0(tad2cptmt_dt$hicds_region)]





legTitle <- ""

tad2cptmt_dt$tad_binaryCptmtLab <- factor(tad2cptmt_dt$tad_binaryCptmtLab,
                                          levels=as.character(sort(unique(as.character(
                                            tad2cptmt_dt$tad_binaryCptmtLab)))))
stopifnot(!is.na(tad2cptmt_dt$tad_binaryCptmtLab))


tad2cptmt_dt$tad_eightCptmtLab <- factor(tad2cptmt_dt$tad_eightCptmtLab,
                                          levels=as.character(sort(unique(as.character(
                                            tad2cptmt_dt$tad_eightCptmtLab)))))
stopifnot(!is.na(tad2cptmt_dt$tad_eightCptmtLab))



##################3 START PLOTTING

plotTit <- ""
mysub <- ""


all_cptmt_vars <- c("tad_binaryCptmtLab","tad_eightCptmtLab")

all_plot_vars <- c("mean_geneExpr", "mean_geneExprLog10")

for(cptmt_var in all_cptmt_vars){
  
  for(plot_var in all_plot_vars) {
    
    plotTit <- paste0(plot_var, " by ", cptmt_var)
    
    mySub <- paste0("# DS = ", length(unique(file.path(tad2cptmt_dt$hicds,tad2cptmt_dt$exprds))),
                                 "; # TADs = ", nrow(tad2cptmt_dt))
    
    pbox <- ggboxplot(tad2cptmt_dt, 
                      x=paste0(cptmt_var),
                      y=paste0(plot_var),
                      add="jitter") + 
      mytheme +
      ggtitle(plotTit, subtitle = mySub)+
      # scale_color_manual(values=my_cols)+
      # scale_fill_manual(values=my_cols)  +
      labs(color=paste0(legTitle),fill=paste0(legTitle), x="", y=paste0("TAD ", plot_var))+
      # guides(color=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
    
    outFile <- file.path(outFolder,paste0(plot_var, "_byCptmt_", cptmt_var, "_boxplot.", plotType))
    ggsave(pbox, filename = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}

############## DENSPLOT WITH THE CONTINUOUS RANKS

all_cptmt_vars <- c("tadCptmtNormRank")
all_plot_vars <- c("mean_geneExpr", "mean_geneExprLog10")

cptmt_var="tadCptmtNormRank"
plot_var="meanCorr"

for(cptmt_var in all_cptmt_vars){
  
  for(plot_var in all_plot_vars) {
    
    
    myx <- tad2cptmt_dt[,paste0(cptmt_var)]
    myy <- tad2cptmt_dt[,paste0(plot_var)]

        
    plotTit <- paste0( plot_var, " by ", cptmt_var)
    
    mySub <- paste0("# DS = ", length(unique(file.path(tad2cptmt_dt$hicds,tad2cptmt_dt$exprds))),
                    "; # TADs = ", nrow(tad2cptmt_dt))
    
    outFile <- file.path(outFolder,paste0(plot_var, "_byCptmt_", cptmt_var, "_densplot.", plotType))
    do.call(plotType, list(outFile, height=myWidth, width=myWidth))
    densplot(x=myx,
             y=myy,
             xlab = paste0(cptmt_var),
             ylab = paste0(plot_var),
             main  = plotTit,
             plot.main=plotCex,
             plot.main=plotCex,
             plot.main=plotCex
    )
    mtext(side=3, text=mySub)
    addCorr(x=myx,y=myy, legPos="topleft", bty="n", corMet="spearman")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    

  }
}


##### pie chart signif by cptmts
