require(ggsci)
require(ggpubr)
require(ggplot2)
require(ggrepel)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("revision_settings.R")

plotType <- "png"
myWidthGG <- 7
myHeightGG <- 5
myHeight <- 400
myWidth <- 400

plotCex <- 1.2

tadSignifThresh <- 0.01

# Rscript revision_purityFilter_cptmts.R

outFolder <- file.path("REVISION_PURITYFILTER_CPTMTS")
dir.create(outFolder, recursive=TRUE)

runFolder <- "."

#### retrieve the purity tagged tads
purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- paste0("Aran - ", pm)

corMet <- "pearson"
transfExpr <- "log10"
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)


purity_file <- file.path(runFolder,"ALLTADS_AND_PURITY_FINAL", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt$region_ID <- file.path(merge_dt$dataset, merge_dt$region)
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

purity_flagged_tads <- merge_dt$region_ID[merge_dt$purityCorr <= purityCorrThresh]

#### retrieve the purity tagged tads

inFolder <- "REVISION_CHANGES_CPTMTLABELS_ALLDS"
outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
tad2cptmt_dt <- get(load(outFile))

stopifnot(tad2cptmt_dt$tad_binaryCptmtLab == tad2cptmt_dt$start_binaryCptmtLab)
stopifnot(tad2cptmt_dt$tad_eightCptmtLab == tad2cptmt_dt$start_eightCptmtLab)
stopifnot(!duplicated(tad2cptmt_dt$region_ID))

stopifnot(merge_dt$region_ID %in% tad2cptmt_dt$region_ID)

tad_binaryCptmts <- setNames(tad2cptmt_dt$tad_binaryCptmtLab, tad2cptmt_dt$region_ID)
tad_eightCptmts <- setNames(tad2cptmt_dt$tad_eightCptmtLab, tad2cptmt_dt$region_ID)
tad_normRank <- setNames(tad2cptmt_dt$tadCptmtNormRank, tad2cptmt_dt$region_ID)

stopifnot(merge_dt$region_ID %in% names(tad_eightCptmts))
stopifnot(merge_dt$region_ID %in% names(tad_binaryCptmts))
stopifnot(merge_dt$region_ID %in% names(tad_normRank))

merge_dt$tad_eightCptmts <- tad_eightCptmts[paste0(merge_dt$region_ID)]
merge_dt$tad_binaryCptmts <- tad_binaryCptmts[paste0(merge_dt$region_ID)]
merge_dt$tad_normRank <- tad_normRank[paste0(merge_dt$region_ID)]


stopifnot(!is.na(merge_dt$purityCorr))

legTitle <- ""

plotTit <- ""
mysub <- ""


all_cptmt_vars <- c("tad_binaryCptmts","tad_eightCptmts")

all_plot_vars <- c("purityCorr")


for(cptmt_var in all_cptmt_vars){
  
  for(plot_var in all_plot_vars) {
    
    plotTit <- paste0(plot_var, " by ", cptmt_var)
    
    mySub <- paste0("# DS = ", length(unique(file.path(merge_dt$hicds,merge_dt$exprds))),
                    "; # TADs = ", nrow(merge_dt))
    
    pbox <- ggboxplot(merge_dt, 
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

all_cptmt_vars <- c("tad_normRank")
all_plot_vars <- c("purityCorr")

cptmt_var="tadCptmtNormRank"
plot_var="purityCorr"

for(cptmt_var in all_cptmt_vars){
  
  for(plot_var in all_plot_vars) {
    
    
    myx <- merge_dt[,paste0(cptmt_var)]
    myy <- merge_dt[,paste0(plot_var)]
    
    
    plotTit <- paste0( plot_var, " by ", cptmt_var)
    
    mySub <- paste0("# DS = ", length(unique(file.path(merge_dt$hicds,merge_dt$exprds))),
                    "; # TADs = ", nrow(merge_dt))
    
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

##### pie chart purity tagged by cptmts

stopifnot(purity_flagged_tads %in% merge_dt$region_ID)

pf_dt <- merge_dt[merge_dt$region_ID %in% purity_flagged_tads,]

all_cptmt_vars <- c("tad_binaryCptmts","tad_eightCptmts")

all_plot_vars <- c("purityCorr")

cptmt_var="tad_binaryCptmts"


for(cptmt_var in all_cptmt_vars){

  tmp_dt <- aggregate(as.formula(paste0("region_ID ~ ", cptmt_var)), data=merge_dt, FUN=length)
  colnames(tmp_dt)[colnames(tmp_dt) == "region_ID"] <- "nTADs"
  tmp_dt$ratioTADs <- tmp_dt$nTADs/nrow(merge_dt)
  tmp_dt$ratioTADs_lab <- paste0(round(tmp_dt$ratioTADs*100, 2), "%")
  
  myTit <- paste0("Dist. of all PF annotated TADs across ", cptmt_var)
  mysub <- paste0("# DS = ", length(unique(file.path(merge_dt$hicds, merge_dt$exprds))) , "; # TADs = ", nrow(merge_dt))
  
  p_dist <- ggplot(tmp_dt, aes_string(x="1", y="ratioTADs", fill=cptmt_var)) +
    geom_col() +
    geom_text_repel(aes(label = ratioTADs_lab), position = position_stack(vjust = 0.5))+
    coord_polar(theta = "y") + 
    ggtitle(myTit, subtitle=mysub) +
    theme_void() +    
    labs(fill="") +
    blank_theme
  
  outFile <- file.path(outFolder, paste0( "distAllTADs_by_", cptmt_var, "_pie.", plotType))
  ggsave(p_dist, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  agg_dt <- aggregate(as.formula(paste0("region_ID ~ ", cptmt_var)), data=pf_dt, FUN=length)
  colnames(agg_dt)[colnames(agg_dt) == "region_ID"] <- "nTADs"
  agg_dt$ratioTADs <- agg_dt$nTADs/nrow(pf_dt)
  
  agg_dt$ratioTADs_lab <- paste0(round(agg_dt$ratioTADs*100, 2), "%")
  
  myTit <- paste0("Dist. purity-tagged TADs across ", cptmt_var)
  mysub <- paste0("# DS = ", length(unique(file.path(pf_dt$hicds, pf_dt$exprds))) , "; # TADs = ", nrow(pf_dt))
  
  p_dist <- ggplot(agg_dt, aes_string(x="1", y="ratioTADs", fill=cptmt_var)) +
    geom_col() +
    geom_text_repel(aes(label = ratioTADs_lab), position = position_stack(vjust = 0.5))+
    coord_polar(theta = "y") + 
    ggtitle(myTit, subtitle=mysub) +
    theme_void() +    
    labs(fill="") +
    blank_theme
  
  outFile <- file.path(outFolder, paste0( "distPurityTaggedTADs_by_", cptmt_var, "_pie.", plotType))
  ggsave(p_dist, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
}



# 
# 
# 
# 
# 
# 
# 
# final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
# stopifnot(file.exists(final_table_file))
# final_table_DT <- get(load(final_table_file))
# final_table_DT$region_ID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
# stopifnot(!duplicated(final_table_DT$regionID))
# regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)
# signif_tads <- final_table_DT$regionID[final_table_DT$adjPvalComb <= tadSignifThresh]
# 
# logFC_values <- setNames(final_table_DT$meanLogFC, final_table_DT$region_ID)
# corr_values <- setNames(final_table_DT$meanCorr, final_table_DT$region_ID)
# 
# stopifnot(setequal(names(logFC_values), tad2cptmt_dt$region_ID))
# stopifnot(setequal(names(corr_values), tad2cptmt_dt$region_ID))
# 
# tad2cptmt_dt[,paste0("meanCorr")] <- corr_values[tad2cptmt_dt$region_ID]
# tad2cptmt_dt[,paste0("meanLogFC")] <- logFC_values[tad2cptmt_dt$region_ID]
# 
# 
# legTitle <- ""
# 
# tad2cptmt_dt$tad_binaryCptmtLab <- factor(tad2cptmt_dt$tad_binaryCptmtLab,
#                                           levels=as.character(sort(unique(as.character(
#                                             tad2cptmt_dt$tad_binaryCptmtLab)))))
# stopifnot(!is.na(tad2cptmt_dt$tad_binaryCptmtLab))
# 
# 
# tad2cptmt_dt$tad_eightCptmtLab <- factor(tad2cptmt_dt$tad_eightCptmtLab,
#                                          levels=as.character(sort(unique(as.character(
#                                            tad2cptmt_dt$tad_eightCptmtLab)))))
# stopifnot(!is.na(tad2cptmt_dt$tad_eightCptmtLab))
# 
# 
