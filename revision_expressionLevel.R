require(foreach)
require(doMC)
registerDoMC(40)
require(ggsci)
require(ggpubr)
require(ggrepel)
require(ggplot2)
require(patchwork)
require(stringr)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


# Rscript revision_expressionLevel.R

setDir <- "/media/electron"
setDir <- ""

plotType <- "png"
myHeightGG <- 5
myWidthGG <- 7
myHeight <- 400
myWidth <- 500

plotCex <- 1.2

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- "REVISION_EXPRESSION_LEVEL"
dir.create(outFolder)

tadSignifThresh <- 0.01

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_dt <- get(load(final_table_file))
final_table_DT <- final_dt
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)
final_table_DT$signif_lab <- ifelse(final_table_DT$adjPvalComb <= tadSignifThresh, "signif.", "not signif.")
regionID_signif <- setNames(final_table_DT$signif_lab,final_table_DT$regionID )

ngenes1 <- str_count(final_table_DT$region_genes, pattern=",")+1
names(ngenes1) <- final_table_DT$regionID

expr_level_dt <- get(load("REVISION_EXPRESSION_LEVELDATA/all_exprLevel_dt.Rdata"))
expr_level_dt$regionID <- file.path(expr_level_dt$hicds, expr_level_dt$exprds, expr_level_dt$region)

ngenes2 <- setNames(as.numeric(table(expr_level_dt$regionID)), names(table(expr_level_dt$regionID)))

stopifnot(length(setdiff(names(ngenes1), names(ngenes2))) == 0)

stopifnot(ngenes1 == ngenes2[names(ngenes1)])

# > okay, now it is checked, in expr_level_dt i have only the genes used in the pipeline, i can aggregate directly this table

legTitle <- ""

mytheme <-     theme(
  # text = element_text(family=fontFamily),
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

all_vars <- c("aggLog10Expr", "zscoreAggExpr", "qqnormAggExpr", "histqt1AggExpr" ,"histqt2AggExpr")

expr_var <- "aggLog10Expr"
aggFun <- "mean"

for(expr_var in all_vars) {
  
  agg_dt <- aggregate(as.formula(paste0(expr_var, "~regionID")), data = expr_level_dt, FUN=aggFun)
  out_dt <- agg_dt
  stopifnot(agg_dt$regionID %in% names(regionID_signif))
  stopifnot(agg_dt$regionID %in% names(regionID_pvals))
  agg_dt$signif_lab <- regionID_signif[paste0(agg_dt$regionID)]
  agg_dt$pval <- regionID_pvals[paste0(agg_dt$regionID)]
  out_dt$adjPval <- regionID_pvals[paste0(out_dt$regionID)]
  
  outFile <- file.path(outFolder, paste0(expr_var, "_aggByTAD_", aggFun, ".Rdata"))
  save(out_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  stopifnot(!is.na(agg_dt$signif_lab))
  
  plotTit <- paste0("Expression level and TAD signif.: ", expr_var)
  
  mySub <- paste0("# DS = ", length(unique(dirname(agg_dt$regionID))), "; # TADs = ", length(unique(agg_dt$regionID)),
                  " (# signif. = ", sum(agg_dt$signif_lab == "signif."), ")")
  
  my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(agg_dt$signif_lab))
  
  
  p3 <- ggdensity(agg_dt,
                  x = paste0(expr_var),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0(expr_var, " (TAD ", aggFun, ")"),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "signif_lab",
                  fill = "signif_lab",
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = mySub)+
    scale_color_manual(values=my_cols)+
    scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    mytheme
  
  outFile <- file.path(outFolder, paste0("tad_", aggFun, "_", expr_var, "_signif_notsignif_density.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  plotTit <- paste0("Expression level and TAD adjPval: ", expr_var)
  
  mySub <- paste0("# DS = ", length(unique(dirname(agg_dt$regionID))), "; # TADs = ", length(unique(agg_dt$regionID)),
                  " (# signif. = ", sum(agg_dt$pval <= tadSignifThresh), ")")
  
  myx <- agg_dt[,paste0(expr_var)]
  myy <- -log10(agg_dt$pval)
  outFile <- file.path(outFolder, paste0("tad_", aggFun, "_", expr_var, "_vs_adjPval.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myHeight))
  densplot(
    x= myx,
    y=myy,
    xlab =paste0(expr_var, " (TAD ", aggFun, ")"),
    ylab = "adj. pval [-log10]",
    cex.main=plotCex,
    cex.main=plotCex,
    cex.main=plotCex,
    main=plotTit
  )
  mtext(side=3, text=mySub)
  addCorr(x=myx,y=myy, legPos="topright", bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

all_vars <- c("histqt1AggExpr" ,"histqt2AggExpr")
all_vars <- c("histqt1AggExpr")
expr_var <- "histqt1AggExpr"


for(expr_var in all_vars) {
  
  # for each TAD, count how many genes in the quantiles (take the ratio)
  agg_tad_dt <- aggregate(as.formula(paste0("symbol~regionID + ", expr_var)), data=expr_level_dt, FUN=length)
  colnames(agg_tad_dt)[colnames(agg_tad_dt) == "symbol"] <- "nGenes"
  stopifnot(agg_tad_dt$regionID %in% names(regionID_signif))
  agg_tad_dt$signif_lab <- regionID_signif[paste0(agg_tad_dt$regionID)]
  
  stopifnot(agg_tad_dt$regionID %in% names(ngenes2))
  agg_tad_dt$nTotGenes <- ngenes2[paste0(agg_tad_dt$regionID)]
  agg_tad_dt$ratioGenes <- agg_tad_dt$nGenes/agg_tad_dt$nTotGenes
  tmpfoo <- aggregate(ratioGenes~regionID, FUN=sum, data=agg_tad_dt)
  stopifnot(round(tmpfoo$ratioGenes,10) == 1)
  
  stopifnot(range(agg_tad_dt$ratioGenes) >= 0)
  stopifnot(range(agg_tad_dt$ratioGenes) <= 1)
        
  # on average, how many genes by quantile ?
  ## !!! cannot take the mean in aggregate, because I don't have the rows when 0 gene in a given hist break
  # agg_dt <- aggregate(as.formula(paste0("ratioGenes~signif_lab + ", expr_var)), data=agg_tad_dt, FUN=mean)
  # colnames(agg_dt)[colnames(agg_dt) == "ratioGenes"] <- "mean_ratioGenes"
  tmp_agg_dt <- aggregate(as.formula(paste0("ratioGenes~signif_lab + ", expr_var)), data=agg_tad_dt, FUN=sum)
  colnames(tmp_agg_dt)[colnames(tmp_agg_dt) == "ratioGenes"] <- "sum_ratioGenes"
  
  aggfoo <- aggregate(regionID~signif_lab, FUN=function(x)length(unique(x)), data=agg_tad_dt)
  nTADs <- setNames(aggfoo$regionID, aggfoo$signif_lab)
  tmp_agg_dt$nTot <- nTADs[paste0(tmp_agg_dt$signif_lab)]
  
  tmp_agg_dt$mean_ratioGenes <- tmp_agg_dt$sum_ratioGenes/tmp_agg_dt$nTot
  
  stopifnot(round(sum(tmp_agg_dt$mean_ratioGenes[tmp_agg_dt$signif_lab == "not signif."]), 6) == 1)
  stopifnot(round(sum(tmp_agg_dt$mean_ratioGenes[tmp_agg_dt$signif_lab == "signif."]), 6) == 1)
  
  agg_dt <- tmp_agg_dt
  
  agg_dt[, paste0(expr_var)] <- factor(as.character(agg_dt[, paste0(expr_var)]), 
                                       levels = as.character(min(agg_dt[, paste0(expr_var)]): max(agg_dt[, paste0(expr_var)])))
  stopifnot(!is.na(  agg_dt[, paste0(expr_var)] ))
  
  
  ############################################################## DRAW DOUBLE PIEPLOT
  
  plot_dt <- agg_dt[agg_dt$signif_lab == "signif.",]
  # plot_dt <- plot_dt[order(as.numeric(plot_dt$histqt1AggExpr)),]
  plot_dt$freq <- 100*plot_dt$mean_ratioGenes
  plot_dt$freq_rd <- paste0(round(100*plot_dt$mean_ratioGenes,2), "%")
  
  p_signif <- ggplot(plot_dt, aes_string(x="1", y="freq", fill=expr_var)) +
    ggtitle( "signif. TADs")+
    geom_col() +
    # geom_text(aes(label = freq_rd), position = position_stack(vjust = 0.5))+
    geom_text_repel(aes(label = freq_rd), position = position_stack(vjust = 0.5))+
    coord_polar(theta = "y") + 
    theme_void() +    
    labs(fill="hist. break") +
    blank_theme

  plot_dt <- agg_dt[agg_dt$signif_lab == "not signif.",]
  # plot_dt <- plot_dt[order(as.numeric(plot_dt$histqt1AggExpr)),]
  plot_dt$freq <- 100*plot_dt$mean_ratioGenes
  plot_dt$freq_rd <- paste0(round(100*plot_dt$mean_ratioGenes,2), "%")
  
  p_notsignif <- ggplot(plot_dt, aes_string(x="1", y="freq", fill=expr_var)) +
    ggtitle( "not signif. TADs")+
    geom_col() +
    geom_text_repel(aes(label = freq_rd), position = position_stack(vjust = 0.5))+
    coord_polar(theta = "y") + 
    theme_void() +    
    labs(fill="hist. break") +
    blank_theme
  
  
  p_tmp <- (p_signif + p_notsignif)
  
  stopifnot(nTADs["signif."] == sum(final_table_DT$adjPvalComb<=tadSignifThresh))
           
  
  plot_tit <- "mean ratios"
  subtit <- paste0("# TADs: ", paste0(names(nTADs), "=", nTADs, collapse="; "), " adjPval<=", tadSignifThresh)

    
  p_out <- p_tmp + 
    plot_layout(guides = 'collect')+
    plot_annotation(
      title = plot_tit,
      subtitle = paste0(subtit) 
    ) & theme(plot.title=element_text(hjust=0.5, size=14, face = "bold"), 
              plot.subtitle=element_text(hjust=0.5, size=12, face = "italic"),
              legend.position = 'top')
  
  outFile <- file.path(outFolder,paste0("mean_ratioGenes_byHistBreak_", expr_var, "_signif_notsignif_pieplot.", plotType))
  ggsave(p_out, filename = outFile, height=myHeightGG, width=myWidthGG*1.5)
  cat(paste0("... written: ", outFile, "\n"))
  
  ############################################################## DRAW BARPLOT
  legTitle <- paste0("hist. break")
    
  ggbar_p <-  ggbarplot(agg_dt, 
                        y="mean_ratioGenes",
                        x="signif_lab", 
                        fill=paste0(expr_var)) + 
    ggtitle(plotTit, subtitle=mySub)+
    mytheme +
    labs(x="" , y ="ratio of TADs", color=paste0(legTitle),fill=paste0(legTitle)) + 
    theme(
      axis.text.x = element_text(hjust=1, vjust=0.5,size=10,angle=90)
    )
  outFile <- file.path(outFolder,paste0("mean_ratioGenes_byHistBreak_", expr_var, "_barplot.", plotType))
  ggsave(ggbar_p, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
}





