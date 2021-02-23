
# Rscript revision_annot_hks_tfs.R

require(foreach)
require(doMC)
registerDoMC(40)

require(ggplot2)
require(ggpubr)
require(ggsci)

minGenes <- 3

buildTable <- TRUE


plotType <- "png"
myHeightGG <- 6
myWidthGG <- 7

outFolder <- file.path("REVISION_ANNOT_HKS_TFS")
dir.create(outFolder, recursive=TRUE)

# data from https://www.tau.ac.il/~elieis/HKG/HK_genes.txt (Eisenberg and Levanon 2013, updated)
hk_dt <- read.delim("HK_genes.txt", header=FALSE, col.names=c("symbol", "id"), stringsAsFactors = FALSE, sep=" ")
# a bit weird sep: " \t"... will put the \t in the id but ok not used
nrow(hk_dt)
# 3804
all_hks <- as.character(hk_dt$symbol)

# list of transcription factors:  http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt, 23.02.21
tf_dt <- read.delim("TF_names_v_1.01.txt", header=FALSE, col.names="symbol")
nrow(tf_dt)
# 1639
all_tfs <- as.character(tf_dt$symbol)

purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- paste0("Aran - ", pm)

corMet <- "pearson"
transfExpr <- "log10"
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)

runFolder <- "."

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

stopifnot(any(all_hks %in% gff_dt$symbol))
mean(all_hks %in% gff_dt$symbol)
# [1] 0.995531
stopifnot(any(all_tfs %in% gff_dt$symbol))
mean(all_tfs %in% gff_dt$symbol)
# [1] 0.9853569

purity_file <- file.path(runFolder,"ALLTADS_AND_PURITY_FINAL", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt$region_id <- file.path(merge_dt$dataset, merge_dt$region)
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))
purity_flagged_tads <- merge_dt$region_id[merge_dt$purityCorr <= purityCorrThresh]

# discard the flagged -> means that I keep those for which no puritycorr annot
merge_dt_all <- merge(agg_purity, resultData, by=c("dataset", "region"), all=TRUE)  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt_all$signif <- merge_dt_all$adjPvalComb <= signifThresh
merge_dt_all$region_id <- file.path(merge_dt_all$dataset, merge_dt_all$region)
stopifnot(purity_flagged_tads %in%  merge_dt_all$region_id)

plot_dt <- merge_dt_all
plot_dt$tad_purity_label <- ifelse(! plot_dt$region_id %in% merge_dt$region_id, "not known",
                                   ifelse(plot_dt$region_id %in% purity_flagged_tads, "purity-flagged", "not PF"))
stopifnot(sum(plot_dt$tad_purity_label == "not known") == nrow(merge_dt_all) - nrow(merge_dt))

stopifnot(nrow(plot_dt) > 0)

my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(plot_dt$tad_purity_label))


i=1
if(buildTable) {
  gene_annot <- foreach(i = 1:nrow(plot_dt)) %dopar% {
    region_symbols <- as.character(unlist(strsplit(plot_dt$region_genes[i], split=",")))
    stopifnot(region_symbols %in% gff_dt$symbol)
    nTFs <- sum(region_symbols %in% all_tfs)
    nHKs <- sum(region_symbols %in% all_hks)
    nGenes <- length(region_symbols)
    stopifnot(nGenes >= minGenes)
    list(
      nGenes=nGenes,
      nTFs=nTFs,
      nHKs=nHKs
    )
  }
  stopifnot(length(gene_annot) == nrow(plot_dt))
  
  plot_dt$nGenes <- unlist(lapply(gene_annot, function(x) x[["nGenes"]]))
  plot_dt$nTFs <- unlist(lapply(gene_annot, function(x) x[["nTFs"]]))
  plot_dt$nHKs <- unlist(lapply(gene_annot, function(x) x[["nHKs"]]))
  
  stopifnot(plot_dt$nTFs <= plot_dt$nGenes)
  stopifnot(plot_dt$nHKs <= plot_dt$nGenes)
  
  plot_dt$ratioTFs <- plot_dt$nTFs/plot_dt$nGenes
  plot_dt$ratioHKs <- plot_dt$nHKs/plot_dt$nGenes
  
  # stopifnot(!is.na(plot_dt)) # not true because NA for purity corr missing
  outFile <-  file.path(outFolder, "plot_dt.Rdata")
  save(plot_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <-  file.path(outFolder, "plot_dt.Rdata")
  plot_dt <- get(load(outFile))
  
}
all_vars <- c("TFs", "HKs")
plot_var = "TFs"


plot_dt$signif_lab <- ifelse(plot_dt$signif, "signif.", "not signif.")

signif_vars <- unique(plot_dt$signif_lab)
signif_var=signif_vars[1]

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

for(plot_var in all_vars) {
  
  plotTit <- paste0("TAD genes annotation: ", plot_var)
  
  mySub <- paste0("# DS = ", length(unique(plot_dt$dataset)), "; # TADs = ", length(unique(plot_dt$region_id)) )
  
  p3 <- ggdensity(plot_dt,
                  x = paste0("ratio", plot_var),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0("Ratio of ", plot_var, " in  TAD genes"),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "tad_purity_label",
                  fill = "tad_purity_label",
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
  
  outFile <- file.path(outFolder, paste0("tad_genes_annot_", plot_var, "_density.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  for(signif_var in signif_vars) {
    
    sub_dt <- plot_dt[plot_dt$signif_lab == signif_var,]
    stopifnot(nrow(sub_dt) > 0)
    plotTit <- paste0(signif_var, " TAD genes annotation: ", plot_var)
    
    mySub <- paste0("# DS = ", length(unique(sub_dt$dataset)), "; # TADs = ", length(unique(sub_dt$region_id)) )
    
    p3 <- ggdensity(sub_dt,
                    x = paste0("ratio", plot_var),
                    y = "..density..",
                    # combine = TRUE,                  # Combine the 3 plots
                    xlab = paste0("Ratio of ", plot_var, " in  TAD genes"),
                    # add = "median",                  # Add median line.
                    rug = FALSE,                      # Add marginal rug
                    color = "tad_purity_label",
                    fill = "tad_purity_label",
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
    
    outFile <- file.path(outFolder, paste0(gsub(" ", "", signif_var), "TADs_tad_genes_annot_", plot_var, "_density.", plotType))
    ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
  
  
}
