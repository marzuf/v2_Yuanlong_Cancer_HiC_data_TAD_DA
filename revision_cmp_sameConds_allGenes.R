# Rscript revision_cmp_sameConds_allGenes.R

require(ggpubr)
require(ggsci)
require(doMC)
require(foreach)
require(reshape2)

registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myWidth <- 500
myWidthGG <- 7
myHeightGG <- 5
myHeight <- 400

plotCex <- 1.2

buildTable <-TRUE

source("revision_settings.R")
all_pairs_dt2 <- all_pairs_dt[,c("hicds2", "hicds1", "exprds")]
colnames(all_pairs_dt2) <- c("hicds1", "hicds2", "exprds")
all_pairs_dt_all <- rbind(all_pairs_dt,all_pairs_dt2) 

cmp_dt <- get(load("CMP_SAME_CONDITIONS_V3/all_conds_dt.Rdata"))
cmp_dt$entrezID <- as.character(cmp_dt$entrezID)

outFolder <- file.path("REVISION_CMP_SAMECONDS_ALLGENES")
dir.create(outFolder, recursive=TRUE)

# signifThresh <- 0.01
signif_cmp_dt <- cmp_dt

# signif_cmp_dt <- cmp_dt[cmp_dt$tad_adjCombPval_hicds1 <= signifThresh &
#                           cmp_dt$tad_adjCombPval_hicds2 <= signifThresh ,
#                           ]
signif_cmp_dt$ds_pair <- file.path(signif_cmp_dt$hicds_hicds1, signif_cmp_dt$hicds_hicds2, signif_cmp_dt$exprds)

all_dspairs <- unique(signif_cmp_dt$ds_pair)

cmp  = all_dspairs[1]

if(buildTable){
  
  tad_genepairs_dt <- foreach(cmp  = all_dspairs, .combine='rbind') %dopar% {
    
    ds_signif_cmp_dt <- signif_cmp_dt[signif_cmp_dt$ds_pair == cmp,]
    
    intersect_signif_genes <- ds_signif_cmp_dt$entrezID
    
    stopifnot(length(intersect_signif_genes) > 0)
    
    tad=ds_signif_cmp_dt$region_hicds1[1]
    
    ds1_pairs_dt <- foreach(tad = unique(ds_signif_cmp_dt$region_hicds1), .combine='rbind') %dopar% {
      tad_genes <- sort(ds_signif_cmp_dt$entrezID[ds_signif_cmp_dt$region_hicds1 == tad])
      stopifnot(length(tad_genes) > 0)
      if(length(tad_genes) < 2) return(NULL)
      gene_combs <- combn(tad_genes, m=2)
      stopifnot(nrow(gene_combs) == 2)
      data.frame(
        hicds1 = unique(ds_signif_cmp_dt$hicds_hicds1[ds_signif_cmp_dt$region_hicds1 == tad]),
        tad1 = tad,
        gene1 = gene_combs[1,],
        gene2 = gene_combs[2,],
        stringsAsFactors = FALSE
      )
    }
    stopifnot(ds1_pairs_dt$gene1 < ds1_pairs_dt$gene2)
    stopifnot(!duplicated(file.path(ds1_pairs_dt$gene1, ds1_pairs_dt$gene2)))
    ds1_pairs_dt$gene_pair <- file.path(ds1_pairs_dt$gene1, ds1_pairs_dt$gene2)
    
    ds2_pairs_dt <- foreach(tad = unique(ds_signif_cmp_dt$region_hicds2), .combine='rbind') %dopar% {
      tad_genes <- sort(ds_signif_cmp_dt$entrezID[ds_signif_cmp_dt$region_hicds2 == tad])
      stopifnot(length(tad_genes) > 0)
      if(length(tad_genes) < 2) return(NULL)
      gene_combs <- combn(tad_genes, m=2)
      stopifnot(nrow(gene_combs) == 2)
      data.frame(
        hicds2 = unique(ds_signif_cmp_dt$hicds_hicds2[ds_signif_cmp_dt$region_hicds2 == tad]),
        tad2 = tad,
        gene1 = gene_combs[1,],
        gene2 = gene_combs[2,],
        stringsAsFactors = FALSE
      )
    }
    stopifnot(ds2_pairs_dt$gene1 < ds2_pairs_dt$gene2)
    stopifnot(!duplicated(file.path(ds2_pairs_dt$gene1, ds2_pairs_dt$gene2)))
    ds2_pairs_dt$gene_pair <- file.path(ds2_pairs_dt$gene1, ds2_pairs_dt$gene2)
    
    
    tot_nbrPairs <- length(unique(c(ds1_pairs_dt$gene_pair, ds2_pairs_dt$gene_pair)))
    onlyHicds1_nbrPairs <- sum(! ds1_pairs_dt$gene_pair %in% ds2_pairs_dt$gene_pair)
    onlyHicds2_nbrPairs <- sum(! ds2_pairs_dt$gene_pair %in% ds1_pairs_dt$gene_pair)
    intersect_nbrPairs <- length(intersect(ds2_pairs_dt$gene_pair,ds1_pairs_dt$gene_pair))
    stopifnot(intersect_nbrPairs+onlyHicds1_nbrPairs+onlyHicds2_nbrPairs == tot_nbrPairs)
    
    
    data.frame(
      hicds1 = dirname(dirname(cmp)),
      hicds2 = basename(dirname(cmp)),
      exprds = basename(cmp),
      tot_nbrPairs=tot_nbrPairs,
      onlyHicds1_nbrPairs=onlyHicds1_nbrPairs,
      onlyHicds2_nbrPairs=onlyHicds2_nbrPairs,
      intersect_nbrPairs=intersect_nbrPairs,
      stringsAsFactors = FALSE
    )
  }
  outFile <- file.path(outFolder, "tad_genepairs_dt.Rdata")
  save(tad_genepairs_dt, file =outFile, version=2)
  cat("... written: ", outFile, "\n")
  
} else { # if/else
  outFile <- file.path(outFolder, "tad_genepairs_dt.Rdata")
  tad_genepairs_dt <- get(load(outFile))
  # load("REVISION_CMP_SAMECONDS/tad_genepairs_dt.Rdata")
}

tad_genepairs_dt$ratioShared <- tad_genepairs_dt$intersect_nbrPairs/tad_genepairs_dt$tot_nbrPairs
tad_genepairs_dt$ratioOnlyHicds1 <- tad_genepairs_dt$onlyHicds1_nbrPairs/tad_genepairs_dt$tot_nbrPairs
tad_genepairs_dt$ratioOnlyHicds2 <- tad_genepairs_dt$onlyHicds2_nbrPairs/tad_genepairs_dt$tot_nbrPairs

## sub is to take only the norm vs tumor DS

for(plot_type in c("all", "sub")) {
  
  if(plot_type == "all") {
    sub_tad_genepairs_dt <- tad_genepairs_dt
  }else if(plot_type == "sub") {
    sub_tad_genepairs_dt <- tad_genepairs_dt[tad_genepairs_dt$hicds1 %in% all_pairs_dt_all$hicds1 &
                                               tad_genepairs_dt$hicds2 %in% all_pairs_dt_all$hicds2 &
                                               tad_genepairs_dt$exprds %in% all_pairs_dt_all$exprds,
                                               ]
  }
  
  plot_dt <- melt(sub_tad_genepairs_dt[,c("ratioShared", "ratioOnlyHicds1", "ratioOnlyHicds2")])
  
  plotTit <- paste0("Gene pairs from signif. TADs (", plot_type, " DS)")
  
  all_cmps <- unique(file.path(sub_tad_genepairs_dt$hicds1, sub_tad_genepairs_dt$hicds2, sub_tad_genepairs_dt$exprds))
  
  # mySub <- paste0("# DS comparisons = ", length(all_cmps), "; signif: p-val <= ", signifThresh)
  mySub <- paste0("# DS comparisons = ", length(all_cmps), "; all genes")
  
  legTitle <- ""
  
  p3 <- ggdensity(plot_dt,
                  x = paste0("value"),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0("Ratio of gene pairs"),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "variable",
                  fill = "variable",
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = mySub)+
    # scale_color_manual(values=my_cols)+
    # scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    mytheme
  
  outFile <- file.path(outFolder, paste0("ratio_shared_gene_pairs_fromSignifTADs_", plot_type, "ds_density.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  p3 <- ggboxplot(plot_dt,
                  x = paste0("variable"),
                  y = "value",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0("Ratio of gene pairs"),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "variable",
                  #fill = "variable",
                  add="jitter",
                  outlier.shape=NA,
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = mySub)+
    # scale_color_manual(values=my_cols)+
    # scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
   # guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    mytheme
  
  outFile <- file.path(outFolder, paste0("ratio_shared_gene_pairs_fromSignifTADs_", plot_type, "ds_boxplot.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
}

