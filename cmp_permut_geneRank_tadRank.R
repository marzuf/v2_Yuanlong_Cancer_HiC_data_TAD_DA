options(scipen=100)

SSHFS=F

setDir <- "/media/electron"
setDir <- ""

# Rscript cmp_permut_geneRank_tadRank.R


script_name <- "cmp_permut_geneRank_tadRank.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(RColorBrewer)
require(reshape2)
require(ggplot2)

require(reshape2)

require(ggpubr)
require(ggsci)

# require(gplots)
registerDoMC(4)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "png"
myWidth <- myHeight <- 400
myWidthGG <- 9
myHeightGG <- 7
plotTypeGG <- "svg"

outFolder <- "CMP_PERMUT_GENERANK_TADRANK"
dir.create(outFolder, recursive = TRUE)

obs_rank_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
obs_rank_dt$dataset <- file.path(obs_rank_dt$hicds, obs_rank_dt$exprds)

permut_rank_dt <- get(load("GENE_RANK_TAD_RANK_PERMUT/all_gene_tad_signif_dt.Rdata"))
permut_rank_dt$dataset <- file.path(permut_rank_dt$hicds, permut_rank_dt$exprds)
permut_rank_dt$hicds_lab <- gsub("(.+)_.+_40kb", "\\1_40kb", permut_rank_dt$hicds)
permut_rank_dt$dataset_lab <- file.path(permut_rank_dt$hicds_lab, permut_rank_dt$exprds)


tad_signif_thresh <- 0.01
gene_signif_thresh <- 0.01

ds_withPerm <- intersect(permut_rank_dt$dataset_lab, obs_rank_dt$dataset)
plot_ds = ds_withPerm[1]

all_dt <- foreach(plot_ds = ds_withPerm, .combine='rbind') %dopar% {
  
  sub_obs_dt <- obs_rank_dt[obs_rank_dt$dataset == plot_ds,]
  stopifnot(nrow(sub_obs_dt)>0)
  
  sub_perm_dt <- permut_rank_dt[permut_rank_dt$dataset_lab == plot_ds,]
  all_perms <- unique(sub_perm_dt$dataset)
  rd=all_perms[1]
  
  permut_dt <- foreach(rd = all_perms, .combine='rbind') %do% {
    
    stopifnot(basename(plot_ds) == basename(rd))
    
    rd_type <- gsub(".+_(.+?_40kb)", "\\1", dirname(rd))
    
    rd_sub_perm_dt <- sub_perm_dt[sub_perm_dt$dataset == rd,]
    stopifnot(nrow(rd_sub_perm_dt) > 0)
    
    stopifnot(!duplicated(rd_sub_perm_dt$entrezID))
    stopifnot(!duplicated(sub_obs_dt$entrezID))
    
    merge_dt <- merge(sub_obs_dt[,c("entrezID", "tad_rank")],
                      rd_sub_perm_dt[,c("entrezID", "tad_rank")], 
                      by="entrezID",
                      all.x=FALSE, all.y=FALSE,
                      suffixes = c("_obs", "_perm")
                      )
    
    source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
    
    my_x <- merge_dt$tad_rank_obs
    my_y <- merge_dt$tad_rank_perm
    
    outFile <- file.path(outFolder, paste0(dirname(plot_ds), "_", basename(plot_ds), "_TADrank_", rd_type, "_vs_obs_densplot.", plotType))    
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    
    
    densplot(
      x=my_x,
      xlab = "TAD rank: obs.",
      ylab = paste0("TAD rank: ", rd_type),
      main = paste0(plot_ds),
      y=my_y,
      pch=16,
      cex=0.7
    )
    mtext(side=3, text=paste0("with perm. ", rd_type))
    addCorr(x=my_x, y=my_y, legPos = "topleft", bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))  
    
        
    
    merge_dt <- merge(sub_obs_dt[,c("entrezID", "tad_rank", "tad_adjCombPval", "adj.P.Val")],
                      rd_sub_perm_dt[,c("entrezID", "tad_rank", "tad_adjCombPval", "adj.P.Val")], 
                      by="entrezID",
                      all.x=FALSE, all.y=FALSE,
                      suffixes = c("_obs", "_perm")
    )
    merge_dt$tad_signif_perm <- merge_dt$tad_adjCombPval_perm <= tad_signif_thresh
    merge_dt$tad_signif_obs <- merge_dt$tad_adjCombPval_obs <= tad_signif_thresh
    
    signifPermOnly <- sum(merge_dt$tad_signif_perm & ! merge_dt$tad_signif_obs)
    signifObsOnly <- sum(! merge_dt$tad_signif_perm &  merge_dt$tad_signif_obs)
    signifPermObs <- sum(merge_dt$tad_signif_perm &  merge_dt$tad_signif_obs)
    
    stopifnot(sum(merge_dt$tad_signif_obs  | merge_dt$tad_signif_perm) == signifObsOnly+signifPermObs+signifPermOnly)
    
    data.frame(
      obs_ds = plot_ds,
      permut_dt = rd,
      signifPermOnly=signifPermOnly,
      signifObsOnly=signifObsOnly,
      signifPermObs=signifPermObs,
      stringsAsFactors = FALSE
    )
    
  }
  
  permut_dt
  
}

outFile <- file.path(outFolder, "all_dt.Rdata")
save(all_dt, file= outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
load(outFile)    


plot_dt <- melt(all_dt, id=c("obs_ds", "permut_dt"))

plot_dt$permut_type <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                         gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                              gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                   gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", dirname(plot_dt$permut_dt)))))

bar_plot <- ggbarplot(data = plot_dt, y = "value", x = "permut_type" , fill = "variable", 
          title ="# genes in signif. TADs obs/permut",
          position=position_dodge(),
          xlab = "", ylab="# genes in signif. TADs")+ 
  labs(fill ="")+
  scale_fill_nejm()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  facet_grid(~obs_ds, switch = "x") +
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"))

outFile <- file.path(outFolder, paste0("cmp_signif_genes_obs_permut_barplot.", plotTypeGG))
ggsave(bar_plot, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

