require(ggpubr)
require(foreach)
require(doMC)
registerDoMC(40)

# Rscript  dist_signif_to_borders.R

outFolder <- "DIST_SIGNIF_TO_BORDERS"
dir.create(outFolder, recursive = TRUE)

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 9

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]

myHicds <- "ENCSR489OCU_NCI-H460"

all_hicds <- all_hicds[grepl(myHicds, all_hicds)]


hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"

signifThresh <- 0.01

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  
  tad_dt_file <- file.path(hicds,  "genes2tad", "all_assigned_regions.txt")
  tad_dt <- read.delim(tad_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo", "region", "region_start", "region_end"))
  
  g2t_dt_file <- file.path(hicds,  "genes2tad", "all_genes_positions.txt")
  g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "gene_start", "gene_end", "region"))
  g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
  g2t_tad_dt <- g2t_dt[grep("_TAD", g2t_dt$region),]
  
  g2t_tad_pos_dt <- merge(g2t_tad_dt, tad_dt, by=c("chromo", "region"), all.x=T, all.y=F)
  stopifnot(!is.na(g2t_tad_pos_dt))
  
  g2t_tad_pos_dt$minDist_geneStart_border <- pmin(abs(g2t_tad_pos_dt$region_start - g2t_tad_pos_dt$gene_start),
                                                  abs(g2t_tad_pos_dt$region_end - g2t_tad_pos_dt$gene_start))
  g2t_tad_pos_dt$minDist_geneStart_border_log10 <- log10(g2t_tad_pos_dt$minDist_geneStart_border)
  
  de_dt <- get(load(file.path("PIPELINE/OUTPUT_FOLDER",  hicds, exprds, "1_runGeneDE/DE_topTable.Rdata")))
  de_dt$entrezID <- as.character(de_dt$genes)
  
  g2t_signif_dt <- merge(g2t_tad_pos_dt, de_dt[,c("entrezID", "adj.P.Val")], all=FALSE, by=c("entrezID"))
  
  # plot(log10(g2t_signif_dt$minDist_geneStart_border), -log10(g2t_signif_dt$adj.P.Val))
  
  g2t_signif_dt$signif <- ifelse(g2t_signif_dt$adj.P.Val <= signifThresh, "signif.", "not signif.")
  
  # boxplot(
  #   minDist_geneStart_border_log10~signif, 
  #   ylab = "min dist. gene start <-> TAD border [bp - log10]",
  #   xlab = "",
  #   data=g2t_signif_dt)
  
  g2t_signif_dt$hicds <- hicds
  g2t_signif_dt$exprds <- exprds
  g2t_signif_dt
}

all_dt$hicds_lab <- gsub(myHicds, "", all_dt$hicds)

p_box <- ggboxplot(data=all_dt, x="signif", y="minDist_geneStart_border_log10", xlab="", ylab="min dist. gene start <-> TAD border [bp log10]")+
  ggtitle(myHicds, subtitle = paste0("gene adj. p-val <= ", signifThresh))+
  facet_grid(~hicds_lab, switch="x") +
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


outFile <- file.path(outFolder, paste0(myHicds, "_obs_permut_minDist_geneStart_border.", plotType))
ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path(outFolder, paste0(myHicds, "_obs_permut_minDist_geneStart_border_density.", plotType))
do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG*1.2))
par(bty="L")
plot_multiDens(

  split(all_dt$minDist_geneStart_border_log10, all_dt$hicds_lab),
  legPos = "topleft",
  plotTit = paste0(myHicds))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



