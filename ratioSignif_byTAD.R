require(ggpubr)
require(foreach)
require(doMC)
registerDoMC(40)

# Rscript  ratioSignif_byTAD.R

outFolder <- "RATIOSIGNIF_BYTAD"
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

limma_signif <- 0.01
tad_signif <- 0.01

result_dt_1 <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
result_dt_2 <- get(load("CREATE_FINAL_TABLE_RANDOM//all_result_dt.Rdata"))
stopifnot(!result_dt_2$hicds %in% result_dt_1$hicds)
result_dt <- rbind(result_dt_1,result_dt_2)

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  
  tad_dt_file <- file.path(hicds,  "genes2tad", "all_assigned_regions.txt")
  tad_dt <- read.delim(tad_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo", "region", "region_start", "region_end"))
  
  g2t_dt_file <- file.path(hicds,  "genes2tad", "all_genes_positions.txt")
  g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "gene_start", "gene_end", "region"))
  g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
  g2t_tad_dt <- g2t_dt[grep("_TAD", g2t_dt$region),]
  
  geneList <- get(load(file.path("PIPELINE/OUTPUT_FOLDER",  hicds, exprds, "0_prepGeneData/pipeline_geneList.Rdata")))
  stopifnot( geneList %in% g2t_tad_dt$entrezID )
  g2t_tad_dt <- g2t_tad_dt[g2t_tad_dt$entrezID %in% geneList,]
  
  de_dt <- get(load(file.path("PIPELINE/OUTPUT_FOLDER",  hicds, exprds, "1_runGeneDE/DE_topTable.Rdata")))
  de_dt$entrezID <- geneList[as.character(de_dt$genes)]
  de_dt <- na.omit(de_dt)
  
  stopifnot(geneList %in% de_dt$entrezID)
  
  g2t_tad_de_dt <- merge(g2t_tad_dt, de_dt[,c("entrezID", "adj.P.Val")], all.x=TRUE, all.y=FALSE, by="entrezID")
  
  stopifnot(geneList %in% g2t_tad_de_dt$entrezID)
  
  
  tad_limmaSignif_dt_1 <- aggregate(adj.P.Val ~ region, data=g2t_tad_de_dt, FUN=function(x) mean(x < limma_signif) ) 
  colnames(tad_limmaSignif_dt_1) [colnames(tad_limmaSignif_dt_1) == "adj.P.Val"] <- "ratioLimmaSignif"
  
  tad_limmaSignif_dt_2 <- aggregate(adj.P.Val ~ region, data=g2t_tad_de_dt, FUN=function(x) mean(-log10(x)) ) 
  colnames(tad_limmaSignif_dt_2) [colnames(tad_limmaSignif_dt_2) == "adj.P.Val"] <- "meanLimmaPvalLog10"
  
  tad_limmaSignif_dt_12 <- merge(tad_limmaSignif_dt_1, tad_limmaSignif_dt_2, by="region")
  stopifnot(!is.na(tad_limmaSignif_dt_12))
  stopifnot(nrow(tad_limmaSignif_dt_1) == nrow(tad_limmaSignif_dt_2))
  
  
  tad_result_dt <- result_dt[result_dt$hicds == hicds & result_dt$exprds == exprds,]
  
  stopifnot(!duplicated(tad_result_dt$region))
  
  tad_result_dt$region_rank <- rank(tad_result_dt$adjPvalComb, ties="min")
  
  
  tad_limmaSignif_dt <- merge(tad_limmaSignif_dt_12, tad_result_dt[,c("region", "adjPvalComb", "region_rank")], all.x=TRUE, all.y=FALSE, by="region")
  stopifnot(!is.na(tad_limmaSignif_dt))
  stopifnot(nrow(tad_limmaSignif_dt) == nrow(tad_limmaSignif_dt_12))
  
  tad_limmaSignif_dt$hicds <- hicds
  tad_limmaSignif_dt$exprds <- exprds
  tad_limmaSignif_dt
  
}

all_dt$signif <- ifelse(all_dt$adjPvalComb <= tad_signif, "signif.", "not signif.")
all_dt$signif <- factor(all_dt$signif, levels = c("signif.", "not signif."))
stopifnot(!is.na(all_dt$signif))

all_dt$hicds_lab <- gsub(myHicds, "", all_dt$hicds)

p_box <- ggboxplot(data=all_dt, x="signif", y="ratioLimmaSignif", xlab="TAD signif.", ylab="ratio limma signif. genes by TAD")+
  ggtitle(paste0(myHicds, " - ", exprds), subtitle = paste0("gene p-val thresh <= ", limma_signif, "; TAD p-val thresh <= ", tad_signif))+
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


outFile <- file.path(outFolder, paste0(myHicds, "_obs_permut_ratioSignifByTAD.", plotType))
ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


p_box <- ggboxplot(data=all_dt, x="signif", y="meanLimmaPvalLog10", xlab="TAD signif.", ylab="mean limma p-val [-log10]")+
  ggtitle(paste0(myHicds, " - ", exprds), subtitle = paste0("TAD p-val thresh <= ", tad_signif))+
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


outFile <- file.path(outFolder, paste0(myHicds, "_obs_permut_meanLimmaPvalByTAD.", plotType))
ggsave(p_box, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



#############################################################################################################################

save(all_dt, file="all_dt.Rdata", version=2)

# nFAMs of topTADs

nTopRank <- 50

plot_dt <- all_dt
plot_dt <- plot_dt[plot_dt$region_rank <= nTopRank,]

sub <- paste0(nTopRank, " top-ranking TADs only")

outFile <- file.path(outFolder, paste0(myHicds, "_", exprds, "_meanLimmaPvalLog10_", nTopRank, "nTopTADs_boxplot.", plotType))

p_box <- ggboxplot(data=plot_dt, x="hicds_lab", y=paste0("meanLimmaPvalLog10")) +
  ggtitle(paste0(myHicds, " - ", exprds), subtitle = paste0(sub))+
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


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path(outFolder, paste0(myHicds, "_", exprds, "_meanLimmaPvalLog10_", nTopRank, "nTopTADs_density.", plotType))

do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
plot_multiDens(split(plot_dt$meanLimmaPvalLog10, plot_dt$hicds_lab),
               plotTit=paste0("TAD mean limma p-val [log10]"), legPos = "topright")
mtext(side=3, text = paste0(hicds, " - ", exprds, "; ", sub))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(paste0(txt))
cat(paste0("*** DONE: ", script_name, "\n"))











