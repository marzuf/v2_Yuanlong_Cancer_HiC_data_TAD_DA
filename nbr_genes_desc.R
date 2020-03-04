
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)



# Rscript nbr_genes_desc.R

plotType <- "png"
myHeight <- myWidth <- 400

plotTypeGG <- "svg"
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "NBR_GENES_DESC"
dir.create(outFolder, recursive = TRUE)
  
all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
hicds = all_hicds[1]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))


all_nbr_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    gl <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
    
    data.frame(
      hicds =hicds,
      exprds=exprds,
      nGenes = length(gl),
      stringsAsFactors = FALSE
    )
  }
  exprds_dt
}

all_nbr_dt$dataset <- file.path(all_nbr_dt$hicds, all_nbr_dt$exprds)
all_nbr_dt <- all_nbr_dt[order(all_nbr_dt$nGenes, decreasing = TRUE),]
all_nbr_dt$dataset <- factor(all_nbr_dt$dataset, levels = all_nbr_dt$dataset)

nDS <- length(unique(all_nbr_dt$dataset))

p_nbr <- ggbarplot(all_nbr_dt, 
                    x = "dataset", y = "nGenes", col = "darkblue", fill ="darkblue",
                    xlab = "Datasets ordered by decreasing # of genes")+
  labs(title = "# Genes", subtitle = paste0("all datasets - n = ", nDS))+
  scale_y_continuous(name=paste0("# genes used in pipeline"),
                     breaks = scales::pretty_breaks(n = 10))+
  theme( 
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    axis.line.x= element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent")
  )
outFile <- file.path(outFolder, paste0("allDS_nbrGenes_barplot.", plotTypeGG))
ggsave(p_nbr, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

