
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)



# Rscript tad_size_desc.R

plotType <- "png"
myHeight <- myWidth <- 400

plotTypeGG <- "svg"
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "TAD_SIZE_DESC"
dir.create(outFolder, recursive = TRUE)
  
all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
hicds = all_hicds[1]

all_sizes_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  tad_dt <- read.delim(file.path(hicds, "genes2tad", "all_assigned_regions.txt"), header=FALSE, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
  stopifnot(!duplicated(tad_dt$region))
  tad_dt$hicds <- hicds
  tad_dt
}


all_sizes_dt$size <- all_sizes_dt$end - all_sizes_dt$start + 1
all_sizes_dt$size_kb <- all_sizes_dt$size/1000

nTADs_dt <- aggregate(region~hicds, data=all_sizes_dt, FUN=length)
nTADs_dt <- nTADs_dt[order(nTADs_dt$region, decreasing = TRUE),]
ds_ord <- nTADs_dt$hicds
nTADs_dt$hicds <- factor(nTADs_dt$hicds, levels = ds_ord)
stopifnot(!is.na(nTADs_dt$hicds))

nDS <- length(unique(all_sizes_dt$hicds))

p_nbr <- ggbarplot(nTADs_dt, 
                    x = "hicds", y = "region", col = "darkblue", fill ="darkblue",
                    xlab = "Datasets ordered by decreasing # of TADs")+
  labs(title = "# TADs", subtitle = paste0("all datasets - n = ", nDS))+
  scale_y_continuous(name=paste0("# TADs"),
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
outFile <- file.path(outFolder, paste0("allDS_tad_nbr_barplot.", plotTypeGG))
ggsave(p_nbr, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))




all_sizes_dt$hicds <- factor(all_sizes_dt$hicds, levels = ds_ord)
stopifnot(!is.na(all_sizes_dt$hicds))

nDS <- length(unique(all_sizes_dt$hicds))

p_size <- ggboxplot(all_sizes_dt, 
       x = "hicds", y = "size_kb",
       xlab = "Datasets ordered by decreasing # of TADs")+
  labs(title = "TAD size", subtitle = paste0("all datasets - n = ", nDS))+
  scale_y_continuous(name=paste0("TAD size [kb]"),
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
outFile <- file.path(outFolder, paste0("allDS_tad_size_boxplot.", plotTypeGG))
ggsave(p_size, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("allDS_tad_size_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot(density(
  all_sizes_dt$size_kb
), main = paste0("all hicds - n=", length(all_hicds), " - TAD size [kb]"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("allDS_tad_nbr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot(density(
  nTADs_dt$region
), main = paste0("all hicds - n=", nrow(nTADs_dt), " - # TADs"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
