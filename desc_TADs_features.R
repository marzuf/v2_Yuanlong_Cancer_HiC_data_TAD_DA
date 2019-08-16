options(scipen=100)

# Rscript desc_TADs_features.R

script_name <- "desc_TADs_features.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(lattice)
require(ggplot2)

registerDoMC(1)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 8

script0_name <- "0_prepGeneData"


mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

hicds = all_hicds[1]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
exprds = all_exprds[1]
outFolder <- "DESC_TADS_FEATURES"
dir.create(outFolder, recursive = TRUE)

hicds = all_hicds[1]

####################################################################################################################################### >>> collect # of TADs and TAD size
all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(hicds_file))
  
  tadpos_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
  tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]
  stopifnot(nrow(tadpos_dt) > 0 )
  
  all_cols <- colnames(tadpos_dt)
  tadpos_dt$hicds <- hicds
  tadpos_dt <- tadpos_dt[, c("hicds", all_cols)]
  
  tadpos_dt
  
}
outFile <- file.path(outFolder, "all_dt.Rdata")
save(all_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
######################################################### >>> plot TAD size

all_dt$size <- all_dt$end - all_dt$start + 1

all_dt$size_log10 <- log10(all_dt$size)

#### PLOT THE SIGNIF PVAL FEATURES => compare dist Pval signif. vs. not signif.

# length(unique(all_dt$hicds)) # 18

nRows <- 5
nCols <- 6

all_vars <- c("size", "size_log10")
plot_var=all_vars[1]
for(plot_var in all_vars) {
  outFile <- file.path(outFolder, paste0("all_hicds_dist_TAD_", plot_var, "_density_lattice.", plotType))
  do.call(plotType, list(outFile, height=myHeight*3, width=myWidth*3))
  # myplot <- densityplot( formula(paste0("~",plot_var)), groups = hicds, data = all_dt, #auto.key = TRUE,
  myplot <- densityplot( formula(paste0("~",plot_var, "| hicds")), data = all_dt, #auto.key = TRUE,
                         # par.strip.text=list(cex=1), # width of the strip bar
                         par.strip.text = list(cex = 0.8, font = 4, col = "brown"),
                         layout = c(nCols, nRows), # column,row
                         scales=list(y=list(relation="free"),
                                     x=list(relation="free")
                         ),
                         main = paste0(plot_var))
  print(myplot)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}



######################################################### >>> boxplot size

all_dt$datasetName <- gsub("_40kb", "", all_dt$hicds)
stopifnot(is.numeric(all_dt$size))
stopifnot(is.numeric(all_dt$size_log10))

mean_dt <- aggregate(size ~ datasetName, FUN=mean, data=all_dt)
mean_dt <- mean_dt[order(mean_dt$size, decreasing = TRUE),]

all_dt$datasetName <- factor(as.character(all_dt$datasetName), levels=as.character(mean_dt$datasetName))

plot_vars <- c("size", "size_log10")

plot_var=plot_vars[1]

for(plot_var in plot_vars) {
  
  p_var <-  ggplot(all_dt, aes_string(x = "datasetName", y = plot_var)) + 
    geom_boxplot()+
    # coord_cartesian(expand = FALSE) +
    ggtitle(paste0("TAD ", plot_var))+
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(plot_var),
                       breaks = scales::pretty_breaks(n = 20))+
    theme( # Increase size of axis lines
      strip.text = element_text(size = 12),
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"),
      strip.text.x = element_text(size = 10),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.text.x =  element_text(color="black", hjust=1,vjust = 0.5, angle=90),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank(),
      legend.title = element_text(face="bold")
    )
  
  outFile <- file.path(outFolder, paste0("TAD", "_", plot_var, "_all_hicds_boxplot.", plotType))
  ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
}

######################################################### >>> plot # TADs by chromo

all_dt_chr <- aggregate(region ~ hicds + chromo, FUN=length, data = all_dt)

colnames(all_dt_chr)[colnames(all_dt_chr) == "region"] <- "chromo_nbrTADs"

all_vars <- c("chromo_nbrTADs")
plot_var=all_vars[1]
for(plot_var in all_vars) {
  outFile <- file.path(outFolder, paste0("all_hicds_dist_TAD_", plot_var, "_density_lattice.", plotType))
  do.call(plotType, list(outFile, height=myHeight*3, width=myWidth*3))
  # myplot <- densityplot( formula(paste0("~",plot_var)), groups = hicds, data = all_dt, #auto.key = TRUE,
  myplot <- densityplot( formula(paste0("~",plot_var, "| hicds")), data = all_dt_chr, #auto.key = TRUE,
                         # par.strip.text=list(cex=1), # width of the strip bar
                         par.strip.text = list(cex = 0.8, font = 4, col = "brown"),
                         layout = c(nCols, nRows), # column,row
                         scales=list(y=list(relation="free"),
                                     x=list(relation="free")
                         ),
                         main = paste0(plot_var))
  print(myplot)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

######################################################### >>> boxplot # TADs

all_dt_chr$datasetName <- gsub("_40kb", "", all_dt_chr$hicds)
stopifnot(is.numeric(all_dt_chr$chromo_nbrTADs))

mean_dt <- aggregate(chromo_nbrTADs ~ datasetName, FUN=mean, data=all_dt_chr)
mean_dt <- mean_dt[order(mean_dt$chromo_nbrTADs, decreasing = TRUE),]

all_dt_chr$datasetName <- factor(as.character(all_dt_chr$datasetName), levels=as.character(mean_dt$datasetName))

plot_vars <- c("chromo_nbrTADs")

plot_var=plot_vars[1]

for(plot_var in plot_vars) {
  
  p_var <-  ggplot(all_dt_chr, aes_string(x = "datasetName", y = plot_var)) + 
    geom_boxplot()+
    # coord_cartesian(expand = FALSE) +
    ggtitle(paste0("TAD ", plot_var))+
    scale_x_discrete(name="")+
    scale_y_continuous(name=paste0(plot_var),
                       breaks = scales::pretty_breaks(n = 20))+
    theme( # Increase size of axis lines
      strip.text = element_text(size = 12),
      # top, right, bottom and left
      # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"),
      strip.text.x = element_text(size = 10),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
      axis.text.x =  element_text(color="black", hjust=1,vjust = 0.5, angle=90),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(color="black", size=12),
      axis.title.x = element_text(color="black", size=12),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank(),
      legend.title = element_text(face="bold")
    )
  
  outFile <- file.path(outFolder, paste0("TAD", "_", plot_var, "_all_hicds_boxplot.", plotType))
  ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
}


all_dt_chr_size <- aggregate(size_log10 ~ datasetName + chromo, data = all_dt, FUN=mean)

colnames(all_dt_chr_size)[colnames(all_dt_chr_size) == "size_log10"] <- "chromo_meanSizeLog10"

merged_dt <- merge(all_dt_chr[, c("datasetName", "chromo", "chromo_nbrTADs"),],
                   all_dt_chr_size[, c("datasetName", "chromo", "chromo_meanSizeLog10"),],
                   by=c("datasetName", "chromo"),
                   all=TRUE
                   )

stopifnot(nrow(all_dt_chr) == nrow(all_dt_chr_size))
stopifnot(nrow(all_dt_chr) == nrow(merged_dt))

totDS <- length(unique(merged_dt$datasetName))

x_var <- "chromo_nbrTADs"
y_var <- "chromo_meanSizeLog10"

myx <- merged_dt[,paste0(x_var)]
myy <- merged_dt[,paste0(y_var)]

outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_all_hicds_densplot", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x=myx,
  y=myy,
  xlab = paste0(x_var),
  ylab = paste0(y_var),
  main = paste0(gsub("chromo_", "", y_var), " vs. ", gsub("chromo_", "", x_var))
)
mtext(side=3, text = paste0("all hicds - n=", totDS))
addCorr(x = myx, y = myy, bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


####################################################################################################################################### >>> collect # of pipeline genes & pipeline regions

all_dt_pipeline <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %dopar% {
    
    gene_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(gene_file))  
    
    region_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(region_file))  
    
    geneList <- get(load(gene_file))
    regionList <- get(load(region_file))
    
    
    data.frame(
      hicds= hicds,
      exprds = exprds,
      nPipelineGenes = length(geneList),
      nPipelineRegions = length(regionList),
      stringsAsFactors = FALSE
    )
  }
  exprds_dt
}
outFile <- file.path(outFolder, "all_dt_pipeline.Rdata")
save(all_dt_pipeline, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

all_dt_pipeline$dataset <- paste0(all_dt_pipeline$hicds, "\n", all_dt_pipeline$exprds)

p_var <-  ggplot(all_dt_pipeline, aes(x = dataset, y = nPipelineGenes)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("")+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0(""),
                     breaks = scales::pretty_breaks(n = 20))+
  # scale_fill_brewer(palette="YlOrRd")+
  # labs(fill  = "meanCorr type") +
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    strip.text.x = element_text(size = 10),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .3, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("nPipelineGenes_barplot.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
