options(scipen=100)

# Rscript actual_TAD_maxSize.R

script_name <- "actual_TAD_maxSize.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)

registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")



plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 8


script0_name <- "0_prepGeneData"

setDir <- "/media/electron"
setDir <- ""


mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- "ACTUAL_TAD_MAXSIZE"
dir.create(outFolder, recursive = TRUE)

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

minTADsize <- 3

hicds = all_hicds[1]

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    size_file <- file.path(pipFolder, hicds, exprds, script0_name, "rangeTADgenes.Rdata")
    stopifnot(file.exists(size_file))
    size_info <- get(load(size_file))
    stopifnot(size_info[1] == minTADsize)
    data.frame(hicds=hicds, exprds=exprds,maxSize=size_info[2], stringsAsFactors = FALSE)
  } # end-foreach iterating over exprds
  exprds_dt
} # end-foreach iterating over hicds

all_dt$dataset <- paste0(all_dt$hicds, "\n", all_dt$exprds)
all_dt <- all_dt[order(all_dt$maxSize, decreasing=TRUE),]
all_dt$dataset <- factor(all_dt$dataset, levels=as.character(all_dt$dataset))

all_dt$exprds_type_col <- all_cols[all_cmps[all_dt$exprds]]

p_var <-  ggplot(all_dt, aes(x = dataset, y = maxSize)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE, ylim=c(minTADsize, max(all_dt$maxSize))) +
  ggtitle("Max TAD size (# of genes)")+
  scale_x_discrete(name="")+
  scale_y_continuous(name=paste0("max # genes"),
                     breaks = scales::pretty_breaks(n = 10))+
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
    axis.text.x = element_text(color=all_dt$exprds_type_col, hjust=1,vjust = 0.5, size=8, angle=90),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_text(color="black", size=12),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("all_ds_TAD_max_size_barplot.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



