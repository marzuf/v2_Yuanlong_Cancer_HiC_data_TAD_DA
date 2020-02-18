########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript cmp_permutTADs_obsTADs.R\n"))

# Rscript cmp_permutTADs_obsTADs.R

script_name <- "cmp_permutTADs_obsTADs.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)


plotType <- "svg"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4
myHeightGG <- 6
myWidthGG <- 7.5

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- "CMP_PERMUTTADS_OBSTADS"
dir.create(outFolder, recursive = TRUE)

hicds="ENCSR489OCU_NCI-H460_40kb"

mainFolder <- "."

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]

all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))


rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T" )


buildData <- FALSE

if(buildData){
  
  hicds=all_hicds[1]
  cmp_obs_permut_TADs_dt <- foreach(hicds=all_hicds, .combine='rbind') %dopar% {
    
    cat(paste0("... start - ", hicds, "\n"))
    
    g2t_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_file, col.names=c("entrezID", "chromo", "start", "end", "region"), header=FALSE, stringsAsFactors = FALSE)
    g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
    obs_g2t <- by(g2t_dt, g2t_dt$region, function(x) sort(x$entrezID))
    
    rd_patt = rd_patterns[1]
    hicds_perm_dt <-foreach(rd_patt = rd_patterns, .combine='rbind') %dopar% {
      
      cat(paste0("... start - ", hicds, " - ", rd_patt, "\n"))
      
      rd_hicds <- file.path(dirname(hicds), gsub("_40kb", paste0("_", rd_patt, "_40kb"), basename(hicds)))
      rd_g2t_file <- file.path(mainFolder, rd_hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(rd_g2t_file))
      rd_g2t_dt <- read.delim(rd_g2t_file, col.names=c("entrezID", "chromo", "start", "end", "region"), header=FALSE, stringsAsFactors = FALSE)
      
      rd_g2t_dt <- rd_g2t_dt[grepl("_TAD", rd_g2t_dt$region),]
      rd_g2t_dt$entrezID <- as.character(rd_g2t_dt$entrezID)
      
      rd_g2t <- by(rd_g2t_dt, rd_g2t_dt$region, function(x) sort(x$entrezID))
      
      rd_obs_match <- rd_g2t  %in% obs_g2t
      
      stopifnot(sum(rd_obs_match) == sum(!is.na(match(rd_g2t, obs_g2t))))
      
      data.frame(
        hicds = hicds,
        randomData= rd_patt,
        nObsTADs = length(obs_g2t),
        nRdTADs = length(rd_g2t),
        nSameTADs = sum(rd_obs_match),
        stringsAsFactors = FALSE
      )
    }
    hicds_perm_dt
  }
  outFile <- file.path(outFolder, "cmp_obs_permut_TADs_dt.Rdata")
  save(cmp_obs_permut_TADs_dt, file = outFile, version=2)
  cat(paste0("... written: ",outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "cmp_obs_permut_TADs_dt.Rdata")
  cmp_obs_permut_TADs_dt <- get(load(outFile))
}
stopifnot(cmp_obs_permut_TADs_dt$nSameTADs <= cmp_obs_permut_TADs_dt$nObsTADs)
cmp_obs_permut_TADs_dt$ratioSameTAD <- cmp_obs_permut_TADs_dt$nSameTADs/cmp_obs_permut_TADs_dt$nRdTADs


outFile <- file.path(outFolder, paste0("allDS_ratioSameTAD_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(cmp_obs_permut_TADs_dt$ratioSameTAD, cmp_obs_permut_TADs_dt$randomData),
  plotTit = paste0("all DS - n=", length(unique(cmp_obs_permut_TADs_dt$hicds)), " - ratio same TAD"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("allDS_ratioSameTAD_density_noG2T.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(cmp_obs_permut_TADs_dt$ratioSameTAD[cmp_obs_permut_TADs_dt$randomData!="PERMUTG2T"], cmp_obs_permut_TADs_dt$randomData[cmp_obs_permut_TADs_dt$randomData!="PERMUTG2T"]),
  plotTit = paste0("all DS - n=", length(unique(cmp_obs_permut_TADs_dt$hicds)), " - ratio same TAD"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


box_meanFC <- ggboxplot(data= cmp_obs_permut_TADs_dt, x="randomData", y= "ratioSameTAD", xlab="") +
  ggtitle("ratio same TAD", subtitle = paste0("all ds - n=", length(unique(cmp_obs_permut_TADs_dt$hicds))))+
  # facet_grid(~hicds_lab, switch="x") + 
  scale_y_continuous(name=paste0("ratio same TAD"),
                     breaks = scales::pretty_breaks(n = 10))+
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"))

outFile <- file.path(outFolder, paste0("allDS_ratioSameTAD_boxplot.", plotType))
ggsave(box_meanFC, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))








#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))


