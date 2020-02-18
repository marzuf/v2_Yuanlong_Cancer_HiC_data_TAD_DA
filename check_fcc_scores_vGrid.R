

options(scipen=100)

# Rscript check_fcc_scores_vGrid.R

script_name <- "check_fcc_scores_vGrid.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildTable <- F

require(flux)
require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)
registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("CHECK_FCC_SCORES_VGRID")
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.4
plotType <- "svg"
myHeight <- 8
myWidth <- 10

all_fcc_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern = "all_obs_prodSignedRatio.Rdata", recursive = TRUE, full.names = TRUE)
all_fcc_obs_files <- all_fcc_obs_files[!grepl("_RANDOM", all_fcc_obs_files)]
all_fcc_obs_files <- all_fcc_obs_files[!grepl("_PERMUT", all_fcc_obs_files)]
all_rd_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern="all_obs_ratioDown.Rdata", recursive=TRUE, full.names = TRUE)
all_rd_obs_files <- all_rd_obs_files[!grepl("_RANDOM", all_rd_obs_files)]
all_rd_obs_files <- all_rd_obs_files[!grepl("_PERMUT", all_rd_obs_files)]
all_negFC_obs_files <- list.files("OBS_TAD_NEGATIVE_FC", pattern="all_obs_negFC.Rdata", recursive=TRUE, full.names = TRUE)
all_ratioFC_obs_files <- list.files("OBS_TAD_FC_RATIO", pattern="all_obs_ratioFC.Rdata", recursive=TRUE, full.names = TRUE)

stopifnot(length(all_rd_obs_files) == length(all_fcc_obs_files) )
stopifnot(length(all_rd_obs_files) == length(all_negFC_obs_files) )
stopifnot(length(all_rd_obs_files) == length(all_ratioFC_obs_files) )

pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_hicds <- all_hicds[!grepl("_RANDOM", all_hicds)]
all_hicds <- all_hicds[!grepl("_PERMUT", all_hicds)]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

ratio_vect <- seq(from=0, to=1, by=0.1)
ratio_fract_names <- paste0("ratio > ", ratio_vect[1:(length(ratio_vect)-1)], " and ratio <= ",ratio_vect[2:length(ratio_vect)])
ratio_fract_names <- paste0("ratio \u2208 ]", ratio_vect[1:(length(ratio_vect)-1)], ", ",ratio_vect[2:length(ratio_vect)], "]")
ratio_fract_names <- paste0("]", ratio_vect[1:(length(ratio_vect)-1)], ", ",ratio_vect[2:length(ratio_vect)], "]")
ratio_fract_names[ratio_fract_names == "]0, 0.1]"] <- "[0, 0.1]"


if(buildTable) {
  
  all_values_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    ds_values <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      fcc_value <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata")))
      
      negFC_value <- get(load(file.path("OBS_TAD_NEGATIVE_FC", hicds, exprds, "all_obs_negFC.Rdata")))
      
      rd_value <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown/all_obs_ratioDown.Rdata")))
      
      ratioFC_value <- get(load(file.path("OBS_TAD_FC_RATIO", hicds, exprds, "all_obs_ratioFC.Rdata")))
      
      stopifnot(setequal(names(fcc_value), names(negFC_value)))
      stopifnot(setequal(names(fcc_value), names(rd_value)))
      stopifnot(setequal(names(fcc_value), names(ratioFC_value)))
      stopifnot(length(fcc_value) == length(negFC_value))
      stopifnot(length(fcc_value) == length(rd_value))
      stopifnot(length(fcc_value) == length(ratioFC_value))
      
      all_tads <- names(fcc_value)
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        FCCscore = fcc_value[all_tads],
        negFC = negFC_value[all_tads],
        ratioFC = ratioFC_value[all_tads],
        rD = rd_value[all_tads],
        stringsAsFactors = FALSE
      )
    } # end foreach iterating exprds
    ds_values
  } # end foreach iterating hicds
  
  
  all_values_dt$ratioFC_fract <- sapply(all_values_dt$ratioFC, function(x) which(hist(x, breaks=ratio_vect, plot=F)$counts == 1))
  all_values_dt$ratioDown_fract <- sapply(all_values_dt$rD, function(x) which(hist(x, breaks=ratio_vect, plot=F)$counts == 1))
  all_values_dt$ratioFC_fractName <- ratio_fract_names[all_values_dt$ratioFC_fract]
  all_values_dt$ratioDown_fractName <- ratio_fract_names[all_values_dt$ratioDown_fract]
  
  outFile <- file.path(outFolder, "all_values_dt.Rdata")
  save(all_values_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  stopifnot(!is.na(all_values_dt$ratioFC_fractName))
  stopifnot(!is.na(all_values_dt$ratioDown_fractName))
  
  outFile <- file.path(outFolder, "all_values_dt.Rdata")
  save(all_values_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  # outFile="PIPELINE/OUTPUT_FOLDER/all_values_dt.Rdata"
  # outFile="CHECK_FCC_SCORES_VGRID/all_values_dt.Rdata"
  outFile <- file.path(outFolder, "all_values_dt.Rdata")
  all_values_dt <- get(load(outFile))
}

all_values_dt <- all_values_dt[order(all_values_dt$ratioDown_fract, all_values_dt$ratioFC_fract),]

all_values_dt$ratioDown_ratioFC_fractNameInteraction <- interaction(all_values_dt$ratioDown_fractName, all_values_dt$ratioFC_fractName)
all_values_dt$ratioDown_ratioFC_fractInteraction <- interaction(all_values_dt$ratioDown_fract, all_values_dt$ratioFC_fract)


nCount_ratioDown_ratioFC_fract_dt <- data.frame(
  ratioDown_ratioFC_fract = names(table(all_values_dt$ratioDown_ratioFC_fractInteraction)),
  ratioDown_ratioFC_nCout = as.numeric(table(all_values_dt$ratioDown_ratioFC_fractInteraction)),
  stringsAsFactors = FALSE
)

nCount_ratioDown_ratioFC_fract_dt$ratioDown_fract <- gsub("(.+)\\..+", "\\1", nCount_ratioDown_ratioFC_fract_dt$ratioDown_ratioFC_fract)
nCount_ratioDown_ratioFC_fract_dt$ratioFC_fract <- gsub(".+\\.(.+)", "\\1", nCount_ratioDown_ratioFC_fract_dt$ratioDown_ratioFC_fract)
nCount_ratioDown_ratioFC_fract_dt$ratioDown_ratioFC_fract <- NULL

nCount_ratioDown_ratioFC_fract_dt$ratioDown_fract_lab <- ratio_fract_names[as.numeric(nCount_ratioDown_ratioFC_fract_dt$ratioDown_fract)]
nCount_ratioDown_ratioFC_fract_dt$ratioFC_fract_lab <- ratio_fract_names[as.numeric(nCount_ratioDown_ratioFC_fract_dt$ratioFC_fract)]


nCount_ratioDown_ratioFC_fract_dt$ratioDown_fract_lab <- factor(nCount_ratioDown_ratioFC_fract_dt$ratioDown_fract_lab, levels = ratio_fract_names)
nCount_ratioDown_ratioFC_fract_dt$ratioFC_fract_lab <- factor(nCount_ratioDown_ratioFC_fract_dt$ratioFC_fract_lab, levels = ratio_fract_names)

stopifnot(!is.na(nCount_ratioDown_ratioFC_fract_dt$ratioDown_fract_lab))
stopifnot(!is.na(nCount_ratioDown_ratioFC_fract_dt$ratioFC_fract_lab))


grid_plot <- ggplot(data = nCount_ratioDown_ratioFC_fract_dt, aes(x=ratioDown_fract_lab, y=ratioFC_fract_lab, fill=ratioDown_ratioFC_nCout)) + 
  ggtitle(paste0("# TADs by ratioDown and ratioFC"),   subtitle = paste0("obs. data - all datasets (n=", length(unique(file.path(all_values_dt$hicds, all_values_dt$exprds))), ")"))+
  scale_x_discrete(name="ratioDown")  + 
  scale_y_discrete(name="ratioFC")  + 
  geom_tile() +
  labs(fill = "# TADs")+
  scale_fill_gradient( trans = 'log', na.value = "white" )  +
  theme(
    axis.text.x = element_text(colour = "black", size=12),
    axis.text.y= element_text(colour = "black", size=12),
    axis.title.x = element_text(colour = "black", size=14, face="bold"),
    axis.title.y = element_text(colour = "black", size=14, face="bold"),
    plot.title = element_text(hjust=0.5, size=16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
    panel.background = element_rect(fill = "transparent")
    # legend.background =  element_rect()
  )

outFile <- file.path(outFolder, paste0("ratioDown_ratioFC_gridPlot.", plotType))
ggsave(grid_plot, file = outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))


grid_plot_nbr <- grid_plot + geom_text(aes(label=ratioDown_ratioFC_nCout))
outFile <- file.path(outFolder, paste0("ratioDown_ratioFC_gridPlot_withNbr.", plotType))
ggsave(grid_plot_nbr, file = outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))




mat_dt <- acast(nCount_ratioDown_ratioFC_fract_dt, ratioDown_fract~ratioFC_fract, value.var="ratioDown_ratioFC_nCout")
mat_dt <- mat_dt[order(as.numeric(rownames(mat_dt))), order(as.numeric(rownames(mat_dt)))]
rownames(mat_dt) <- paste0("ratioDown_", rownames(mat_dt))
colnames(mat_dt) <- paste0("ratioFC_", colnames(mat_dt))

outFile <- file.path(outFolder, "mat_dt.Rdata")
save(mat_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))





stopifnot( (2*all_values_dt$ratioFC - 1) * (2* all_values_dt$rD- 1) == all_values_dt$FCCscore)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))












