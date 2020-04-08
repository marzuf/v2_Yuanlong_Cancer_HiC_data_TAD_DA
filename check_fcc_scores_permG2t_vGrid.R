

options(scipen=100)

# Rscript check_fcc_scores_permG2t_vGrid.R

script_name <- "check_fcc_scores_permG2t_vGrid.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildTable <- TRUE

require(flux)
require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)
registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("CHECK_FCC_SCORES_PERMG2T_VGRID")
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.4
plotType <- "svg"
myHeight <- 8
myWidth <- 10

# all_fcc_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern = "all_obs_prodSignedRatio.Rdata", recursive = TRUE, full.names = TRUE)
# all_fcc_obs_files <- all_fcc_obs_files[!grepl("_RANDOM", all_fcc_obs_files)]
# all_fcc_obs_files <- all_fcc_obs_files[!grepl("_PERMUT", all_fcc_obs_files)]
# all_rd_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern="all_obs_ratioDown.Rdata", recursive=TRUE, full.names = TRUE)
# all_rd_obs_files <- all_rd_obs_files[!grepl("_RANDOM", all_rd_obs_files)]
# all_rd_obs_files <- all_rd_obs_files[!grepl("_PERMUT", all_rd_obs_files)]
# all_negFC_obs_files <- list.files("OBS_TAD_NEGATIVE_FC", pattern="all_obs_negFC.Rdata", recursive=TRUE, full.names = TRUE)
# all_ratioFC_obs_files <- list.files("OBS_TAD_FC_RATIO", pattern="all_obs_ratioFC.Rdata", recursive=TRUE, full.names = TRUE)
# 
# stopifnot(length(all_rd_obs_files) == length(all_fcc_obs_files) )
# stopifnot(length(all_rd_obs_files) == length(all_negFC_obs_files) )
# stopifnot(length(all_rd_obs_files) == length(all_ratioFC_obs_files) )

pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_hicds <- all_hicds[!grepl("_RANDOM", all_hicds)]
all_hicds <- all_hicds[!grepl("_PERMUT", all_hicds)]

all_hicds <- all_hicds[!grepl("NCI-H460_", all_hicds)]

all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

ratio_vect <- seq(from=0, to=1, by=0.1)
ratio_fract_names <- paste0("ratio > ", ratio_vect[1:(length(ratio_vect)-1)], " and ratio <= ",ratio_vect[2:length(ratio_vect)])
ratio_fract_names <- paste0("ratio \u2208 ]", ratio_vect[1:(length(ratio_vect)-1)], ", ",ratio_vect[2:length(ratio_vect)], "]")
ratio_fract_names <- paste0("]", ratio_vect[1:(length(ratio_vect)-1)], ", ",ratio_vect[2:length(ratio_vect)], "]")
ratio_fract_names[ratio_fract_names == "]0, 0.1]"] <- "[0, 0.1]"

hicds = "ENCSR504OTV_transverse_colon_40kb"
exprds ="TCGAcoad_msi_mss"

keepPermut <- 1000

if(buildTable) {
  
  all_values_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    ds_values <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      

      # load("PERMG2T_TAD_FC_RATIO/ENCSR504OTV_transverse_colon_40kb/TCGAcoad_msi_mss/ratioFC_1000Permut_permDT.Rdata")
      
      ratioDown_file <- file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown/ratioDown_permDT.Rdata")
      
      if(!file.exists(ratioDown_file)) return(NULL)
      
      ratioFC_file <- file.path("PERMG2T_TAD_FC_RATIO", hicds, exprds, paste0("ratioFC_", keepPermut, "Permut_permDT.Rdata"))
      
      if(!file.exists(ratioFC_file)) return(NULL)
      
      cat("... ", hicds, " - " ,exprds, "\n")
      
      ratioDown_value_dt <- get(load(ratioDown_file))
      rD_value <- ratioDown_value_dt[,1:keepPermut]
      
      
      ratioFC_value <- get(load(ratioFC_file))
      
      stopifnot(dim(rD_value) == dim(ratioFC_value))
      
      stopifnot(rownames(rD_value) == rownames(ratioFC_value))
      
      
      data.frame(
        hicds = hicds,
        exprds = exprds,

        ratioFC = as.numeric(ratioFC_value),
        rD = as.numeric(rD_value),
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



hicds_toplot <- c("ENCSR504OTV_transverse_colon_40kb")
exprds_toplot <- c("TCGAcoad_msi_mss")


hicds_toplot <- c("ENCSR312KHQ_SK-MEL-5_40kb")
exprds_toplot <- c( "TCGAskcm_wt_mutCTNNB1")

init_all_values_dt <- all_values_dt

for(i in 1:length(hicds_toplot)){
  
  curr_hicds <- hicds_toplot[i]
  curr_exprds <- exprds_toplot[i]
  
  all_values_dt <- init_all_values_dt[order(all_values_dt$ratioDown_fract, all_values_dt$ratioFC_fract),]
  
  all_values_dt <- all_values_dt[all_values_dt$hicds == curr_hicds & all_values_dt$exprds == curr_exprds, ]
  
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
  
  
  stopifnot(nrow(nCount_ratioDown_ratioFC_fract_dt) > 0)
  
  grid_plot <- ggplot(data = nCount_ratioDown_ratioFC_fract_dt, aes(x=ratioDown_fract_lab, y=ratioFC_fract_lab, fill=ratioDown_ratioFC_nCout)) +
    ggtitle(paste0("# TADs by ratioDown and ratioFC"),   
            subtitle = paste0("", curr_hicds, " - ", curr_exprds))+
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
  
  outFile <- file.path(outFolder, paste0(curr_hicds, "_", curr_exprds, "_ratioDown_ratioFC_gridPlot.", plotType))
  ggsave(grid_plot, file = outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
  
}





#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))












