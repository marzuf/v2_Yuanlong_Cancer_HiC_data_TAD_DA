

options(scipen=100)

# Rscript check_fcc_scores_permut_vGrid.R

script_name <- "check_fcc_scores_permut_vGrid.R"

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

outFolder <- file.path("CHECK_FCC_SCORES_PERMUT_VGRID")
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.4
plotType <- "svg"
myHeight <- 8
myWidth <- 10
all_sfx <- c("all", "right", "left")

ratio_vect <- seq(from=0, to=1, by=0.1)
ratio_fract_names <- paste0("ratio > ", ratio_vect[1:(length(ratio_vect)-1)], " and ratio <= ",ratio_vect[2:length(ratio_vect)])
ratio_fract_names <- paste0("ratio \u2208 ]", ratio_vect[1:(length(ratio_vect)-1)], ", ",ratio_vect[2:length(ratio_vect)], "]")
ratio_fract_names <- paste0("]", ratio_vect[1:(length(ratio_vect)-1)], ", ",ratio_vect[2:length(ratio_vect)], "]")
ratio_fract_names[ratio_fract_names == "]0, 0.1]"] <- "[0, 0.1]"


pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_hicds <- all_hicds[!grepl("_RANDOM", all_hicds)]
all_hicds <- all_hicds[!grepl("_PERMUT", all_hicds)]

all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
permut_script <- "5sameNbr_runPermutationsCorr"


get_fcc <- function(fc_vect) {
  (2* sum(fc_vect < 0)/length(fc_vect) -1) *  (2* sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect)) -1)
}

get_ratioDown <- function(fc_vect) {
  sum(fc_vect < 0)/length(fc_vect) 
}

get_ratioFC <- function(fc_vect) {
  sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect))
}




if(buildTable) {
  all_permut_ratios <- foreach(hicds = all_hicds, .combine='rbind') %do%{
    exprds_ratios <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat(paste0("... start ", hicds, " - ", exprds, "\n"))
      
      permut_data <- get(load(file.path(pipFolder, hicds, exprds, permut_script, "sample_around_TADs_sameNbr.Rdata")))
      
      de_dt <- get(load(file.path(pipFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")))
      
      geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
      all(names(geneList) %in% de_dt$genes)
      
      rd_tad = names(permut_data)[1]
      
      ds_all_permut <- foreach(rd_tad = names(permut_data), .combine='rbind') %dopar% {
        
        # right FCC
        tad_entrez_right <- permut_data[[paste0(rd_tad)]][["genes_right"]]
        stopifnot(tad_entrez_right %in% geneList)
        tad_entrez_right_de <- names(geneList)[geneList %in% tad_entrez_right]
        stopifnot(!duplicated(tad_entrez_right_de))
        tad_entrez_right_de <- unique(tad_entrez_right_de)
        stopifnot(tad_entrez_right_de %in% de_dt$genes)
        all_tad_right_fc <- de_dt$logFC[de_dt$genes %in% tad_entrez_right_de]
        fcc_right <- get_fcc(all_tad_right_fc)
        
        all_tad_right_fc_txt <- paste0(all_tad_right_fc, collapse=",")
        
        # left FCC
        tad_entrez_left <- permut_data[[paste0(rd_tad)]][["genes_left"]]
        stopifnot(tad_entrez_left %in% geneList)
        tad_entrez_left_de <- names(geneList)[geneList %in% tad_entrez_left]
        stopifnot(!duplicated(tad_entrez_left_de))
        tad_entrez_left_de <- unique(tad_entrez_left_de)
        stopifnot(tad_entrez_left_de %in% de_dt$genes)
        all_tad_left_fc <- de_dt$logFC[de_dt$genes %in% tad_entrez_left_de]
        
        all_tad_left_fc_txt <- paste0(all_tad_left_fc, collapse=",")
        
        fcc_left <- get_fcc(all_tad_left_fc)
        
        fcc_mean <- mean(c(fcc_right, fcc_left), na.rm=TRUE)
        
        ratioDown_right <- get_ratioDown(all_tad_right_fc)
        ratioFC_right <- get_ratioFC(all_tad_right_fc)
        
        ratioDown_left <- get_ratioDown(all_tad_left_fc)
        ratioFC_left <- get_ratioFC(all_tad_left_fc)
        
        
        # all FCC
        tad_entrez <- union(permut_data[[paste0(rd_tad)]][["genes_right"]], permut_data[[paste0(rd_tad)]][["genes_left"]])
        stopifnot(tad_entrez %in% geneList)
        tad_entrez_de <- names(geneList)[geneList %in% tad_entrez]
        stopifnot(!duplicated(tad_entrez_de))
        tad_entrez_de <- unique(tad_entrez_de)
        stopifnot(tad_entrez_de %in% de_dt$genes)
        all_tad_fc <- de_dt$logFC[de_dt$genes %in% tad_entrez_de]
        fcc_all <- get_fcc(all_tad_fc)
        
        ratioDown_all <- get_ratioDown(all_tad_fc)
        ratioFC_all <- get_ratioFC(all_tad_fc)
        
        data.frame(
          
          hicds = hicds, 
          exprds = exprds,
          
          region = rd_tad,
          
          fcc_all = fcc_all,
          fcc_right = fcc_right,
          fcc_left = fcc_left,
          
          ratioDown_right=ratioDown_right,
          ratioFC_right=ratioFC_right,
          
          ratioDown_left=ratioDown_left,
          ratioFC_left=ratioFC_left,
          
          ratioDown_all=ratioDown_all,
          ratioFC_all=ratioFC_all,
          
          all_tad_left_fc = all_tad_left_fc_txt,
          all_tad_right_fc = all_tad_right_fc_txt,
          
          stringsAsFactors = FALSE
        )
      } # end-foreach TAD
      
      ds_all_permut
      
    } # end foreach iterating exprds
    exprds_ratios
  } # end foreach iterating hicds
  
  
  outFile <- file.path(outFolder, "all_permut_ratios.Rdata")
  save(all_permut_ratios, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  for(side in all_sfx) {
    all_permut_ratios[,paste0("ratioFC_", side, "_fract")] <- sapply(all_permut_ratios[,paste0("ratioFC_", side)], function(x) {
      if(is.na(x)) return(NA)
      which(hist(x, breaks=ratio_vect, plot=F)$counts == 1)
    })
    all_permut_ratios[,paste0("ratioDown_", side, "_fract")] <- sapply(all_permut_ratios[,paste0("ratioDown_", side)], function(x) {
      if(is.na(x)) return(NA)
      which(hist(x, breaks=ratio_vect, plot=F)$counts == 1)
    })
    all_permut_ratios[,paste0("ratioFC_", side, "_fractName")] <- ratio_fract_names[all_permut_ratios[,paste0("ratioFC_", side, "_fract")]]
    all_permut_ratios[,paste0("ratioDown_", side, "_fractName")] <- ratio_fract_names[all_permut_ratios[,paste0("ratioDown_", side, "_fract")]]
    
  }
  
  outFile <- file.path(outFolder, "all_permut_ratios.Rdata")
  save(all_permut_ratios, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  # stopifnot(!is.na(all_permut_ratios[,paste0("ratioFC_", side, "_fractName")]))
  # stopifnot(!is.na(all_permut_ratios[,paste0("ratioDown_", side, "_fractName")]))
  
  outFile <- file.path(outFolder, "all_permut_ratios.Rdata")
  save(all_permut_ratios, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  # stop("ok")
  

} else {
  # outFile="CHECK_FCC_SCORES_PERMUT/all_permut_ratios.Rdata"
  outFile <- file.path(outFolder, "all_permut_ratios.Rdata")
  all_permut_ratios <- get(load(outFile))
}


for(side in all_sfx) {
  
  all_values_dt <- all_permut_ratios[,c( paste0("ratioFC_", side), paste0("ratioDown_", side) ,
                                         paste0("ratioDown_", side, "_fractName"), paste0("ratioDown_", side, "_fract"),
                                         paste0("ratioFC_", side, "_fractName"), paste0("ratioFC_", side, "_fract"))]
  
  all_values_dt <- all_values_dt[order(all_values_dt[,paste0("ratioDown_", side, "_fract")], all_values_dt[,paste0("ratioDown_", side, "_fract")]),]
  
  all_values_dt$ratioDown_ratioFC_fractNameInteraction <- interaction(all_values_dt[,paste0("ratioDown_", side, "_fractName")], all_values_dt[,paste0("ratioFC_", side, "_fractName")])
  all_values_dt$ratioDown_ratioFC_fractInteraction <- interaction(all_values_dt[,paste0("ratioDown_", side, "_fract")], all_values_dt[,paste0("ratioFC_", side, "_fract")])
  
  
  nCount_ratioDown_ratioFC_fract_dt <- data.frame(
    ratioDown_ratioFC_fract = names(table(all_values_dt$ratioDown_ratioFC_fractInteraction)),
    ratioDown_ratioFC_nCout = as.numeric(table(all_values_dt$ratioDown_ratioFC_fractInteraction)),
    stringsAsFactors = FALSE
  )
  
  nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioDown_", side, "_fract")] <- gsub("(.+)\\..+", "\\1", nCount_ratioDown_ratioFC_fract_dt$ratioDown_ratioFC_fract)
  nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioFC_", side, "_fract")] <- gsub(".+\\.(.+)", "\\1", nCount_ratioDown_ratioFC_fract_dt$ratioDown_ratioFC_fract)
  nCount_ratioDown_ratioFC_fract_dt$ratioDown_ratioFC_fract <- NULL
  
  nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioDown_", side, "_fract_lab")] <- ratio_fract_names[as.numeric(nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioDown_", side, "_fract")])]
  nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioFC_", side, "_fract_lab")] <- ratio_fract_names[as.numeric(nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioFC_", side, "_fract")])]
  
  
  nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioDown_", side, "_fract_lab")] <- factor(nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioDown_", side, "_fract_lab")], levels = ratio_fract_names)
  nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioFC_", side, "_fract_lab")] <- factor(nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioFC_", side, "_fract_lab")], levels = ratio_fract_names)
  
  stopifnot(!is.na(nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioDown_", side, "_fract_lab")]))
  stopifnot(!is.na(nCount_ratioDown_ratioFC_fract_dt[,paste0("ratioFC_", side, "_fract_lab")]))
  
  
  grid_plot <- ggplot(data = nCount_ratioDown_ratioFC_fract_dt, aes_string(x=paste0("ratioDown_", side, "_fract_lab"), 
                                                                           y=paste0("ratioFC_", side, "_fract_lab"), 
                                                                           fill=paste0("ratioDown_ratioFC_nCout"))) + 
    
    
    ggtitle(paste0("# TADs by ratioDown and ratioFC"),   
            subtitle = paste0(side, " permut. data - all datasets (n=", length(unique(file.path(all_values_dt$hicds, all_values_dt$exprds))), ")"))+
    
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
  
  outFile <- file.path(outFolder, paste0("ratioDown_ratioFC_", side, "permut_gridPlot.", plotType))
  ggsave(grid_plot, file = outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  grid_plot_nbr <- grid_plot + geom_text(aes(label=ratioDown_ratioFC_nCout))
  outFile <- file.path(outFolder, paste0("ratioDown_ratioFC_", side, "permut_gridPlot_withNbr.", plotType))
  ggsave(grid_plot_nbr, file = outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  

  mat_dt <- acast(nCount_ratioDown_ratioFC_fract_dt, as.formula(paste0("ratioDown_",side, "_fract~ratioFC_", side, "_fract")), value.var="ratioDown_ratioFC_nCout")
  mat_dt <- mat_dt[order(as.numeric(rownames(mat_dt))), order(as.numeric(rownames(mat_dt)))]
  rownames(mat_dt) <- paste0("ratioDown_", rownames(mat_dt))
  colnames(mat_dt) <- paste0("ratioFC_", colnames(mat_dt))
  
  outFile <- file.path(outFolder, paste0(side, "permut_mat_dt.Rdata"))
  save(mat_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))

}


# check some bins of the grid
# left permut: RatioDown [0-0.1]; ratioFC[0.5-0.6]

dt <- get(load(file.path("CHECK_FCC_SCORES_PERMUT_VGRID", "all_permut_ratios.Rdata")))
dt <- na.omit(dt)
dt[dt$ratioDown_left_fractName == "[0, 0.1]" & dt$ratioFC_left_fractName == "]0.5, 0.6]", c("hicds", "exprds", "region", "all_tad_left_fc")]
fc <- dt[dt$ratioDown_left_fractName == "[0, 0.1]" & dt$ratioFC_left_fractName == "]0.5, 0.6]", c("all_tad_left_fc")]

dt[dt$ratioDown_left_fractName == "]0.8, 0.9]" & dt$ratioFC_left_fractName == "]0.1, 0.2]", c("hicds", "exprds", "region", "all_tad_left_fc")]
fc <- dt[dt$ratioDown_left_fractName == "]0.8, 0.9]" & dt$ratioFC_left_fractName == "]0.1, 0.2]", c("all_tad_left_fc")][1]
round(as.numeric(unlist(strsplit(fc, ","))),3)


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))












