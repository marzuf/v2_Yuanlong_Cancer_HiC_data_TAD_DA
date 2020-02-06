

options(scipen=100)

# Rscript check_fcc_scores_permut.R

script_name <- "check_fcc_scores_permut.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildTable <- FALSE

require(flux)
require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)
registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("CHECK_FCC_SCORES_PERMUT")
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.4
plotType <- "png"
myHeight <- 400
myWidth <- 400
myWidthLeg <- 500

all_fcc_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern = "all_obs_prodSignedRatio.Rdata", recursive = TRUE, full.names = TRUE)
all_rd_obs_files <- list.files("PIPELINE/OUTPUT_FOLDER", pattern="all_obs_ratioDown.Rdata", recursive=TRUE, full.names = TRUE)
all_negFC_obs_files <- list.files("OBS_TAD_NEGATIVE_FC", pattern="all_obs_negFC.Rdata", recursive=TRUE, full.names = TRUE)
all_ratioFC_obs_files <- list.files("OBS_TAD_FC_RATIO", pattern="all_obs_ratioFC.Rdata", recursive=TRUE, full.names = TRUE)

stopifnot(length(all_rd_obs_files) == length(all_fcc_obs_files) )
stopifnot(length(all_rd_obs_files) == length(all_negFC_obs_files) )
stopifnot(length(all_rd_obs_files) == length(all_ratioFC_obs_files) )

pipFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
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
        
        # left FCC
        tad_entrez_left <- permut_data[[paste0(rd_tad)]][["genes_left"]]
        stopifnot(tad_entrez_left %in% geneList)
        tad_entrez_left_de <- names(geneList)[geneList %in% tad_entrez_left]
        stopifnot(!duplicated(tad_entrez_left_de))
        tad_entrez_left_de <- unique(tad_entrez_left_de)
        stopifnot(tad_entrez_left_de %in% de_dt$genes)
        all_tad_left_fc <- de_dt$logFC[de_dt$genes %in% tad_entrez_left_de]
        
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
} else {
  # outFile="CHECK_FCC_SCORES_PERMUT/all_permut_ratios.Rdata"
  outFile <- file.path(outFolder, "all_permut_ratios.Rdata")
  all_permut_ratios <- get(load(outFile))
}



all_sfx <- c("all", "right", "left")
all_sfx <- c("right", "left")

rbPal <- colorRampPalette(c('red','blue'))

for(sfx in all_sfx) {
  
  
  all_permut_ratios[,paste0("dotCols_", sfx)] <- rev(rbPal(10))[as.numeric(cut(all_permut_ratios[,paste0("fcc_", sfx)],breaks = 10))]
  
  nDS <- length(unique(file.path(all_permut_ratios$hicds, all_permut_ratios$exprds)))
  
  outFile <- file.path(outFolder, paste0("ratioDown_ratioFC_colFCC_", sfx, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidthLeg))
  
  par(xpd = T, mar = par()$mar + c(0,0,0,6))
  
  plot(
    x = all_permut_ratios[,paste0("ratioFC_", sfx)],
    y = all_permut_ratios[,paste0("ratioDown_", sfx)],
    xlab=paste0("ratioNegFC_", sfx),
    ylab=paste0("ratioDown_", sfx),
    col = all_permut_ratios[,paste0("dotCols_", sfx)],
    pch = 16,
    cex = 0.7,
    cex.lab = plotCex,
    cex.axis = plotCex,
    main = paste0("all datasets (permut ", sfx, ")")
  )
  mtext(side=3, text = paste0("# DS = ", nDS, "; # TADs = ", nrow(all_permut_ratios)))
  
  
  # legend("bottomright",
  legend(1.1,1,
         title="FCC score",
         legend=rev(levels(cut(all_permut_ratios[,paste0("fcc_", sfx)],breaks = 10))),
         col =rbPal(10),
         pch=20,
         cex = 0.8,
         ncol=1,
         horiz = F,
         bty="n")
  
  
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0("ratioNegFC_FCCscore_", sfx, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = all_permut_ratios[,paste0("ratioFC_", sfx)],
    y = all_permut_ratios[,paste0("fcc_", sfx)],
    xlab=paste0("ratioNegFC_", sfx),
    ylab=paste0("FCC score_", sfx),
    # col = all_permut_ratios$dotCols,
    pch = 16,
    cex = 0.7,
    cex.lab = plotCex,
    cex.axis = plotCex,
    main = paste0("all datasets (permut ", sfx, ")")
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0("ratioDown_FCCscore_", sfx, "_densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = all_permut_ratios[,paste0("ratioDown_", sfx)],
    y = all_permut_ratios[,paste0("fcc_", sfx)],
    xlab="ratioDown",
    ylab="FCC score",
    # col = all_permut_ratios$dotCols,
    pch = 16,
    cex = 0.7,
    cex.lab = plotCex,
    cex.axis = plotCex,
    main = paste0("all datasets (permut ", sfx, ")")
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  tmp_dt <- all_permut_ratios[,c(paste0("ratioFC_", sfx), paste0("ratioDown_", sfx), paste0("fcc_", sfx))]
  tmp_dt <- na.omit(tmp_dt)
  stopifnot( (2*tmp_dt[,paste0("ratioFC_", sfx)] - 1) * (2* tmp_dt[,paste0("ratioDown_", sfx)]- 1) == tmp_dt[,paste0("fcc_", sfx)]) 
  
  
  
  
  
  
  
  sub_permut_ratios <- all_permut_ratios[!is.na(all_permut_ratios[,paste0("fcc_", sfx)]),]
  
  
  fcc_fract <- seq(from=-1, to=1, by=0.25)
  # fcc_fract_names <- paste0("FCC > ", fcc_fract[1:(length(fcc_fract)-1)], " and FCC <= ",fcc_fract[2:length(fcc_fract)])
  fcc_fract_names <- paste0("FCC \u2208 ]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
  fcc_fract_names <- paste0("]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
  fcc_fract_names[fcc_fract_names == "]-1, -0.75]"] <- "[-1, -0.75]"
  
  
  
  fcc_scoreFract <- sapply(1:nrow(sub_permut_ratios) ,function(x) which(hist(sub_permut_ratios[x,paste0("fcc_", sfx)], breaks=fcc_fract, plot=F)$counts == 1))
  
  
  save(sub_permut_ratios,file= "sub_permut_ratios.Rdata", version=2)
  save(fcc_scoreFract, file ="fcc_scoreFract.Rdata", version=2)
  
  
  stopifnot(length(fcc_scoreFract) == nrow(sub_permut_ratios))
  
  require(ggsci)
  ggsci_pal <- "lancet"
  ggsci_subpal <- ""
  myPals <-  eval(parse(text=paste0("pal_", ggsci_pal, "(", ggsci_subpal, ")")))(length(unique(fcc_fract_names)))
  myPals <- rev(myPals)
  
  outFile <- file.path(outFolder, paste0("ratioDown_ratioFC_colFCCfract_", sfx, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidthLeg))
  
  par(xpd = T, mar = par()$mar + c(0,0,0,6))
  
  save(myPals, file= "myPals.Rdata", version=2)
  
  
  
  plot(
    x = sub_permut_ratios[,paste0("ratioFC_", sfx)],
    y = sub_permut_ratios[,paste0("ratioDown_", sfx)],
    xlab=paste0("ratioNegFC_", sfx),
    ylab=paste0("ratioDown_", sfx),
    col = myPals[fcc_scoreFract],
    pch = 16,
    cex = 0.7,
    cex.lab = plotCex,
    cex.axis = plotCex,
    main = paste0("all datasets (permut ", sfx, ")")
  )
  mtext(side=3, text = paste0("# DS = ", nDS, "; # TADs = ", nrow(sub_permut_ratios)))
  
  # legend("bottomright",
  legend(1.1,1,
         title="FCC score",
         legend= rev(fcc_fract_names),
         col =rev(myPals),
         pch=20,
         cex = 0.8,
         ncol=1,
         horiz = F,
         bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
}


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))












