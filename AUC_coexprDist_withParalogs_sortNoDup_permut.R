startTime <- Sys.time()
cat(paste0("> Rscript AUC_coexprDist_withParalogs_sortNoDup_permut.R\n"))

options(scipen=100)

buildTable <- TRUE

printAndLog <- function(text, logFile = ""){
  cat(text)
  cat(text, append =T , file = logFile)
}

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
#suppressPackageStartupMessages(library(ggstatsplot, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
#suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
# update 21.08.2019
# to avoid requiring flux
auc <- function (x, y, thresh = NULL, dens = 100, sort.x = TRUE) {
    x <- x[!is.na(x)]
    y <- y[!is.na(x)]
    x <- x[!is.na(y)]
    y <- y[!is.na(y)]
    if (sort.x) {
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
    }
    idx = 2:length(x)
    x <- as.vector(apply(cbind(x[idx - 1], x[idx]), 1, function(x) seq(x[1], 
        x[2], length.out = dens)))
    y <- as.vector(apply(cbind(y[idx - 1], y[idx]), 1, function(x) seq(x[1], 
        x[2], length.out = dens)))
    if (!is.null(thresh)) {
        y.0 <- y <= thresh
        y[y.0] <- thresh
    }
    idx = 2:length(x)
    integral <- as.double((x[idx] - x[idx - 1]) %*% (y[idx] + 
        y[idx - 1]))/2
    integral
}


axisLabSize <- 12
legendSize <- 10
plotTitSize <- 14

mytheme <- theme(
  # top, right, bottom and left
  plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "lines"),
  plot.title = element_text(hjust = 0.5, face = "bold", size=plotTitSize, vjust=1),
  plot.subtitle = element_text(hjust = 0.5, face = "bold", size=plotTitSize-2, vjust=1),
  panel.background = element_rect(fill = "white", colour = NA), 
  panel.border = element_rect(fill = NA, colour = "grey20"), 
  panel.grid.major = element_line(colour = "grey92"), 
  panel.grid.minor = element_line(colour = "grey92", size = 0.25), 
  strip.background = element_rect(fill = "grey85", colour = "grey20"), 
  #legend.key = element_rect(fill = "white", colour = NA), 
  axis.line.x = element_line(size = .3, color = "black"),
  axis.line.y = element_line(size = .3, color = "black"),
  axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=axisLabSize),
  axis.text.x = element_text(color="black", hjust=0.5,vjust = 1, size=axisLabSize),
  axis.title.y = element_text(color="black", size=axisLabSize+1),
  axis.title.x = element_text(color="black", size=axisLabSize+1),
  legend.text = element_text(size=legendSize),
  legend.key.height = unit(1.5,"cm"),
  legend.key = element_blank()
)


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

### HARD CODED
corMethod <- "pearson"
# for plotting:
# look at coexpression ~ distance up to distLimit bp
distLimit <- 500 * 10^3
fitMeth <- "loess"

# nbr of points for loess fit to take the AUC
nbrLoessPoints <- 1000

scatterFontSizeLabel <- 14
scatterFontSizeTitle <- 12

# UPDATE 30.06.2018:
# -> check that always $gene1 < $gene2 before left_join !!!



# Rscript AUC_coexprDist_withParalogs_sortNoDup_permut.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad

### RETRIEVE FROM COMMAND LINE

#### !!! change for otherTADfile: <<< GENE DATA DO NOT CHANGE
# retrieve the sameTAD data frame from:
# file.path("CREATE_SAME_TAD_SORTNODUP", curr_TADlist, "all_TAD_pairs.Rdata")

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)


txt <- paste0("> Parameters retrieved from command line:\n")
stopifnot(length(args) == 2)
curr_TADlist <- args[1]
curr_dataset <- args[2]

outFold <- file.path("AUC_COEXPRDIST_WITHPARALOGS_SORTNODUP_PERMUT",  curr_TADlist, paste0(curr_dataset))
dir.create(outFold, recursive = TRUE)

logFile <- file.path(outFold, paste0("coexpr_dist_withParalogs_logFile.txt"))  
file.remove(logFile)

printAndLog(txt, logFile)
txt <- paste0("... curr_TADlist = ",  curr_TADlist, "\n")
printAndLog(txt, logFile)
txt <- paste0("... curr_dataset = ",  curr_dataset, "\n")
printAndLog(txt, logFile)



txt <- paste0("> ! Hard-coded parameters:\n")
printAndLog(txt, logFile)
txt <- paste0("... corMethod = ",  corMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... buildTable = ",  as.character(buildTable), "\n")
printAndLog(txt, logFile)
txt <- paste0("... distLimit = ",  distLimit, "\n")
printAndLog(txt, logFile)
txt <- paste0("... fitMeth = ",  fitMeth, "\n")
printAndLog(txt, logFile)

mycols <- c("same TAD" ="darkorange1" , "diff. TAD"="darkslateblue",  "paralogs + same TAD"="violetred1", "paralogs + diff. TAD" = "lightskyblue")

sameTADcol <- mycols["same TAD"]
diffTADcol <- mycols["diff. TAD"]
paraSameTADcol <- mycols["paralogs + same TAD"]
paraDiffTADcol <- mycols["paralogs + diff. TAD"]

plotType <- "png"
# myHeight <- ifelse(plotType == "png", 400, 7)
# myWidth <- ifelse(plotType == "png", 600, 10)
myHeight <- ifelse(plotType == "png", 400, 5)
myWidth <- ifelse(plotType == "png", 600, 6)

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

utilsDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg")
source(file.path(utilsDir, "coreg_utils_ggscatterhist.R"))

### UPDATE 08.01.19 => USE TAD LIST SPECIFIC FILES = TISSUE SPECIFIC FILES
distFile <- file.path("CREATE_DIST_SORTNODUP", curr_TADlist, "all_dist_pairs.Rdata")
stopifnot(file.exists(distFile))

paralogFile <- file.path("CREATE_PARALOG_SORTNODUP", "all_paralog_pairs.Rdata")
stopifnot(file.exists(paralogFile))

# CHANGED 08.01.19 !!!
coexprFile <- file.path("CREATE_COEXPR_SORTNODUP", curr_TADlist,  paste0(curr_dataset), corMethod, "coexprDT.Rdata")
stopifnot(file.exists(coexprFile))

## HERE -> FAMILY FROM PERMUT
sameTADfile <- file.path("CREATE_SAME_TAD_PERMUTG2T_SORTNODUP", curr_TADlist, curr_dataset, "all_TAD_pairs.Rdata")
stopifnot(file.exists(sameTADfile))

dataset_pipDir <- file.path("PIPELINE", "OUTPUT_FOLDER", curr_TADlist, curr_dataset) # used to retrieve gene list
stopifnot(dir.exists(dataset_pipDir))

txt <- paste0("... distFile = ",  distFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... coexprFile = ",  coexprFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... sameTADfile = ",  sameTADfile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... sameParalogFile = ",  paralogFile, "\n")
printAndLog(txt, logFile)

################################################ DATA PREPARATION

if(buildTable) {
  
  cat(paste0("... load data\n"))
  all_paralog_pairs <- get(load(paralogFile))
  all_paralog_pairs$gene1 <- as.character(all_paralog_pairs$gene1)
  all_paralog_pairs$gene2 <- as.character(all_paralog_pairs$gene2)
  stopifnot(all_paralog_pairs$gene1 < all_paralog_pairs$gene2)
  
  cat(paste0("... load DIST data\t", distFile, "\t", Sys.time(), "\t"))
  load(distFile)
  cat(paste0(Sys.time(), "\n"))
  head(all_dist_pairs)
  nrow(all_dist_pairs)
  all_dist_pairs$gene1 <- as.character(all_dist_pairs$gene1)
  all_dist_pairs$gene2 <- as.character(all_dist_pairs$gene2)
  # UPDATE 30.06.2018
  stopifnot(all_dist_pairs$gene1 < all_dist_pairs$gene2)
  
  cat(paste0("... load TAD data\t", sameTADfile, "\t", Sys.time(), "\t"))
  ### =>>> CHANGED HERE FOR OTHER TAD FILE !!!
  load(sameTADfile)
  cat(paste0(Sys.time(), "\n"))
  head(all_TAD_pairs)
  nrow(all_TAD_pairs)
  all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
  all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
  # UPDATE 30.06.2018
  stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
  
  cat(paste0("... load COEXPR data\t",coexprFile, "\t", Sys.time(), "\t"))
  load(coexprFile)
  cat(paste0(Sys.time(), "\n"))
  head(coexprDT)
  nrow(coexprDT)
  coexprDT$gene1 <- as.character(coexprDT$gene1)
  coexprDT$gene2 <- as.character(coexprDT$gene2)
  all_TAD_pairs$gene2
  # UPDATE 30.06.2018
  stopifnot(coexprDT$gene1 < coexprDT$gene2)
  
  #============================== RETRIEVE PIPELINE DATA FOR THIS DATASET - USED ONLY FOR GENE LIST
  script0_name <- "0_prepGeneData"
  geneFile <- file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata")
  cat(paste0("... load GENELIST file\t",geneFile, "\t", Sys.time(), "\t"))
  pipeline_geneList <- eval(parse(text = load(geneFile)))
  pipeline_geneList <- as.character(pipeline_geneList)
  
  dataset_dist_pair <- all_dist_pairs[all_dist_pairs$gene1 %in% pipeline_geneList & 
                                        all_dist_pairs$gene2 %in% pipeline_geneList,]
  
  dataset_dist_pairs_limit <- dataset_dist_pair[dataset_dist_pair$dist <= distLimit,]
  head(dataset_dist_pairs_limit)
  nrow(dataset_dist_pairs_limit)
  
  dataset_TAD_pairs <- all_TAD_pairs[all_TAD_pairs$gene1 %in% pipeline_geneList & 
                                       all_TAD_pairs$gene2 %in% pipeline_geneList,]
  head(dataset_TAD_pairs)
  nrow(dataset_TAD_pairs)
  
  # START MERGING DATA 
  cat(paste0("... merge DIST - TAD data\t", Sys.time(), "\t"))
  dataset_dist_TAD_DT <- left_join(dataset_dist_pairs_limit, dataset_TAD_pairs, by=c("gene1", "gene2"))
  cat(paste0(Sys.time(), "\n"))
  
  dataset_dist_TAD_DT$sameTAD <- ifelse(is.na(dataset_dist_TAD_DT$region), 0, 1)
  
}





  if(buildTable){
    
    

    
    cat(paste0(Sys.time(), "\n"))
    head(all_paralog_pairs)
    nrow(all_paralog_pairs)
    all_paralog_pairs$gene1 <- as.character(all_paralog_pairs$gene1)
    all_paralog_pairs$gene2 <- as.character(all_paralog_pairs$gene2)
    stopifnot(all_paralog_pairs$gene1 < all_paralog_pairs$gene2)
    
    dataset_paralog_pairs <- all_paralog_pairs[all_paralog_pairs$gene1 %in% pipeline_geneList & 
                                               all_paralog_pairs$gene2 %in% pipeline_geneList,]
    head(dataset_paralog_pairs)
    
    cat(paste0("... merge PARALOG data\t", Sys.time(), "\t"))
    dataset_dist_TAD_para_DT <- left_join(dataset_dist_TAD_DT, dataset_paralog_pairs, by=c("gene1", "gene2"))
    cat(paste0(Sys.time(), "\n"))
    dataset_dist_TAD_para_DT$paralogs <- ifelse(is.na(dataset_dist_TAD_para_DT$paralogs), 0, 1)
    
    cat(paste0("... merge COEXPR data\t", Sys.time(), "\t"))
    dataset_dist_TAD_paralog_coexpr_DT <- left_join(dataset_dist_TAD_para_DT, coexprDT, by=c("gene1", "gene2"))
    cat(paste0(Sys.time(), "\n"))
    
    allData_dt <- dataset_dist_TAD_paralog_coexpr_DT
    allData_dt$region <- NULL
    allData_dt <- na.omit(allData_dt)
    
    outFile <-file.path(outFold, paste0( "allData_dt.Rdata"))
    save(allData_dt, file = outFile, version=2)
    cat(paste0("... written: ", outFile, "\n"))    
    
  } else{
    outFile <-file.path(outFold, paste0( "allData_dt.Rdata"))
    load(outFile)

  }
  
  nrow(allData_dt)
  allData_dt$dist_kb <- allData_dt$dist/1000
  
  allData_dt$curve1 <-  ifelse(allData_dt$sameTAD == "0", "diff. TAD", "same TAD")
  
  allData_dt$curve2 <-  ifelse(allData_dt$paralogs == "0", NA,
                             ifelse(allData_dt$sameTAD == "0", "paralogs + diff. TAD", "paralogs + same TAD"))
  
  sameTAD_DT <- allData_dt[allData_dt$sameTAD == 1,c("gene1", "gene2", "coexpr", "dist", "dist_kb")]
  sameTAD_DT <- na.omit(sameTAD_DT)
  sameTAD_DT <- sameTAD_DT[order(sameTAD_DT$dist_kb),]
  sameTAD_DT$nPair <- 1:nrow(sameTAD_DT)
  sameTAD_DT$label <- "same TAD"
  
  diffTAD_DT <- allData_dt[allData_dt$sameTAD == 0,c("gene1", "gene2",  "coexpr", "dist", "dist_kb")]
  diffTAD_DT <- na.omit(diffTAD_DT)
  diffTAD_DT <- diffTAD_DT[order(diffTAD_DT$dist_kb),]
  diffTAD_DT$nPair <- 1:nrow(diffTAD_DT)
  diffTAD_DT$label <- "diff. TAD"
  
  para_sameTAD_DT <- allData_dt[allData_dt$paralogs == 1 & allData_dt$sameTAD == 1 ,c("gene1", "gene2", "coexpr", "dist",  "dist_kb")]
  para_sameTAD_DT <- na.omit(para_sameTAD_DT)
  para_sameTAD_DT <- para_sameTAD_DT[order(para_sameTAD_DT$dist_kb),]
  para_sameTAD_DT$nPair <- 1:nrow(para_sameTAD_DT)
  para_sameTAD_DT$label <- "paralogs + same TAD"
  
  para_diffTAD_DT <- allData_dt[allData_dt$paralogs == 1 & allData_dt$sameTAD == 0 ,c("gene1", "gene2",  "coexpr", "dist", "dist_kb")]
  para_diffTAD_DT <- na.omit(para_diffTAD_DT)
  para_diffTAD_DT <- para_diffTAD_DT[order(para_diffTAD_DT$dist_kb),]
  para_diffTAD_DT$nPair <- 1:nrow(para_diffTAD_DT)
  para_diffTAD_DT$label <- "paralogs + diff. TAD"
  
  stopifnot(is.numeric(sameTAD_DT$dist[1]))
  stopifnot(is.numeric(sameTAD_DT$coexpr[1]))
  stopifnot(is.numeric(diffTAD_DT$dist[1]))
  stopifnot(is.numeric(diffTAD_DT$coexpr[1]))
  stopifnot(is.numeric(para_sameTAD_DT$dist[1]))
  stopifnot(is.numeric(para_sameTAD_DT$coexpr[1]))
  stopifnot(is.numeric(para_diffTAD_DT$dist[1]))
  stopifnot(is.numeric(para_diffTAD_DT$coexpr[1]))
  
  #***
  
  if(fitMeth == "loess") {
    my_ylab <- paste0("Gene pair coexpression (", corMethod, ", qqnormDT)")
    my_xlab <- paste0("Distance between the 2 genes (kb)")
    my_sub <- paste0(curr_dataset)
    
    # PREDICT WITH ORIGINAL DISTANCE VALUES
    my_xlab <- paste0("Distance between the 2 genes (bp)")
    
    diffTAD_mod <- loess(coexpr ~ dist, data = diffTAD_DT)
    sameTAD_mod <- loess(coexpr ~ dist, data = sameTAD_DT)
    
    smooth_vals_sameTAD <- predict(sameTAD_mod, sort(sameTAD_DT$dist))
    smooth_vals_diffTAD <- predict(diffTAD_mod, sort(diffTAD_DT$dist))
    
    auc_diffTAD_obsDist <- auc(x = sort(diffTAD_DT$dist), y = smooth_vals_diffTAD)
    auc_sameTAD_obsDist <- auc(x = sort(sameTAD_DT$dist), y = smooth_vals_sameTAD)
    
    outFile <- file.path(outFold, paste0(curr_TADlist, "_", curr_dataset, "_sameTAD_diffTAD_loessFit_originalDist", ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(allData_dt$dist), 
         ylim = range(c(smooth_vals_sameTAD, smooth_vals_diffTAD)),
         xlab=my_xlab, 
         ylab=my_ylab,
         main=paste0(curr_TADlist, " - ", curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = "observed distance values", side = 3)
    lines( x = sort(sameTAD_DT$dist), y = smooth_vals_sameTAD, col = sameTADcol)
    lines( x = sort(diffTAD_DT$dist), y = smooth_vals_diffTAD, col = diffTADcol)
    legend("topright", 
           legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_obsDist, 2), ")"), paste0("diffTAD\n(AUC=", round(auc_diffTAD_obsDist, 2))),  
           col = c(sameTADcol, diffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # PREDICT WITH DISTANCE VECTOR
    distVect <- seq(from=0, to = distLimit, length.out = nbrLoessPoints)
    smooth_vals_sameTAD_distVect <- predict(sameTAD_mod, distVect)
    smooth_vals_diffTAD_distVect <- predict(diffTAD_mod, distVect)
    
    auc_diffTAD_distVect <- auc(x = distVect, y = smooth_vals_diffTAD_distVect)
    auc_sameTAD_distVect <- auc(x = distVect, y = smooth_vals_sameTAD_distVect)
    
    outFile <- file.path(outFold, "auc_diffTAD_distVect.Rdata")
    save(auc_diffTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "auc_sameTAD_distVect.Rdata")
    save(auc_sameTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "smooth_vals_sameTAD_distVect.Rdata")
    save(smooth_vals_sameTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "smooth_vals_diffTAD_distVect.Rdata")
    save(smooth_vals_diffTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "distVect.Rdata")
    save(distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(curr_TADlist, "_", curr_dataset, "_sameTAD_diffTAD_loessFit_vectDist.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(distVect), 
         ylim = range(c(na.omit(smooth_vals_sameTAD_distVect), na.omit(smooth_vals_diffTAD_distVect))),
         xlab=my_xlab,
         ylab=my_ylab,
         main=paste0(curr_TADlist, " - ", curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = paste0("distance values seq from 0 to ", distLimit, " (# points = ", nbrLoessPoints, ")"), side = 3)
    lines( x = distVect, y = smooth_vals_sameTAD_distVect, col = sameTADcol)
    lines( x = distVect, y = smooth_vals_diffTAD_distVect, col = diffTADcol)
    legend("topright", 
           legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_distVect, 2), ")"), paste0("diffTAD\n(AUC=", round(auc_diffTAD_distVect, 2))), 
           col = c(sameTADcol, diffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ################################ DO THE SAME FOR PARALOGS SAME TAD VS. PARALOGS DIFF TAD
    
    
    paraDiffTAD_mod <- loess(coexpr ~ dist, data = para_diffTAD_DT)
    paraSameTAD_mod <- loess(coexpr ~ dist, data = para_sameTAD_DT)
    
    
    smooth_vals_paraSameTAD <- predict(paraSameTAD_mod, sort(para_sameTAD_DT$dist))
    smooth_vals_paraDiffTAD <- predict(paraDiffTAD_mod, sort(para_diffTAD_DT$dist))
    
    auc_paraDiffTAD_obsDist <- auc(x = sort(para_diffTAD_DT$dist), y = smooth_vals_paraDiffTAD)
    auc_paraSameTAD_obsDist <- auc(x = sort(para_sameTAD_DT$dist), y = smooth_vals_paraSameTAD)
    
    outFile <- file.path(outFold, paste0(curr_TADlist, "_", curr_dataset, "_paraSameTAD_paraDiffTAD_loessFit_originalDist", ".", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(allData_dt$dist), 
         ylim = range(c(smooth_vals_paraSameTAD, smooth_vals_paraDiffTAD)),
         xlab=my_xlab, 
         ylab=my_ylab,
         main=paste0(curr_TADlist, " - ", curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = "observed distance values", side = 3)
    lines( x = sort(para_sameTAD_DT$dist), y = smooth_vals_paraSameTAD, col = paraSameTADcol)
    lines( x = sort(para_diffTAD_DT$dist), y = smooth_vals_paraDiffTAD, col = paraDiffTADcol)
    legend("topright", 
           legend=c(paste0("paraSameTAD\n(AUC=", round(auc_paraSameTAD_obsDist, 2), ")"), 
                    paste0("paraDiffTAD\n(AUC=", round(auc_paraDiffTAD_obsDist, 2))),  
           col = c(paraSameTADcol, paraDiffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    # PREDICT WITH DISTANCE VECTOR
    smooth_vals_paraSameTAD_distVect <- predict(paraSameTAD_mod, distVect)
    smooth_vals_paraDiffTAD_distVect <- predict(paraDiffTAD_mod, distVect)
    
    auc_paraDiffTAD_distVect <- auc(x = distVect, y = smooth_vals_paraDiffTAD_distVect)
    auc_paraSameTAD_distVect <- auc(x = distVect, y = smooth_vals_paraSameTAD_distVect)
    
    
    outFile <- file.path(outFold, "diffTAD_mod.Rdata")
    save(diffTAD_mod, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "paraDiffTAD_mod.Rdata")
    save(paraDiffTAD_mod, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "sameTAD_mod.Rdata")
    save(sameTAD_mod, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "paraSameTAD_mod.Rdata")
    save(paraSameTAD_mod, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    paraDiffTAD_obsDist <- para_diffTAD_DT$dist
    outFile <- file.path(outFold, "paraDiffTAD_obsDist.Rdata")
    save(paraDiffTAD_obsDist, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    paraSameTAD_obsDist <- para_sameTAD_DT$dist
    outFile <- file.path(outFold, "paraSameTAD_obsDist.Rdata")
    save(paraSameTAD_obsDist, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    outFile <- file.path(outFold, "auc_paraDiffTAD_distVect.Rdata")
    save(auc_paraDiffTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "auc_paraSameTAD_distVect.Rdata")
    save(auc_paraSameTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "smooth_vals_paraSameTAD_distVect.Rdata")
    save(smooth_vals_paraSameTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, "smooth_vals_paraDiffTAD_distVect.Rdata")
    save(smooth_vals_paraDiffTAD_distVect, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFold, paste0(curr_TADlist, "_", curr_dataset, "_paraSameTAD_paraDiffTAD_loessFit_vectDist.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(distVect), 
         ylim = range(c(na.omit(smooth_vals_paraSameTAD_distVect), na.omit(smooth_vals_paraDiffTAD_distVect))),
         # xlab="", 
         # ylab="",
         xlab=my_xlab,
         ylab=my_ylab,
         main=paste0(curr_TADlist, " - ", curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = paste0("distance values seq from 0 to ", distLimit, " (# points = ", nbrLoessPoints, ")"), side = 3)
    lines( x = distVect, y = smooth_vals_paraSameTAD_distVect, col = paraSameTADcol)
    lines( x = distVect, y = smooth_vals_paraDiffTAD_distVect, col = paraDiffTADcol)
    legend("topright", 
           legend=c(paste0("paraSameTAD\n(AUC=", round(auc_paraSameTAD_distVect, 2), ")"), 
                    paste0("paraDiffTAD\n(AUC=", round(auc_paraDiffTAD_distVect, 2))), 
           col = c(paraSameTADcol, paraDiffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ################################################################ 
    outFile <- file.path(outFold, paste0(curr_TADlist, "_", curr_dataset, "_sameTAD_diffTAD_paraSameTAD_paraDiffTAD_loessFit_vectDist.", plotType))
    do.call(plotType, list(outFile, height = myHeight, width = myWidth))
    plot(NULL,
         xlim = range(distVect), 
         ylim = range(c(na.omit(smooth_vals_sameTAD_distVect), 
                        na.omit(smooth_vals_paraSameTAD_distVect),
                        na.omit(smooth_vals_diffTAD_distVect),
                        na.omit(smooth_vals_paraDiffTAD_distVect))),
         # xlab="", 
         # ylab="",
         xlab=my_xlab,
         ylab=my_ylab,
         main=paste0(curr_TADlist, " - ", curr_dataset, ": coexpr ~ dist loess fit"))
    mtext(text = paste0("distance values seq from 0 to ", distLimit, " (# points = ", nbrLoessPoints, ")"), side = 3)
    lines( x = distVect, y = smooth_vals_paraSameTAD_distVect, col = paraSameTADcol)
    lines( x = distVect, y = smooth_vals_paraDiffTAD_distVect, col = paraDiffTADcol)
    lines( x = distVect, y = smooth_vals_sameTAD_distVect, col = sameTADcol)
    lines( x = distVect, y = smooth_vals_diffTAD_distVect, col = diffTADcol)
    legend("topright", 
           legend=c(paste0("sameTAD\n(AUC=", round(auc_sameTAD_distVect, 2), ")"), 
                    paste0("diffTAD\n(AUC=", round(auc_diffTAD_distVect, 2)), 
                     paste0("paraSameTAD\n(AUC=", round(auc_paraSameTAD_distVect, 2), ")"), 
                     paste0("paraDiffTAD\n(AUC=", round(auc_paraDiffTAD_distVect, 2))), 
           col = c(sameTADcol, diffTADcol,
                   paraSameTADcol, paraDiffTADcol),
           lty=1,
           bty = "n")
    
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    ################################################################ 
    auc_values <- list(
      auc_diffTAD_distVect = auc_diffTAD_distVect,
      auc_sameTAD_distVect = auc_sameTAD_distVect,
      auc_ratio_same_over_diff_distVect = auc_sameTAD_distVect/auc_diffTAD_distVect,
      auc_diffTAD_obsDist = auc_diffTAD_obsDist,
      auc_sameTAD_obsDist = auc_sameTAD_obsDist,
      auc_ratio_same_over_diff_obsDist = auc_sameTAD_distVect/auc_diffTAD_obsDist,
      
      auc_paraDiffTAD_distVect = auc_paraDiffTAD_distVect,
      auc_paraSameTAD_distVect = auc_paraSameTAD_distVect,
      auc_ratio_para_same_over_diff_distVect = auc_paraSameTAD_distVect/auc_paraDiffTAD_distVect,
      auc_paraDiffTAD_obsDist = auc_paraDiffTAD_obsDist,
      auc_paraSameTAD_obsDist = auc_paraSameTAD_obsDist,
      auc_ratio_para_same_over_diff_obsDist = auc_paraSameTAD_distVect/auc_paraDiffTAD_obsDist
    )
    
    outFile <- file.path(outFold, paste0("auc_values.Rdata"))
    save(auc_values, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
  } else{
    stop("only loess implemented yet\n")
  }

    
######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


