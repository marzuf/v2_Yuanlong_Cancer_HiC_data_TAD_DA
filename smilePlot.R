#!/usr/bin/Rscript

# Rscript smilePlot.R

script_name <- "smilePlot.R"

options(scipen=100)

startTime <- Sys.time()

suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
require(foreach)
require(doMC)
registerDoMC(40)


# set files and folders
script8_name <- "8cOnlyRatioDownFastSave_runAllDown"
outFolder <- file.path("SMILEPLOT")
dir.create(outFolder, recursive=TRUE)
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

# settings for plotting
thresh1 <- 0.75
thresh2 <- 1-thresh1
nRandomInit <- 100000
keepPermut <- 1000 # because converted as vector -> too big if keep 100000
nRandom <- keepPermut
curr_ratio_type <- "ratioDown"

## TO SET
## hard-coded:
obsCol <- "bisque"
obsColText <- "bisque2"
permutCol <- "mediumseagreen"
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 1028 , 15)
myWidth <- ifelse(plotType == "png", 686, 10)

histBreakStep <- 0.1
stopifnot( (1/histBreakStep) %% 1 == 0 )
#***********************************************************************************

stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
mainFolder <- "."
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"

setDir <-""

all_missing_data <- foreach(hicds = all_hicds) %dopar% {
  missing_data <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {


    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_smilePlot.", plotType))
if(file.exists(outFile)) return(NULL)

    pipOutFold <- file.path(pipFolder, hicds, exprds)
    
    setDir <- ""
    
    cat(paste0("setDir =", setDir, "\n"))
    
    settingF <- file.path(settingFolder,hicds, paste0("run_settings_", exprds, ".R"))
    source(settingF)
    
    curr_ratio_breaks <- seq(0,1,histBreakStep)
    
    cat(paste0("*** START smile plot for: ", hicds, " - " ,exprds, " - ", curr_ratio_type, "\n"))
    
    cat("... load observed data\n")
    
    obs_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
    
    cat("... load permut data\n")
    permut_currDown_DT <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"))))
    
    cat(ncol(permut_currDown_DT), "\n")
    
    if(ncol(permut_currDown_DT) != nRandomInit) return(paste0(hicds, "-", exprds))
    
    stopifnot(ncol(permut_currDown_DT) == nRandomInit)
    
    permut_currDown_DT <- permut_currDown_DT[,1:keepPermut]
    stopifnot(!is.na(permut_currDown_DT))
    stopifnot(grepl("_TAD",  rownames(permut_currDown_DT)))
    
    obsThresh <- mean(obs_curr_down >= thresh1 | obs_curr_down <= thresh2)
    stopifnot( obsThresh >= 0 & obsThresh <= 1)
    allThreshPermut <- unlist(apply(permut_currDown_DT, 2, function(x) mean(x >= thresh1 | x <= thresh2)))
    stopifnot( allThreshPermut >= 0 & allThreshPermut <= 1)
    stopifnot(length(allThreshPermut) == keepPermut)
    
    cat("> Prepare simulation data \n")
    nTAD <- sum(regexpr("_TAD",  rownames(permut_currDown_DT)) > 0)
    stopifnot(nTAD == nrow(permut_currDown_DT))
    cat(paste0("... TADs: ", nTAD, "/", nrow(permut_currDown_DT), "\n"))
    permut_currDown_vect <- as.vector(permut_currDown_DT[grep("_TAD", rownames(permut_currDown_DT)),])
    stopifnot(length(permut_currDown_vect) == ncol(permut_currDown_DT) * nTAD)
    
    cat(paste0("length(permut_currDown_vect) = ", length(permut_currDown_vect), "\n"))
    
    cat("> Prepare observed data \n")
    nTAD <- sum(regexpr("_TAD",  names(obs_curr_down)) > 0)
    cat(paste0("... TADs: ", nTAD, "/", length(obs_curr_down), "\n"))
    obs_curr_down <- obs_curr_down[grep("_TAD", names(obs_curr_down))]
    
    cat("> Prepare data for plotting \n")
    
    obs_rel_hist <- try(hist(obs_curr_down, breaks=curr_ratio_breaks, plot=F))
    #plot(obs_rel_hist)
    # some are not bounded 0-1 => go to next directly !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"WILL NEED REFINEMENT!!!
    if(class(obs_rel_hist) == "try-error") {
      warning(paste0("!!! CANNOT PLOT HISTOGRAM FOR: ", curr_ratio_type, " !!! \n" ))
      next
    }
    cat(paste0("length(permut_currDown_vect) = ", length(permut_currDown_vect), "\n"))
    cat(paste0("length(permut_currDown_vect) = ", length(na.omit(permut_currDown_vect)), "\n"))
    shuff_rel_hist <-  hist(permut_currDown_vect,breaks=curr_ratio_breaks, plot=F)
    # plot(shuff_rel_hist)
    
    cat(paste0("ncol(permut_currDown_DT) = ", ncol(permut_currDown_DT), "\n"))
    cat(paste0("sum(shuff_rel_hist$counts) = ", sum(shuff_rel_hist$counts), "\n"))
    cat(paste0("sum(shuff_rel_hist$counts/nRandom) = ", sum(shuff_rel_hist$counts/nRandom), "\n"))
    cat(paste0("sum(obs_rel_hist$counts) = ", sum(obs_rel_hist$counts),"\n"))
    stopifnot(sum(shuff_rel_hist$counts/nRandom) == sum(obs_rel_hist$counts))   # 830
    stopifnot(sum(shuff_rel_hist$counts) == length(permut_currDown_vect))
    stopifnot(sum(obs_rel_hist$counts) == length(na.omit(obs_curr_down)))
    
    rel_freqValues <- rbind(obs_rel_hist$counts, shuff_rel_hist$counts/nRandom)
    rownames(rel_freqValues) <- c("observed", "randomized")
    
    # FOR THE ERROR BARS:
    # for each randomization, calculate the observed, then divide expected by observed
    # table 1 column = 1 randomization, rows = breaks of the hist.
    rel_HistValues <- apply(permut_currDown_DT[grep("_TAD", rownames(permut_currDown_DT)),],
                            2, function(x) hist(x,breaks=curr_ratio_breaks, plot=F)$counts)
    stopifnot(nrow(rel_HistValues) == length(obs_rel_hist$counts))
    # for each randomization, divide observ./exp.
    rel_obsOverExp <- apply(rel_HistValues, 2, function(x) log2(obs_rel_hist$counts/x))
    # now get the SEM for each break of the hist (i.e. each row)
    rel_obsOverExp_sem <- apply(rel_obsOverExp, 1, function(x) {
      x <- na.omit(x[is.finite(x)])
      sd(x, na.rm=T)/sqrt(length(x))
    })
    
    # Calculate observed/expected ratio
    rel_logRatio <- log2(rel_freqValues["observed",]/rel_freqValues["randomized",])
    # put y coord at the middle of the breaks
    rel_ycoord <- (curr_ratio_breaks * 100)[-1]
    rel_ycoord <- rel_ycoord-(histBreakStep*100)/2
    stopifnot(length(rel_ycoord) == length(rel_logRatio))
    toKeep <- !is.na(rel_logRatio) & abs(rel_logRatio) != "Inf"
    rel_logRatio <- rel_logRatio[toKeep]
    rel_ycoord <- rel_ycoord[toKeep]
    rel_obsOverExp_sem <- rel_obsOverExp_sem[toKeep] 
    my_ylab <- paste0(rel_ycoord - (histBreakStep*100)/2, "-", rel_ycoord + (histBreakStep*100)/2, "%")
    plotTit <- paste0(hicds, " - ", cond1, " vs. ", cond2)
    subtitDir1 <- paste0("kept # permut=", keepPermut, "; >=",  thresh1, "|<=", thresh2, ": ", round(obsThresh*100, 2),
                         "% (avg. permut: ", round(mean(allThreshPermut) * 100, 2), "%)")
    subtitDir2 <- ""
    rel_ycoord <- rel_ycoord/100
    
    ####### PLOT HERE
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_smilePlot.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    par(mfrow = c(2,1))
    # top panel: line plot
    myxlab <-  paste0("% of down-regulated genes per TAD")
    
    # A) TOP PLOT - smile plot
    plot(rel_logRatio ~ rel_ycoord, type="l",axes=F, cex.main=0.9,
         main = plotTit,
         xlab = myxlab, ylab = "log2(Observed/Randomized)")
    mtext(side=3, paste0(subtitDir1))
    box(bty="l")
    axis(2)
    axis(1, at=rel_ycoord, labels = my_ylab)
    abline(h=0, lty=2)
    errbar(x=rel_ycoord, y=rel_logRatio, yplus=rel_logRatio+rel_obsOverExp_sem, yminus=rel_logRatio-rel_obsOverExp_sem, add=T)
    
    # B) BOTTOM PLOT - histogram bar plot
    barplot(rel_freqValues, beside=T, col=c(obsCol, permutCol), ylab="Frequency", 
            xlab=myxlab)
    legend("topright", c( "observed Freq.", "randomized Freq."),
           text.col=c(obsColText, permutCol), bty='n', cex=1)
    mtext(side=3, text = subtitDir2)
    
    rel_ycoord <- (curr_ratio_breaks * 100)[-1]
    rel_ycoord <- rel_ycoord-(histBreakStep*100)/2
    my_ylab <- rel_ycoord
    cat("my_ylab\n", my_ylab, "\n")
    xlabpos <-  seq(2, 1/histBreakStep * 6, by=3)
    cat("xlabpos\n", xlabpos, "\n")
    box(bty="l")
    axis(1)
    axis(2, at=c(0, max(rel_freqValues, na.rm=T)), labels=F)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    return(NULL)
  }
  missing_data
}
outFile <- file.path(outFolder, "all_missing_data.Rdata")
save(all_missing_data, file=outFile, version=2)

###################################################################
txt <- paste0(startTime, "\n", Sys.time(), "\n")
cat(txt)
cat(paste0("*** DONE: ", script_name, "\n"))




