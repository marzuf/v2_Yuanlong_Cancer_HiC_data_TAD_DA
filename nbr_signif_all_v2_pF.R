


require(foreach)
require(doMC)
registerDoMC(40)


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


# Rscript nbr_signif_all_v2_pF.R

outFolder <- "NBR_SIGNIF_ALL_V2_PF"
dir.create(outFolder)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <-  "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4

jitterCol <- "blue"

qt_dt<- read.delim("hicds_sparsity.csv", sep="\t", header = TRUE)

pipOutFolder <- "PIPELINE/OUTPUT_FOLDER"

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[! (grepl("RANDOMSUB", all_hicds) | grepl("RANDOMNBR", all_hicds) | grepl( "RANDOMSHIFT", all_hicds))]
hicds = all_hicds[1]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

tad_signifThresh <- 0.01

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    
    combPval_file <- file.path(pipOutFolder, hicds, exprds, "11sameNbrPF_runEmpPvalCombined", "emp_pval_combined.Rdata")
    
    if(!file.exists(combPval_file))return(NULL)
    
    stopifnot(file.exists(combPval_file))
    pval <- get(load(combPval_file))
    adjPval <- p.adjust(pval, method="BH")
    
    data.frame(
      hicds=hicds,
      exprds=exprds,
      nSignif=sum(adjPval <= tad_signifThresh),
      nTot = length(adjPval),
      stringsAsFactors = FALSE
    )
  }
}

save(all_dt,file=file.path(outFolder, "all_dt.Rdata"), version=2)

all_dt$hicds_lab <- gsub(".+_(.+?)_40kb","\\1",  all_dt$hicds)
all_dt$hicds_lab <- ifelse(grepl("RANDOM",all_dt$hicds_lab) | grepl("PERMUT", all_dt$hicds_lab),all_dt$hicds_lab,  
                           "OBSERVED" )

# all_dt$hicds_lab <- factor(all_dt$hicds_lab, 
#                            levels=c("OBSERVED", "RANDOMMIDPOS", "RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT", "PERMUTG2T"))

subTit <- paste0(paste0("# ", names(table(all_dt$hicds_lab)), "=", as.numeric(table(all_dt$hicds_lab))), collapse="; ")
legText <- paste0("# ", names(table(all_dt$hicds_lab)), "=", as.numeric(table(all_dt$hicds_lab)))


# stopifnot(diff(table(all_dt$hicds_lab)) == 0)  # TO UNCOMMENT LATER!

outFile <- file.path(outFolder, paste0("nbrSignifTADs_densityplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  split(all_dt$nSignif, all_dt$hicds_lab),
  plotTit = paste0("# signif. TADs")
)
legend("topright", legend=legText, bty="n", cex=0.9)

mtext(side=3, text = paste0("adj. comb. p-val <= ", tad_signifThresh), font = 3)
foo <- dev.off()

outFile <- file.path(outFolder, paste0("nbrSignifTADs_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
par(mar = par()$mar + c(10,3,0,0))
boxplot(nSignif~hicds_lab, data = all_dt, outline=FALSE,
        main = "# signif. TADs", xlab="", ylab="", 
        cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)

stripchart(nSignif~hicds_lab, vertical = TRUE, data = all_dt,
           method = "jitter", add = TRUE, pch = 20, col = jitterCol)

legend("topright", legend=legText, bty="n", cex=0.9)

mtext(side=2, text="# signif. TADs", cex=plotCex, line=5)

mtext(side=3, text = paste0("adj. comb. p-val <= ", tad_signifThresh), font = 3)
foo <- dev.off()

all_dt$ratioSignif <- all_dt$nSignif/all_dt$nTot
outFile <- file.path(outFolder, paste0("ratioSignifTADs_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
par(mar = par()$mar + c(10,3,0,0))
boxplot(ratioSignif~hicds_lab, outline=FALSE,
        data = all_dt, main = "Ratio signif. TADs", 
        xlab="", ylab="", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)

stripchart(ratioSignif~hicds_lab, vertical = TRUE, data = all_dt,
           method = "jitter", add = TRUE, pch = 20, col = jitterCol)
legend("topright", legend=legText, bty="n", cex=0.9)

mtext(side=2, text="Ratio signif. TADs", cex=plotCex, line=5)
mtext(side=3, text = paste0("adj. comb. p-val <= ", tad_signifThresh), font = 3)
foo <- dev.off()



# boxplot(NUMS ~ GRP, data = ddf, lwd = 2, ylab = 'NUMS')
# stripchart(NUMS ~ GRP, vertical = TRUE, data = ddf, 
#            method = "jitter", add = TRUE, pch = 20, col = 'blue')



