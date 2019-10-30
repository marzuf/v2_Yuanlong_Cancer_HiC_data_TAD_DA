

options(scipen=100)

# Rscript reportfig2_figure3.R

script_name <- "reportfig2_figure3.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
require(reshape2)

registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

all_cols[all_cols == "red"] <- "brown3"
all_cols[all_cols == "blue"] <- "darkblue"
all_cols[all_cols == "green"] <- "forestgreen"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

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

outFolder <- "REPORTFIG2_FIGURE3"
dir.create(outFolder, recursive = TRUE)


final_DT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))


signifThresh <- 0.01

final_DT$adjPvalComb_log10 <- -log10(final_DT$adjPvalComb)

outFile <- file.path(outFolder, paste0("meanFC_meanCorr_densityplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
par(mfrow=c(1,2))
densplot(y = final_DT$adjPvalComb_log10,
         x = final_DT$meanLogFC,
         xlab = "meanLogFC",
         ylab = "adj. p-val [-log10]",
         main=paste0("TAD signif. and meanLogFC"),
         # sub=paste0("adj. p-val <= ", signifThresh),
         cex.axis =axisCex,
         cex.lab =axisCex
         )
legend("topleft" , legend=paste0("n=", nrow(final_DT)), bty="n")
abline(h=-log10(signifThresh), lty=2, col="darkgrey")
mtext(side=3, text=paste0("all datasets - n=", length(unique(paste0(final_DT$hicds, final_DT$exprds)))))

densplot(y = final_DT$adjPvalComb_log10,
         x = final_DT$meanCorr,
         xlab = "meanIntraCorr",
         ylab = "adj. p-val [-log10]",
         main=paste0("TAD signif. and meanIntraCorr"),
         # sub=paste0("adj. p-val <= ", signifThresh),
         cex.axis =axisCex,
         cex.lab =axisCex
)
legend("topleft" , legend=paste0("n=", nrow(final_DT)), bty="n")
abline(h=-log10(signifThresh), lty=2, col="darkgrey")
mtext(side=3, text=paste0("all datasets - n=", length(unique(paste0(final_DT$hicds, final_DT$exprds)))))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


sigif_dt <- final_DT[final_DT$adjPvalComb <= signifThresh,]

outFile <- file.path(outFolder, paste0("meanFC_meanCorr_densityplot_signifOnly.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
par(mfrow=c(1,2))
densplot(y = sigif_dt$adjPvalComb_log10,
         x = sigif_dt$meanLogFC,
         xlab = "meanLogFC",
         ylab = "adj. p-val [-log10]",
         main=paste0("TAD signif. and meanLogFC"),
         # sub=paste0("adj. p-val <= ", signifThresh),
         cex.axis =axisCex,
         cex.lab =axisCex
)
legend("topleft" , legend=paste0("n=", nrow(final_DT)), bty="n")
mtext(side=3, text=paste0("all datasets - n=", length(unique(paste0(final_DT$hicds, final_DT$exprds)))))

densplot(y = sigif_dt$adjPvalComb_log10,
         x = sigif_dt$meanCorr,
         xlab = "meanIntraCorr",
         ylab = "adj. p-val [-log10]",
         main=paste0("TAD signif. and meanIntraCorr"),
         # sub=paste0("adj. p-val <= ", signifThresh),
         cex.axis =axisCex,
         cex.lab =axisCex
)
legend("topleft" , legend=paste0("n=", nrow(final_DT)), bty="n")
mtext(side=3, text=paste0("all datasets - n=", length(unique(paste0(final_DT$hicds, final_DT$exprds)))))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




