

options(scipen=100)

# Rscript presentation_figures3.R

script_name <- "presentation_figures3.R"

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

outFolder <- "PRESENTATION_FIGURES3"
dir.create(outFolder, recursive = TRUE)


ex_hicds <- "GSE99051_786_O_40kb"
ex_exprds <- "TCGAkich_norm_kich"
hicds_tit <- "786-O"
exprds_tit <- "norm_vs_kich"


plotCex <- 1.2

final_DT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

ex_DT <- final_DT[final_DT$hicds == ex_hicds &
                    final_DT$exprds == ex_exprds,
                  ]


signifThresh <- 0.01
ex_DT_signif <- ex_DT[ex_DT$adjPvalComb <= signifThresh,]


# par(bty="l")
# densplot(
#   main = paste0(hicds_tit, " - ", exprds_tit),
#   x = ex_DT$meanLogFC,
#   y = ex_DT$meanCorr,
#   xlab = "meanLogFC",
#   ylab = "meanIntraCorr",
#   cex.axis=plotCex,
#   cex.lab=plotCex
# )
# mtext(side=3, text = paste0("all TADs - n = ", nrow(ex_DT)))


# 
# par(bty="l")
# densplot(
#   main = paste0(hicds_tit, " - ", exprds_tit),
#   x = ex_DT_signif$meanLogFC,
#   y = ex_DT_signif$meanCorr,
#   xlab = "meanLogFC",
#   ylab = "meanIntraCorr",
#   cex.axis=plotCex,
#   cex.lab=plotCex
# )
# mtext(side=3, text = paste0("signif. TADs - n = ", nrow(ex_DT_signif)))
# legend("topleft", legend=paste0("adj. p-val <= ", signifThresh))




outFile <- file.path(outFolder, paste0(ex_hicds, "_", ex_exprds, "_meanIntraCorr_meanLogFC_dotplot_with_signif.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

ex_DT$dotcol <- ifelse(ex_DT$adjPvalComb <= signifThresh,"red", "darkgrey")
par(bty="l")
plot(
  main = paste0(hicds_tit, " - ", exprds_tit),
  x = ex_DT$meanLogFC,
  y = ex_DT$meanCorr,
  xlab = "meanLogFC",
  ylab = "meanIntraCorr",
  pch=16,
  col = ex_DT$dotcol,
  cex=0.7,
  cex.axis=plotCex,
  cex.lab=plotCex
)
mtext(side=3, text = paste0("# TADs = ", nrow(ex_DT), "; # signif. TADs = ", nrow(ex_DT_signif)))
# legend("topleft", legend=paste0("adj. p-val <= ", signifThresh), pch=16, col="red", bty="n")
legend("bottomright", legend=paste0("adj. p-val <= ", signifThresh), pch=16, col="red", bty="n")



foo <- dev.off()

stop("--ok\n")


final_DT_signif <- final_DT[final_DT$adjPvalComb <=signifThresh,]

nrow(final_DT)# 100884
nrow(final_DT_signif) # 1453
nrow(final_DT_signif)/58
# [1] 25.05172
nrow(final_DT)/58
# [1] 1739.379



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




