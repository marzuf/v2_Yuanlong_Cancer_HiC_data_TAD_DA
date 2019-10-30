

options(scipen=100)

# Rscript reportfig2_figure2.R

script_name <- "reportfig2_figure2.R"

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

plotType <- "svg"
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

outFolder <- "REPORTFIG2_FIGURE2"
dir.create(outFolder, recursive = TRUE)


signifThresh <- 0.01


final_DT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

totNbr <- aggregate(adjPvalComb ~hicds + exprds, data=final_DT, FUN=length)
colnames(totNbr)[colnames(totNbr) == "adjPvalComb"] <- "totNbrTADs"

signifNbr <- aggregate(adjPvalComb ~hicds + exprds, data=final_DT, FUN=function(x) sum(x <= signifThresh))
colnames(signifNbr)[colnames(signifNbr) == "adjPvalComb"] <- "signifNbrTADs"

plot_dt <- merge(totNbr, signifNbr, by=c("hicds", "exprds"))


plot_dt_m <- melt(plot_dt, id=c("hicds", "exprds"))

plot_dt_m$value_log10 <- log10(plot_dt_m$value)

outFile <- file.path(outFolder, paste0("nbrTADs_tot_signif_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))


boxplot(value_log10~variable, data=plot_dt_m,
        main=paste0("# of TADs"),
        names=c("total", "signif."),
        sub=paste0("adj. p-val <= ", signifThresh),
        xlab = "",
        ylab="# of TADs [log10]",
        cex.axis =axisCex,
        cex.lab =axisCex
        )
mtext(side=3, text=paste0("all datasets - n=", length(unique(paste0(plot_dt$hicds, plot_dt$exprds)))))

legend("topright", legend= c(paste0("mean tot. = ", round(mean(plot_dt$totNbrTADs),2), "\n", 
                                    "mean signif. = ", round(mean(plot_dt$signifNbrTADs),2))), bty="n")


foo <- dev.off()


cat(paste0("... written: ", outFile, "\n"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




