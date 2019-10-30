options(scipen=100)

# Rscript reportfig2_figure1.R

script_name <- "reportfig2_figure1.R"

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


barcol <- "darkorange3"


plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <-  ifelse(plotType=="png", 600, 10)
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 8

fcc_col <- "dodgerblue3"
coexpr_col <- "goldenrod"

script17_name <- "170revision2EZH2_score_auc_pval_permGenes"

famType1 <- "hgnc"
famType2 <- "hgnc_family_short"
auc_coexprdist_fold <- file.path("AUC_COEXPRDIST_WITHFAM_SORTNODUP") 
stopifnot(dir.exists(auc_coexprdist_fold))

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

outFolder <- "REPORTFIG2_FIGURE1"
dir.create(outFolder, recursive = TRUE)

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

minTADsize <- 3

hicds = all_hicds[1]

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    fcc_auc_file <- file.path(pipFolder, hicds, exprds, script17_name, "auc_ratios.Rdata")
    stopifnot(file.exists(fcc_auc_file))
    fcc_auc <- as.numeric(get(load(fcc_auc_file))["prodSignedRatio_auc_permGenes"])
    stopifnot(!is.na(fcc_auc))
    coexpr_auc_file <- file.path(auc_coexprdist_fold, hicds, paste0(exprds, "_", famType1), famType2, "auc_values.Rdata")
    stopifnot(file.exists(coexpr_auc_file))
    coexpr_auc <- get(load(coexpr_auc_file))[["auc_ratio_same_over_diff_distVect"]]
    stopifnot(!is.na(coexpr_auc))
    data.frame(hicds=hicds, exprds=exprds,fcc_auc=fcc_auc, coexpr_auc=coexpr_auc, stringsAsFactors = FALSE)
  } # end-foreach iterating over exprds
  exprds_dt
} # end-foreach iterating over hicds

save(all_dt, file=file.path(outFolder, "all_dt.Rdata"),version=2)

# stop("--ok\n")


all_dt <- all_dt[order(all_dt$fcc_auc, decreasing = TRUE),]
labsymbol <- "\u25CF"
exdataset <- paste0("GSE99051_786_O_40kb", "\n", "TCGAkich_norm_kich")


all_dt$plotlab <- paste0(all_dt$hicds, "\n", all_dt$exprds)
all_dt$plotlab_short <- ifelse(all_dt$plotlab == exdataset, exdataset, labsymbol)
all_dt$plotlab_short <- ifelse(all_dt$plotlab == exdataset, labsymbol, labsymbol)
plotlab_color <- all_cols[all_cmps[all_dt$exprds]]


outFile <- file.path(outFolder, paste0("fcc_barplot_fullnames.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
barp <- barplot(all_dt$fcc_auc, col=barcol, axes=F, las=2, cex.names=0.4, ylab="FCC AUC ratio", cex.lab=1.2, xlab="Datasets")
axis(2, cex.axis=axisCex)
axis(1, at=barp, labels=all_dt$plotlab, las=2, cex.lab=0.2, cex.axis=0.4)
legend("topright", pch=16, col=all_cols, legend=names(all_cols),bty="n")
foo <- dev.off()

outFile <- file.path(outFolder, paste0("fcc_barplot_colSymb.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
barp <- barplot(all_dt$fcc_auc-1,
                ylab="FCC AUC ratio", cex.lab=axisCex, xlab="Datasets",
        col=barcol, axes=F,
        main=paste0("FCC AUC ratio"))
axis(2, at = seq(0, 0.8, by=0.1), labels = seq(0, 0.8, by=0.1)+1, cex.axis=axisCex)
mtext(side=1, at = barp, text=all_dt$plotlab_short, col = plotlab_color, las=2)
mtext(paste0("all datasets - n=", nrow(all_dt)))

ds_count <- all_dt[,c("hicds", "exprds")]
ds_count <- unique(ds_count)
ds_count$dataset <- file.path(ds_count$hicds, ds_count$exprds)
ds_count$cmpType <- all_cmps[ds_count$exprds]
legend("topright", pch=16, col=all_cols, legend=paste0(names(all_cols)," (n=", table(ds_count$cmpType)[names(all_cols)], ")"),
       bty="n", cex = axisCex)

foo <- dev.off()

cat(paste0("... written: ", outFile, "\n"))

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




