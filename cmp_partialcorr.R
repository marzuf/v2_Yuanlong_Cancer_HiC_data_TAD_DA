########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript allTADs_and_purity.R\n"))

script_name <- "cmp_partialcorr.R"

# Rscript cmp_partialcorr.R

require(foreach)
require(doMC)
registerDoMC(40)

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

corMet <- "pearson"

plotType <- "png"
myHeight <- 400
myWidth <- myHeight

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
  purity_plot_name <- "EPIC"
} else{
  purity_ds <- ""
  purity_plot_name <- "aran"
}

transfExpr <- "log10"

all_ds_corrPurity_dt <- get(load(file.path("ALLTADS_AND_PURITY", purity_ds, transfExpr, "all_ds_corrPurity_dt.Rdata")))


do_plot <- function(my_x, my_y, ...) {
  plot(x=my_x,
       y=my_y,
       pch=16,
       cex=0.7,
       cex.axis=1.2,
       cex.main=1.2,
       cex.lab=1.2,
       ...)
  addCorr(my_x, my_y, bty="n")
}
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_cols[all_cols=="green"] <- "darkgreen"

outFolder <- file.path("CMP_PARTIALCORR", purity_ds, transfExpr)
dir.create(outFolder, recursive=T)

script4_name <- "4_runMeanTADCorr"
script4partial_name <- "4partial_runMeanTADCorr"

script11_name <- "11sameNbr_runEmpPvalCombined"
script11partial_name <- "11sameNbrPartial_runEmpPvalCombined"

# all_ds=all_ds[1]

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"


# all_ds = file.path(ex_hicds, ex_exprds)
  
all_full_partial_corr_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  cat(paste0("... start: ", ds, "\n"))
  hicds <- dirname(ds)
  exprds <- basename(ds)
  empPval <- get(load(file.path(pipFolder, hicds, exprds, script11_name, "emp_pval_combined.Rdata")))
  partialEmpPval <- get(load(file.path(pipFolder, hicds, exprds, script11partial_name, "emp_pval_combined.Rdata")))
  stopifnot(names(empPval) == names(partialEmpPval))
  
  meanCorr <- get(load(file.path(pipFolder, hicds, exprds,script4_name, "all_meanCorr_TAD.Rdata")))
  partialmeanCorr <- get(load(file.path(pipFolder, hicds, exprds, script4partial_name, "all_meanCorr_TAD.Rdata")))
  stopifnot(names(meanCorr) == names(partialmeanCorr))
  
  
  data.frame(
    dataset =ds,  
    empPval_full_partial_corr = cor(partialEmpPval, empPval),
    meanCorr_full_partial_corr = cor(partialmeanCorr, meanCorr),
    stringsAsFactors = FALSE
  )
   
}
outFile <- file.path(outFolder, "all_full_partial_corr_dt.Rdata")
save(all_full_partial_corr_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

all_full_partial_delta_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  cat(paste0("... start: ", ds, "\n"))
  hicds <- dirname(ds)
  exprds <- basename(ds)
  empPval <- get(load(file.path(pipFolder, hicds, exprds, script11_name, "emp_pval_combined.Rdata")))
  partialEmpPval <- get(load(file.path(pipFolder, hicds, exprds, script11partial_name, "emp_pval_combined.Rdata")))
  stopifnot(names(empPval) == names(partialEmpPval))
  
  meanCorr <- get(load(file.path(pipFolder, hicds, exprds,script4_name, "all_meanCorr_TAD.Rdata")))
  partialmeanCorr <- get(load(file.path(pipFolder, hicds, exprds, script4partial_name, "all_meanCorr_TAD.Rdata")))
  stopifnot(names(meanCorr) == names(partialmeanCorr))
  
  stopifnot(setequal(names(meanCorr) , names(empPval)))
  
  all_regs <- names(meanCorr)
  
  out_dt <- data.frame(
    dataset =ds,  
    region = all_regs,
    empPval_full_partial_delta = as.numeric(empPval[all_regs])-as.numeric(partialEmpPval[all_regs]),
    absMeanCorr_full_partial_delta = abs(as.numeric(meanCorr[all_regs]))-abs(as.numeric(partialmeanCorr[all_regs])),
    stringsAsFactors = FALSE
  )
  stopifnot(!is.na(out_dt))
  out_dt
}
outFile <- file.path(outFolder, "all_full_partial_delta_dt.Rdata")
save(all_full_partial_delta_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))



all_meanCorr_corrPurity_dt <- aggregate(purityCorr~dataset+region,data=all_ds_corrPurity_dt,FUN=mean)
all_dt <- merge(all_full_partial_delta_dt, all_meanCorr_corrPurity_dt, by =c("dataset", "region"), all=TRUE)
if(purity_ds == "EPIC") stopifnot(!is.na(all_dt))
all_dt$dotCols <- all_cols[all_cmps[basename(as.character(all_dt$dataset))]]


my_x <- all_dt[,c("purityCorr")]

myx_lab <- paste0(transfExpr, " expr. and purity correlation (meanTAD)")

all_vars <- c("absMeanCorr_full_partial_delta", "empPval_full_partial_delta")

for(plot_var in all_vars) {
  
  my_y <- all_dt[,c(plot_var)]
  
  outFile <- file.path(outFolder, paste0(plot_var, "_exprPurityCorr_meanTAD.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_plot(my_x= my_x,
          my_y=my_y,
          xlab=myx_lab,
          main=paste0("meanTAD - ", plot_var, " vs. purity corr."),
          ylab=paste0(plot_var),
          col = all_dt$dotCols
  )
  mtext(side=3, text = paste0(corMet, "'s corr.", " - ", purity_plot_name, " data"))
  legend("bottomleft",all_cols, legend=names(all_cols), bty="n", pch=16, col=all_cols)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}
