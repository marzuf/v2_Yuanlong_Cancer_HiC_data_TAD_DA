########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript cmp_pvals_PF.R\n"))

script_name <- "cmp_pvals_PF.R"

# _final -> discussion with Giovanni 04.08.2020 -> take Aran CPE data
# corrected compared to some of the previous versions -> if only non-"A" vial, take the "A" vial
# if multiple vials -> take "A" vials

outFolder <- file.path("CMP_PVALS_PF")
dir.create(outFolder, recursive = TRUE)

# Rscript cmp_pvals_PF.R

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_cols[all_cols=="green"] <- "darkgreen"

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- TRUE

fontFamily <- "Hershey"

myHeightP <- 400 
myWidthP <- 400
plotTypeP <- "png"
plotCex <- 1.4
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


script11pf_name <- "11sameNbrPF_runEmpPvalCombined"
script11_name <- "11sameNbr_runEmpPvalCombined"



mainFolder <- "."

# all_ds=all_ds[1]

ex_hicds <- "ENCSR489OCU_NCI-H460_40kb"
ex_exprds <- "TCGAluad_norm_luad"
ex_TAD <- "chr11_TAD390"

# all_ds = file.path(ex_hicds, ex_exprds)

if(buildTable){
  
  all_pvals_cmp_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
    
    cat(paste0("... start: ", ds, "\n"))
    
    hicds <- dirname(ds)
    exprds <- basename(ds)
    
    pvals_file <- file.path(pipFolder, ds, script11_name, "emp_pval_combined.Rdata")
    pf_pvals_file <- file.path(pipFolder, ds, script11pf_name, "emp_pval_combined.Rdata")
    
    if(!file.exists(pf_pvals_file)) return(NULL)
    
    pvals <- get(load(pvals_file))
    pf_pvals <- get(load(pf_pvals_file))
    
    adj_pvals <- p.adjust(pvals, method="BH")
    adj_pf_pvals <- p.adjust(pf_pvals, method="BH")
    
    stopifnot(names(adj_pf_pvals) %in% names(adj_pvals))
    
    intersectRegs <- intersect(names(adj_pvals), names(adj_pf_pvals))
    
    adj_pvals <- adj_pvals[intersectRegs]
    stopifnot(!is.na(adj_pvals))
    
    adj_pf_pvals <- adj_pf_pvals[intersectRegs]
    stopifnot(!is.na(adj_pf_pvals))
    
    out_dt <- data.frame(
      dataset=ds,
      region=intersectRegs,
      adj_pvals = adj_pvals,
      adj_pf_pvals = adj_pf_pvals,
      stringsAsFactors = FALSE
    )
    rownames(out_dt) <- NULL
    out_dt
  }
  outFile <- file.path(outFolder, "all_pvals_cmp_dt.Rdata")
  save(all_pvals_cmp_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_pvals_cmp_dt.Rdata")
  all_pvals_cmp_dt <- get(load(outFile))
}

nDS <- length(unique(all_pvals_cmp_dt$dataset))
plot_var <- "adjPvals_log10"
my_x <- -log10(all_pvals_cmp_dt[,c("adj_pvals")])
my_y <- -log10(all_pvals_cmp_dt[,c("adj_pf_pvals")])


outFile <- file.path(outFolder, paste0(plot_var,"_PF_vs_all_densplot.", plotTypeP)) 
do.call(plotTypeP, list(outFile, height=myHeightP, width=myWidthP))
densplot(x=my_x, y=my_y,
         main=paste0(plot_var),
         xlab=paste0("adj. pvals all TADs [-log10]"), 
         ylab=paste0("adj. pvals PF TADs [-log10]"),
         cex.main=plotCex,cex.lab=plotCex,cex.axis=plotCex,
         cex=0.7, pch=16)
mtext(side=3, text=paste0("# DS = ", nDS))
curve(1*x, lty=1, col="darkgrey", add=TRUE)
addCorr(x=my_x,y=my_y,bty="n", legPos = "topleft")

foo <- dev.off()
cat(paste0("...written:", outFile,"\n"))

