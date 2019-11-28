########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript pairwise_corr_and_purity_byDS.R\n"))

script_name <- "pairwise_corr_and_purity_byDS.R"


# > similar to figures 4 c) and d) in  Aran et al. 2015

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

#  
# Rscript pairwise_corr_and_purity_byDS.R EPIC

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- TRUE


args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  purity_ds <- "" 
  file_suffix <- ""
} else {
  stopifnot(length(args) == 1)
  purity_ds <- args[1]
  file_suffix <- paste0("_", purity_ds)
}


outFolder <- file.path(paste0("PAIRWISE_CORR_AND_PURITY_BYDS", file_suffix))
dir.create(outFolder, recursive = TRUE)

myHeight <- 7
myWidth <- 9
plotType <- "png"

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

if(purity_ds == "") {
  file_suffix <- ""
} else if(purity_ds == "EPIC") {
  file_suffix <- "_EPIC"
} else{
  stop("---invalid DS\n")
}

coexpr_folder <- file.path(paste0("COEXPR_AND_PURITY_BYDS", file_suffix))
stopifnot(dir.exists(coexpr_folder))

corrpurity_file <- file.path(paste0("CORR_EXPR_AND_PURITY", file_suffix), "corr_expr_purity_dt.Rdata")
stopifnot(file.exists(corrpurity_file))
cat(paste0("load correxpr dt - ", Sys.time(), " - "))
corrpur_dt <- get(load(corrpurity_file))
cat(paste0(Sys.time(), " \n "))


all_hicds <- list.files(pipFolder)
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]


foo <- foreach(ds = all_ds ) %dopar% {

  hicds <- dirname(ds)
  exprds <- basename(ds)
  
  cat(paste0("... start ", hicds, " - ", exprds, "\n"))  
  
  coexpr_corrpur_dt_file <- file.path(coexpr_folder, paste0(hicds, "_", exprds, "_", "coexpr_and_purity_dt.Rdata" ))
  
  if(!file.exists(coexpr_corrpur_dt_file)) return(NULL)
  
  stopifnot(file.exists(coexpr_corrpur_dt_file))
  
  cat(paste0("load coexpr dt - ", Sys.time(), " - "))
  coexpr_dt <- get(load(coexpr_corrpur_dt_file))
  cat(paste0(Sys.time(), " \n "))
  
  coexpr_dt_full <- coexpr_dt
  corrpur_dt_full <- corrpur_dt
  
  coexpr_dt <- coexpr_dt[,c("hicds", "exprds", "gene1", "gene2", "coexpr", "partial_coexpr")]
  corrpur_dt <- corrpur_dt[,c("hicds", "exprds", "entrezID", "all_corr_gene_purity")]
  
  # merge for gene1
  colnames(corrpur_dt)[colnames(corrpur_dt) == "entrezID"] <- "gene1"
  
  coexpr_corrpur_dt_1 <- merge(coexpr_dt, corrpur_dt, by=c("hicds", "exprds", "gene1"), all.x=TRUE, all.y=FALSE)
  colnames(coexpr_corrpur_dt_1)[colnames(coexpr_corrpur_dt_1) == "all_corr_gene_purity"] <- "gene1_corr_gene_purity"
  
  # merge for gene2
  colnames(corrpur_dt)[colnames(corrpur_dt) == "gene1"] <- "gene2"
  coexpr_corrpur_dt <- merge(coexpr_corrpur_dt_1, corrpur_dt, by=c("hicds", "exprds", "gene2"), all.x=TRUE, all.y=FALSE)
  colnames(coexpr_corrpur_dt)[colnames(coexpr_corrpur_dt) == "all_corr_gene_purity"] <- "gene2_corr_gene_purity"
  
  sub_coexpr_corrpur_dt <- coexpr_corrpur_dt
  
  sub_coexpr_corrpur_dt$pairwise_corr_gene_purity <- sub_coexpr_corrpur_dt$gene1_corr_gene_purity * sub_coexpr_corrpur_dt$gene2_corr_gene_purity
  
  sub_coexpr_corrpur_dt$full_partial_corr_diff <- sub_coexpr_corrpur_dt$coexpr - sub_coexpr_corrpur_dt$partial_coexpr
  
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  sub_coexpr_corrpur_dt$puritycorr_col <- rbPal(10)[as.numeric(cut(sub_coexpr_corrpur_dt$pairwise_corr_gene_purity,breaks = 10))]
  
  
  
  # figure 4
  myTit <- paste0(ds)
  mySub <- ""
  
  myylab <- "Pairwise partial correlations"
  myxlab <- "Pairwise correlations"
  
  myx <- "coexpr"
  myy <- "partial_coexpr"
  
  # outFile <- file.path(outFolder, paste0("all_datasets_", myy, "_vs_", myx, ".", plotType))
  outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", myy, "_vs_", myx, ".", plotType))
  do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
  plot(
    as.formula(paste0(myy, "~", myx)),
    data=sub_coexpr_corrpur_dt,
    xlab =myxlab,
    ylab =myylab,
    main=myTit,
    cex=0.7,
    pch=16,
    col=sub_coexpr_corrpur_dt$puritycorr_col,
    cex.lab=plotCex,
    cex.axis=plotCex
  )
  curve(1*x, lty=2, col="grey", add=T)
  mtext(side=3, text = paste0(mySub))
  legend("topleft",legend=levels(cut(sub_coexpr_corrpur_dt$pairwise_corr_gene_purity,breaks = 10)),col =rbPal(10),pch=16, bty="n", cex=0.8)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))    
  
  # figure 5
  myxlab <- "Diff. correlation (regular-controlled)"
  myylab <- "Pairwise correlation with purity"
  
  myy <- "pairwise_corr_gene_purity"
  myx <- "full_partial_corr_diff"
  
  # outFile <- file.path(outFolder, paste0("all_datasets_", myy, "_vs_", myx, ".", plotType))
  outFile <- file.path(outFolder, paste0(hicds, "_", exprds,"_", myy, "_vs_", myx, ".", plotType))
  do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
  plot(
    as.formula(paste0(myy, "~", myx)),
    data=sub_coexpr_corrpur_dt,
    xlab =myxlab,
    ylab =myylab,
    main=myTit,
    cex=0.7,
    pch=16,
    col="black",
    cex.lab=plotCex,
    cex.axis=plotCex
  )
  curve(1*x, lty=2, col="grey", add=T)
  mtext(side=3, text = paste0(mySub))
  legend("topleft",legend=levels(cut(sub_coexpr_corrpur_dt$pairwise_corr_gene_purity,breaks = 10)),col =rbPal(10),pch=16, bty="n", cex=0.8)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))    
  
  
} # end.iterating over ds

##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))




  