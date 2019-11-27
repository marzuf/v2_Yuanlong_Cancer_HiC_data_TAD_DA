########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript coexpr_expr_and_purity_byDS.R\n"))

script_name <- "cexpr_expr_and_purity_byDS.R"

source("pcor.R")

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(psych, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

cor_na_method <- "complete.obs"

get_full_corr_dt <-  function(fpkmdt, cormet, newcol) {
  
  fullcorr <- cor(t(fpkmdt), use=cor_na_method)
  stopifnot(isSymmetric(fullcorr))
  fullcorr[lower.tri(fullcorr, diag = TRUE)] <- NA
  fullcorr_dt <- melt(fullcorr, na.rm = TRUE)
  # fullcorr_dt$gene1 <- as.character(pmin(fullcorr_dt$Var1, fullcorr_dt$Var2))
  fullcorr_dt$gene1 <- fullcorr_dt$Var1
  fullcorr_dt$gene2 <- fullcorr_dt$Var2
  # fullcorr_dt$gene2 <- as.character(pmax(fullcorr_dt$Var1, fullcorr_dt$Var2))
  colnames(fullcorr_dt)[colnames(fullcorr_dt) == "value"] <- newcol
  fullcorr_dt <- fullcorr_dt[,c("gene1", "gene2", newcol)]
  return(fullcorr_dt)
}

get_partial_corr_dt <- function(fpkmdt_with_pur, cormet, newcol) {
  # partial correlations for a set (x) of variables with set (y) removed.
  partialcorr <- as(partial.r(t(fpkmdt_with_pur),
                              x = 1:(nrow(fpkmdt_with_pur)-1),
                              y=nrow(fpkmdt_with_pur),
                              method=cormet,
                              use=cor_na_method), "matrix")
  stopifnot(isSymmetric(partialcorr))
  partialcorr[lower.tri(partialcorr, diag = TRUE)] <- NA
  partialcorr_dt <- melt(partialcorr, na.rm=TRUE)
  # partialcorr_dt$gene1 <- as.character(pmin(partialcorr_dt$Var1, partialcorr_dt$Var2))
  partialcorr_dt$gene1 <- partialcorr_dt$Var1
  # partialcorr_dt$gene2 <- as.character(pmax(partialcorr_dt$Var1, partialcorr_dt$Var2))
  partialcorr_dt$gene2 <- partialcorr_dt$Var2
  colnames(partialcorr_dt)[colnames(partialcorr_dt) == "value"] <- newcol
  partialcorr_dt <- partialcorr_dt[,c("gene1", "gene2", newcol)]
  return(partialcorr_dt)
}

# Rscript coexpr_and_purity_byDS.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- FALSE

myHeight <- 400 
myWidth <- 400
plotType <- "png"
plotCex <- 1.4
mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
} else{
  purity_ds <- ""
}


if(purity_ds == "") {
  file_suffix <- ""
  purity_file <- file.path("tcga_purity_aran2015.csv")
  purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
  pm <- purity_metrics[1]
  # all the ranks are between 1 and 0
  
  
} else if(purity_ds == "EPIC") {
  file_suffix <- "_EPIC"
  purity_file <- file.path("EPIC_PURITY/all_epic_purity_data.Rdata")
  epic_purity_data <- get(load(purity_file))
  purity_dt <- as.data.frame(do.call(rbind, c(lapply(epic_purity_data, function(x) x[["infiltration_fraction"]]))))
  all_pm_metrics <- colnames(purity_dt) #"Bcells"      "CAFs"        "CD4_Tcells"  "CD8_Tcells"  "Endothelial" "Macrophages" "NKcells"     "otherCells" 
  pm <- "otherCells"
  purity_dt$Sample.ID <- rownames(purity_dt)
  purity_dt$Sample.ID <- gsub("\\.", "-", purity_dt$Sample.ID)
  
  
} else {
  stop("--invalid purity_ds\n")
}

outFolder <- file.path(paste0("COEXPR_AND_PURITY_BYDS", file_suffix))
dir.create(outFolder, recursive = TRUE)


cat(paste0("!!! > purity metric: ", pm, "\n"))



all_hicds <- list.files(pipFolder)
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

corMet <- "pearson"

cat(paste0("!!! HARD-CODED !!!\n"))
cat(paste0(">>> corMet\t=\t", corMet, "\n"))
cat(paste0(">>> purity metric\t=\t", pm, "\n"))

if(buildTable) {
  
  foo <- foreach(ds = all_ds, .combine='rbind') %dopar% {
    
    hicds <- dirname(ds)
    exprds <- basename(ds)
    
    settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingFile))
    source(settingFile)
    
    samp1 <- get(load(file.path(setDir, sample1_file)))
    samp2 <- get(load(file.path(setDir, sample2_file)))
    
    pur_samp1 <- samp1[samp1 %in% purity_dt$Sample.ID | paste0(samp1, "A") %in% purity_dt$Sample.ID]
    cat(paste0("For ", cond1, " - available samples:\t", length(pur_samp1), "/", length(samp1), "\n"))
    
    pur_samp2 <- samp2[samp2 %in% purity_dt$Sample.ID | paste0(samp2, "A") %in% purity_dt$Sample.ID]
    cat(paste0("For ", cond2, " - available samples:\t", length(pur_samp2), "/", length(samp2), "\n"))
    
    if(length(pur_samp1) == 0 & length(pur_samp2) == 0) {
      
      ds_coexpr_dt <- data.frame(
        hicds=hicds,
        exprds=exprds,
        gene1=NA,
        gene2=NA,
        coexpr=NA,
        coexpr_samp1=NA,
        coexpr_samp2=NA,
        partial_coexpr=NA,
        partial_coexpr_samp1=NA,
        partial_coexpr_samp2=NA,
        stringsAsFactors = FALSE
      )

  outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_and_purity_dt.Rdata"))
  save(ds_coexpr_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
	
	return(NULL)
    }
    
    pur2_samp1 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp1  | purity_dt$Sample.ID %in% paste0(samp1, "A") ]
    stopifnot(length(pur2_samp1) == length(pur_samp1))
    
    pur2_samp2 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp2  | purity_dt$Sample.ID %in% paste0(samp2, "A") ]
    stopifnot(length(pur2_samp2) == length(pur_samp2))
    
    stopifnot(setequal(gsub("A$", "", pur2_samp1), pur_samp1))
    stopifnot(setequal(gsub("A$", "", pur2_samp2), pur_samp2))
    
    fpkm_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")
    stopifnot(file.exists(fpkm_file))
    fpkm_dt <- get(load(fpkm_file))  
    
    
    
    gene_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")
    stopifnot(file.exists(gene_file))
    geneList <- get(load(gene_file))  
    stopifnot(names(geneList) %in% rownames(fpkm_dt))
    fpkm_dt <- fpkm_dt[rownames(fpkm_dt) %in% names(geneList),]
    newnames <- sapply(rownames(fpkm_dt), function(x) geneList[x])
    stopifnot(!duplicated(newnames))
    stopifnot(length(newnames) == nrow(fpkm_dt))
    rownames(fpkm_dt) <- as.character(newnames)
    
    
    stopifnot(pur_samp1 %in% colnames(fpkm_dt))
    stopifnot(pur_samp2 %in% colnames(fpkm_dt))
    stopifnot(pur2_samp1 %in% purity_dt$Sample.ID)
    stopifnot(pur2_samp2 %in% purity_dt$Sample.ID)
    
    purity_values <- setNames(purity_dt[purity_dt$Sample.ID %in% pur2_samp1 | purity_dt$Sample.ID %in% pur2_samp2,paste0(pm)],
                              purity_dt[purity_dt$Sample.ID %in% pur2_samp1 | purity_dt$Sample.ID %in% pur2_samp2,paste0("Sample.ID")])
    
    names(purity_values) <- gsub("A$", "", names(purity_values))
    pur2_samp1 <- gsub("A$", "", pur2_samp1)
    pur2_samp2 <- gsub("A$", "", pur2_samp2)
    stopifnot(setequal(names(purity_values), c(pur2_samp1, pur2_samp2)))
    
    stopifnot(setequal(names(purity_values), c(pur_samp1, pur_samp2)))
    stopifnot(length(purity_values) == length(c(pur_samp1, pur_samp2)))
    
    purity_values <- purity_values[c(pur2_samp1, pur2_samp2)]
    fpkm_dt <- fpkm_dt[,c(pur_samp1, pur_samp2)]
    
    # updated 27.11.2019 -> avoid pmin and pmax in the function
    tmprownames <- as.character(sort(as.numeric(rownames(fpkm_dt))))
    tmpdim <- dim(fpkm_dt)
    stopifnot(setequal(rownames(fpkm_dt), tmprownames))
    fpkm_dt <- fpkm_dt[paste0(tmprownames),]
    stopifnot(dim(fpkm_dt) == tmpdim)
    
    fpkm_purity_dt <- fpkm_dt
    fpkm_purity_dt[(nrow(fpkm_purity_dt) + 1),] <- purity_values
    rownames(fpkm_purity_dt)[nrow(fpkm_purity_dt)] <- "purity"
    stopifnot(rownames(fpkm_purity_dt)[nrow(fpkm_purity_dt)] == "purity")
    
    cat(paste0("... ", hicds, " - ", exprds, " \t ", "corr. for all samples\n"))
    
    all_partial_corr_dt <- get_partial_corr_dt(fpkmdt_with_pur=fpkm_purity_dt, cormet=corMet, newcol="partial_coexpr")
    all_corr_dt <- get_full_corr_dt(fpkmdt=fpkm_dt, cormet=corMet, newcol="coexpr")
    
    all_dt <- merge(all_partial_corr_dt, all_corr_dt, by=c("gene1", "gene2")) 
    
    cat(paste0("... ", hicds, " - ", exprds, " \t ", "corr. for samp1\n"))    
   if(length(pur_samp1) > 1) {
     samp1_partial_corr_dt <- get_partial_corr_dt(fpkmdt_with_pur=fpkm_purity_dt[,c(pur_samp1)], cormet=corMet, newcol="partial_coexpr_samp1")
     samp1_corr_dt <- get_full_corr_dt(fpkmdt=fpkm_dt[,c(pur_samp1)], cormet=corMet, newcol="coexpr_samp1")
     samp1_dt <- merge(samp1_partial_corr_dt, samp1_corr_dt, by=c("gene1", "gene2")) 
     out_dt <- merge(all_dt, samp1_dt, by=c("gene1", "gene2")) 
   } else {
     out_dt <- all_dt
     out_dt$partial_coexpr_samp1 <- out_dt$coexpr_samp1 <- NA
   }
    
    cat(paste0("... ", hicds, " - ", exprds, " \t ", "corr. for samp2\n"))
    if(length(pur_samp2) > 1) {
      samp2_partial_corr_dt <- get_partial_corr_dt(fpkmdt_with_pur=fpkm_purity_dt[,c(pur_samp2)], cormet=corMet, newcol="partial_coexpr_samp2")
      samp2_corr_dt <- get_full_corr_dt(fpkmdt=fpkm_dt[,c(pur_samp2)], cormet=corMet, newcol="coexpr_samp2")
      samp2_dt <- merge(samp2_partial_corr_dt, samp2_corr_dt, by=c("gene1", "gene2")) 
      out_dt <- merge(out_dt, samp2_dt, by=c("gene1", "gene2")) 
    } else {
      out_dt$partial_coexpr_samp2 <- out_dt$coexpr_samp2 <- NA
    }

  out_dt$hicds <- hicds
  out_dt$exprds <- exprds
      
   ds_coexpr_dt <- out_dt

  outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_coexpr_and_purity_dt.Rdata"))
  save(ds_coexpr_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))

	return(NULL)

  } # end-iterating ds 
    
  

  
  
} else {
cat(paste0("------ done by DS now\n"))
stop("---\n")
  outFile <- file.path(outFolder, "coexpr_and_purity_dt.Rdata")
  # outFile <- file.path(outFolder, "sub_coexpr_and_purity_dt.Rdata")
  cat(paste0("... load coexpr data - ", Sys.time()))
  coexpr_and_purity_dt <- get(load(outFile))
  cat(paste0(" - ", Sys.time(), " - done\n"))
}

# 
# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# 
# all_suffix <- c("", "_samp1", "_samp2")
# 
# suff=""
# for(suff in all_suffix){
#   
#   outFile <- file.path(outFolder, paste0("all_datasets_partial_vs_full_coexpr", suff, ".", plotType))
#   do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
#   densplot(
#     x = coexpr_and_purity_dt[,paste0("coexpr", suff)],
#     y = coexpr_and_purity_dt[,paste0("partial_coexpr", suff)],
# 	xlab = paste0("coexpr", suff),
# 	ylab = paste0("partial_coexpr", suff),
#     main=paste0("partial vs. full coexpr."),
#     cex.axis=plotCex,
#     cex.lab=plotCex,
#     pch=16,
#     cex=0.7
#   )  
#   legend("topleft", legend = paste0("n=", nrow(coexpr_and_purity_dt)), bty="n")
#   
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))    
#   
#   
# }


##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
