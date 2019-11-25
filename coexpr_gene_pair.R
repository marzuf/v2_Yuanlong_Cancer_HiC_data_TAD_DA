########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript coexpr_gene_pair.R\n"))

script_name <- "coexpr_gene_pair.R"


# > similar to figures 4 a) in  Aran et al. 2015

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript coexpr_gene_pair.R
# Rscript coexpr_gene_pair.R <gene1> <gene2> <hicds> <exprds>      # 4
# Rscript coexpr_gene_pair.R EPIC <gene1> <gene2> <hicds> <exprds> # 5

# Rscript coexpr_gene_pair.R <gene1> <gene2>   # 2
# Rscript coexpr_gene_pair.R EPIC <gene1> <gene2> #3

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- TRUE

gene1 <- "GIMAP4"
gene2 <- "GIMAP2"
select_hicds <- "LG1_40kb"
select_exprds <- "TCGAluad_norm_luad"


entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 4 | length(args) == 2) {
  purity_ds <- "" 
  file_suffix <- ""
  gene1 <- args[1]
  gene2 <- args[2]
} else if(length(args) == 5 | length(args) == 3) {
  purity_ds <- args[1]
  file_suffix <- paste0("_", purity_ds)
  gene1 <- args[2]
  gene2 <- args[3]
} else {
  stop("---error args\n")
}
if(length(args) == 4 | length(args) == 5) {
  select_hicds <- args[length(args)-1]
  select_exprds <- args[length(args)]
} else {
  select_hicds <- NA
  select_exprds <- NA
}

stopifnot(gene1 %in% entrez2symb)
stopifnot(gene2 %in% entrez2symb)

gene1_entrez <- names(entrez2symb)[entrez2symb == gene1]
gene2_entrez <- names(entrez2symb)[entrez2symb == gene2]
stopifnot(length(gene1_entrez) == 1)
stopifnot(length(gene2_entrez) == 1)

outFolder <- file.path(paste0("COEXPR_GENE_PAIR", file_suffix))
dir.create(outFolder, recursive = TRUE)

myHeight <- 7
myWidth <- 9
plotType <- "png"

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")


if(purity_ds == "") {
  
  purity_file <- file.path("tcga_purity_aran2015.csv")
  purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
  pm <- purity_metrics[1]
  # all the ranks are between 1 and 0
  
  
} else if(purity_ds == "EPIC") {
  
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


cat(paste0("!!! > purity metric: ", pm, "\n"))



all_hicds <- list.files(pipFolder)
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

stopifnot(is.na(select_exprds) == is.na(select_hicds))

if(!is.na(select_hicds) & !is.na(select_exprds)) {
  all_ds <- all_ds[dirname(all_ds) ==select_hicds & basename(all_ds) ==select_exprds]
}
stopifnot(length(all_ds) > 0)

build_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
  
  settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  
  samp1 <- get(load(file.path(setDir, sample1_file)))
  samp2 <- get(load(file.path(setDir, sample2_file)))
  
  fpkm_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")
  stopifnot(file.exists(fpkm_file))
  fpkm_dt <- get(load(fpkm_file))  
  
  if( ! gene1_entrez %in% rownames(fpkm_dt) & ! gene2_entrez %in% rownames(fpkm_dt)) return(NULL)
  
  
  
  stopifnot(gene1_entrez %in% rownames(fpkm_dt))
  stopifnot(gene2_entrez %in% rownames(fpkm_dt))
  
  pur_samp1 <- samp1[samp1 %in% purity_dt$Sample.ID | paste0(samp1, "A") %in% purity_dt$Sample.ID]
  cat(paste0("For ", cond1, " - available samples:\t", length(pur_samp1), "/", length(samp1), "\n"))
  
  pur_samp2 <- samp2[samp2 %in% purity_dt$Sample.ID | paste0(samp2, "A") %in% purity_dt$Sample.ID]
  cat(paste0("For ", cond2, " - available samples:\t", length(pur_samp2), "/", length(samp2), "\n"))
  
  if( length(pur_samp1) == 0 & length(pur_samp2) == 0) return(NULL)
  
  pur2_samp1 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp1  | purity_dt$Sample.ID %in% paste0(samp1, "A") ]
  stopifnot(length(pur2_samp1) == length(pur_samp1))
  
  pur2_samp2 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp2  | purity_dt$Sample.ID %in% paste0(samp2, "A") ]
  stopifnot(length(pur2_samp2) == length(pur_samp2))
  
  stopifnot(setequal(gsub("A$", "", pur2_samp1), pur_samp1))
  stopifnot(setequal(gsub("A$", "", pur2_samp2), pur_samp2))
  
  purity_values <- setNames(purity_dt[purity_dt$Sample.ID %in% c(pur2_samp1, pur2_samp2), paste0(pm)],
                            purity_dt[purity_dt$Sample.ID %in% c(pur2_samp1, pur2_samp2), paste0("Sample.ID")])
  stopifnot(setequal(c(pur2_samp1, pur2_samp2), names(purity_values)))
  
  names(purity_values) <- gsub("A$", "", names(purity_values))
  stopifnot(setequal(c(pur_samp1, pur_samp2), names(purity_values)))
  purity_values <- purity_values[c(pur_samp1, pur_samp2)]
  
  
  gene1_values <- setNames(unlist(c(fpkm_dt[as.character(gene1_entrez), c(pur_samp1, pur_samp2)])), c(pur_samp1, pur_samp2))
  gene2_values <- setNames(unlist(c(fpkm_dt[as.character(gene2_entrez), c(pur_samp1, pur_samp2)])), c(pur_samp1, pur_samp2))
  
  
  all_samp <- names(purity_values)
  
  stopifnot(length(gene1_values) == length(purity_values))
  stopifnot(length(gene2_values) == length(purity_values))
  
  stopifnot(names(gene1_values) == all_samp)
  stopifnot(names(gene2_values) == all_samp)
  
  data.frame(
    hicds = hicds,
    exprds = exprds,
    samp_id =  all_samp,
    expr_gene1 = as.numeric(gene1_values),
    expr_gene2 = as.numeric(gene2_values),
    samp_purity = as.numeric(purity_values),
    stringsAsFactors = FALSE
  )
  
  
  
}






















coexpr_file <- file.path(paste0("COEXPR_AND_PURITY", file_suffix), "coexpr_and_purity_dt.Rdata")
stopifnot(file.exists(coexpr_file))
coexpr_dt <- get(load(coexpr_file))

corrpurity_file <- file.path(paste0("CORR_EXPR_AND_PURITY", file_suffix), "corr_expr_purity_dt.Rdata")
stopifnot(file.exists(corrpurity_file))
corrpur_dt <- get(load(corrpurity_file))


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

coexpr_corrpur_dt$pairwise_corr_gene_purity <- coexpr_corrpur_dt$gene1_corr_gene_purity * coexpr_corrpur_dt$gene2_corr_gene_purity

coexpr_corrpur_dt$full_partial_corr_diff <- coexpr_corrpur_dt$coexpr - coexpr_corrpur_dt$partial_coexpr

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))
coexpr_corrpur_dt$puritycorr_col <- rbPal(10)[as.numeric(cut(coexpr_corrpur_dt$pairwise_corr_gene_purity,breaks = 10))]



# figure 4
myTit <- ""
mySub <- ""

myylab <- "Pairwise partial correlations"
myxlab <- "Pairwise correlations"

myx <- "coexpr"
myy <- "partial_coexpr"

outFile <- file.path(outFolder, paste0("all_datasets_", myy, "_vs_", myx, ".", plotType))
do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
plot(
  as.formula(paste0(myy, "~", myx)),
  data=coexpr_corrpur_dt,
  xlab =myxlab,
  ylab =myylab,
  main=myTit,
  cex=0.7,
  pch=16,
  col=coexpr_corrpur_dt$puritycorr_col,
  cex.lab=plotCex,
  cex.axis=plotCex
)
curve(1*x, lty=2, col="grey", add=T)
mtext(side=3, text = paste0(mySub))
legend("topleft",legend=levels(cut(coexpr_corrpur_dt$pairwise_corr_gene_purity,breaks = 10)),col =rbPal(10),pch=16, bty="n", cex=0.8)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))    

# figure 5
myxlab <- "Diff. correlation (regular-controlled)"
myylab <- "Pairwise correlation with purity"

myy <- "pairwise_corr_gene_purity"
myx <- "full_partial_corr_diff"

outFile <- file.path(outFolder, paste0("all_datasets_", myy, "_vs_", myx, ".", plotType))
do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
plot(
  as.formula(paste0(myy, "~", myx)),
  data=coexpr_corrpur_dt,
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
legend("topleft",legend=levels(cut(coexpr_corrpur_dt$pairwise_corr_gene_purity,breaks = 10)),col =rbPal(10),pch=16, bty="n", cex=0.8)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))    




##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))




