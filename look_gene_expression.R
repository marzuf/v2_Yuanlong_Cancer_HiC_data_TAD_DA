########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript look_gene_expression.R\n"))

script_name <- "look_gene_expression.R"


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript look_gene_expression.R <hicds> <exprds> <gene_symbol>
# Rscript look_gene_expression.R LG1_40kb TCGAlusc_norm_lusc GIMAP2
# Rscript look_gene_expression.R LG1_40kb TCGAlusc_norm_lusc C1QA

hicds="LG1_40kb"
exprds="TCGAlusc_norm_lusc"
gene_symbol="GIMAP2"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
hicds <- args[1]
exprds <- args[2]
all_gene_symbols <- args[3:length(args)]
init_order <- all_gene_symbols


# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4
myHeightGG <- 7
myWidthGG <- 7

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")


outFolder <- file.path("LOOK_GENE_EXPRESSION")
dir.create(outFolder, recursive = TRUE)

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

all_gene_symbols_0 <- all_gene_symbols
all_gene_symbols <- all_gene_symbols[all_gene_symbols %in% entrez2symb]
if(length(all_gene_symbols) < length(all_gene_symbols_0)) {
  cat(paste0("!!! WARNING: following genes discarded (missing von gff_dt):\t", paste0(all_gene_symbols_0[!all_gene_symbols_0 %in% all_gene_symbols], collapse=","), "\n"))
}


stopifnot(all_gene_symbols %in% entrez2symb)
all_gene_entrez <- names(entrez2symb)[entrez2symb %in% all_gene_symbols]
stopifnot(length(all_gene_entrez) == length(all_gene_symbols))

########################################################################################################################################################################################



settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)


samp1 <- get(load(file.path(setDir, sample1_file)))
samp2 <- get(load(file.path(setDir, sample2_file)))


fpkm_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")
stopifnot(file.exists(fpkm_file))
fpkm_dt <- get(load(fpkm_file))

all_gene_entrez_0 <- all_gene_entrez

all_gene_entrez <- all_gene_entrez[all_gene_entrez %in% rownames(fpkm_dt)]
if(length(all_gene_entrez) < length(all_gene_entrez_0)) {
  missing_entrez <- all_gene_entrez_0[!all_gene_entrez_0 %in% all_gene_entrez]
  cat(paste0("!!! WARNING: following genes discarded (missing von fpkm_dt):\t", paste0(entrez2symb[paste0(missing_entrez)], collapse=","), "\n"))
  all_gene_symbols <- entrez2symb[paste0(all_gene_entrez)]
}

stopifnot(all_gene_entrez %in% rownames(fpkm_dt))
stopifnot(samp1 %in% colnames(fpkm_dt))
stopifnot(samp2 %in% colnames(fpkm_dt))

plotTit <- paste0(hicds, " - ", exprds)
plotSubTit <- paste0("~ ", paste0(all_gene_symbols, collapse=","), " ~")
myylab <- "RNA-seq scaled estimates"

if(length(all_gene_symbols) == 0) {
  stop(cat(paste0("... no data found for this gene\n")))
} else if(length(all_gene_symbols) == 1) {
  
  
  cond1_fpkm <- unlist(fpkm_dt[paste0(all_gene_entrez), samp1])
  cond2_fpkm <- unlist(fpkm_dt[paste0(all_gene_entrez), samp2])
  
  boxplot_dt <- data.frame(
    FPKM = c(cond1_fpkm, cond2_fpkm),
    cond = c(rep(cond1, length(samp1)), rep(cond2, length(samp2))),
    stringsAsFactors = FALSE
  )
  
  outFile <- file.path(outFolder, paste0(all_gene_symbols, "_", hicds, "_", exprds, ".", plotType))
  do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
  par(bty="l")
  boxplot(FPKM~cond, data=boxplot_dt,
          main = plotTit,
          ylab=myylab,
          xlab="",
          cex.axis = plotCex,
          cex.lab=plotCex)
  mtext(side=3, text = plotSubTit, cex=1.4)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))    
  
} else {

  count_dt <- fpkm_dt[paste0(all_gene_entrez), c(samp1, samp2)]
  
  count_dt_t <- as.data.frame(t(count_dt))
  stopifnot(colnames(count_dt_t) %in% names( entrez2symb))
  colnames(count_dt_t) <- entrez2symb[colnames(count_dt_t)]
  
  stopifnot(all_gene_symbols %in% colnames(count_dt_t))
  
  gene_order <- init_order[init_order %in% all_gene_symbols]
  stopifnot(length(gene_order) == length(all_gene_symbols))
  
  count_dt_t$variable <- rownames(count_dt_t)
  count_dt_t$cond <- ifelse(count_dt_t$variable %in% samp1, cond1, 
         ifelse(count_dt_t$variable %in% samp2, cond2, NA))
  
  
  expr_p <- ggboxplot(count_dt_t, x = "cond",
            xlab = "",
            y = gene_order,
            combine = TRUE,
            ylab = myylab,
            color = "cond", palette = "jco") + 
    guides(color=FALSE)+
    # ggtitle(plotTit, sub=plotSubTit)+
    ggtitle(plotTit)+
    theme(
      plot.title = element_text(size=16, hjust=0.5, face = "bold"),
      plot.subtitle = element_text(size=14, hjust=0.5),
      strip.text.x =  element_text(size=12),
      axis.text.x = element_text(size=12),
      axis.title.y = element_text(size=14)
    )
  
  
  outFile <- file.path(outFolder, paste0(paste0(all_gene_symbols, collapse="_"), "_", hicds, "_", exprds, ".", plotType))
  ggsave(filename = outFile, plot = expr_p, height=myHeightGG, width=myWidthGG*length(all_gene_symbols)*0.8)
  cat(paste0("... written: ", outFile, "\n"))    
  
  
}





##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

