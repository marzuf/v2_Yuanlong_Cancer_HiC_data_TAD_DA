########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript look_gene_expression_level.R\n"))

script_name <- "look_gene_expression_level.R"


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript look_gene_expression_level.R <hicds> <exprds> <gene_symbol>
# Rscript look_gene_expression_level.R LG1_40kb TCGAlusc_norm_lusc GIMAP2
# Rscript look_gene_expression_level.R LG1_40kb TCGAlusc_norm_lusc C1QA
# 
# 
# Rscript look_gene_expression_level.R LG1_40kb TCGAlusc_norm_lusc GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R ENCSR079VIJ_40kb TCGAlihc_norm_lihc GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R LG2_40kb TCGAluad_norm_luad GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R PA3_40kb TCGAskcm_wt_mutKRAS GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R LG2_40kb TCGAluad_wt_mutKRAS GIMAP2 GIMAP4 GIMAP6
# Rscript look_gene_expression_level.R K562_40kb TCGAlaml_wt_mutFLT3 GIMAP2 GIMAP4 GIMAP6
# 
# 
# Rscript look_gene_expression_level.R ENCSR444WCZ_A549_40kb TCGAluad_norm_luad MMP1 MMP12
# Rscript look_gene_expression_level.R ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc MMP1 MMP12
# Rscript look_gene_expression_level.R ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss MMP1 MMP12
# 
# 
# Rscript look_gene_expression_level.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF CD244 CD48 LY9 SLAMF7
# Rscript look_gene_expression_level.R ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS CD244 CD48 LY9 SLAMF7
# Rscript look_gene_expression_level.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS CD244 CD48 LY9 SLAMF7
# Rscript look_gene_expression_level.R GSE105381_HepG2_40kb TCGAlihc_norm_lihc CD244 CD48 LY9 SLAMF7
# 

# Rscript look_gene_expression_level.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF LY9 SLAMF7
# Rscript look_gene_expression_level.R ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS LY9 SLAMF7
# Rscript look_gene_expression_level.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS LY9 SLAMF7
# Rscript look_gene_expression_level.R GSE105381_HepG2_40kb TCGAlihc_norm_lihc LY9 SLAMF7

# Rscript look_gene_expression_level.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF CD244 CD48 
# Rscript look_gene_expression_level.R ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS CD244 CD48 
# Rscript look_gene_expression_level.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS CD244 CD48 
# Rscript look_gene_expression_level.R GSE105381_HepG2_40kb TCGAlihc_norm_lihc CD244 CD48 
# 

# 

# Rscript look_gene_expression_level.R LG1_40kb TCGAluad_norm_luad C1QA C1QB C1QC
# Rscript look_gene_expression_level.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf C1QA C1QB C1QC
# Rscript look_gene_expression_level.R GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal C1QA C1QB C1QC
# Rscript look_gene_expression_level.R GSE105194_cerebellum_40kb TCGAgbm_classical_neural C1QA C1QB C1QC
# Rscript look_gene_expression_level.R GSE105318_DLD1_40kb TCGAcoad_msi_mss C1QA C1QB C1QC
# Rscript look_gene_expression_level.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF C1QA C1QB C1QC
# 
# 
# Rscript look_gene_expression_level.R GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural HOXC4 HOXC6
# Rscript look_gene_expression_level.R ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad HOXC4 HOXC6
# Rscript look_gene_expression_level.R GSE118514_22Rv1_40kb TCGAprad_norm_prad HOXC4 HOXC6
# 
# 




hicds="LG1_40kb"
exprds="TCGAlusc_norm_lusc"
gene_symbol="GIMAP2"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
hicds <- args[1]
exprds <- args[2]
all_gene_symbols <- args[3:length(args)]
all_gene_symbols <- unique(all_gene_symbols)
init_order <- all_gene_symbols


# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
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


outFolder <- file.path("LOOK_GENE_EXPRESSION_LEVEL")
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
} else {
  
  rowsOfInterest <- rownames(fpkm_dt) %in% all_gene_entrez
  
  medianExpr <- apply(fpkm_dt, 1, median)
  
  ofinterest_med <- medianExpr[rowsOfInterest]
  other_med <-  medianExpr[! rowsOfInterest]
  
  ofinterest_med_noout <- ofinterest_med[ofinterest_med>= quantile(ofinterest_med, probs=0.05) & ofinterest_med <= quantile(ofinterest_med, probs=0.95)]
  other_med_noout <- other_med[other_med>= quantile(other_med, probs=0.05) & other_med <= quantile(other_med, probs=0.95)]
  
  outFile <- file.path(outFolder, paste0(paste0(all_gene_symbols, collapse="_"), "_", hicds, "_", exprds, "_meds_density_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
  plot(density(log10(other_med_noout)), main=plotTit)
  abline(v=ofinterest_med, lty=2, col="red")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))    
  
  ofinterest_dt <- fpkm_dt[rowsOfInterest,]
  other_dt <-  fpkm_dt[! rowsOfInterest,]
  
  interest_meds_samp1 <- apply(ofinterest_dt[,c(samp1)], 1, median)
  other_meds_samp1 <- apply(ofinterest_dt[,c(samp1)], 1, median)
  interest_meds_samp2<- apply(ofinterest_dt[,c(samp2)], 1, median)
  other_meds_samp2 <- apply(ofinterest_dt[,c(samp2)], 1, median)
  
  outFile <- file.path(outFolder, paste0(paste0(all_gene_symbols, collapse="_"),"_", hicds, "_", exprds, "_samp1_samp2_meds.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="l")
  plot(
    x=interest_meds_samp1,
    y = interest_meds_samp2,
    pch=16,
    cex=0.7,
    col="black",
    cex.axis=plotCex,
    cex.lab=plotCex,
    xlab="med expr. samp1",
    ylab="med expr. samp2",
    main=plotTit)
  points(x=other_meds_samp1,
         y=other_meds_samp2,
         col="red",
         pch=16,
         cex=0.7)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  stopifnot(nrow(ofinterest_dt) + nrow(other_dt) == nrow(fpkm_dt))
  
  cond1_fpkm <- unlist(ofinterest_dt[, samp1])
  cond2_fpkm <- unlist(ofinterest_dt[, samp2])
  all_other_fpkm <- unlist(other_dt[, c(samp1, samp2)])
  
  boxplot_dt <- data.frame(
    FPKM = c(cond1_fpkm, cond2_fpkm, all_other_fpkm),
    cond = c(rep(cond1, length(cond1_fpkm)), rep(cond2, length(cond2_fpkm)), rep("all_other", length(all_other_fpkm))),
    stringsAsFactors = FALSE
  )
  boxplot_dt$FPKM_log10 <- log10(boxplot_dt$FPKM)
  
  outFile <- file.path(outFolder, paste0(paste0(all_gene_symbols, collapse="_"), "_", hicds, "_", exprds, "_all_counts_log10.", plotType))
  p_all <- ggdensity(
      data=boxplot_dt,
      x = "FPKM_log10",
      color="cond",
      title = plotTit, subtitle=plotSubTit
  )
  ggsave(p_all, filename = outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  # ggboxplot(
  #   data=boxplot_dt,
  #   x="cond",
  #   y = "FPKM_log10",
  #   color="cond"
  # )
  # 
  # mtext(side=3, text = plotSubTit, cex=1.4)
  # foo <- dev.off()
  # cat(paste0("... written: ", outFile, "\n"))    
  
} 



##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
