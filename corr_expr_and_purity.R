########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript corr_expr_and_purity.R\n"))

script_name <- "corr_expr_and_purity.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript corr_expr_and_purity.R 
# Rscript corr_expr_and_purity.R EPIC

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- TRUE

myHeightGG <- 7
myWidthGG <- 9
plotType <- "png"

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
  purity_metrics <- colnames(purity_dt) #"Bcells"      "CAFs"        "CD4_Tcells"  "CD8_Tcells"  "Endothelial" "Macrophages" "NKcells"     "otherCells" 
  pm <- "otherCells"
  purity_dt$Sample.ID <- rownames(purity_dt)
  purity_dt$Sample.ID <- gsub("\\.", "-", purity_dt$Sample.ID)
} else {
  stop("--invalid purity_ds\n")
}


outFolder <- file.path(paste0("CORR_EXPR_AND_PURITY", file_suffix))
dir.create(outFolder, recursive = TRUE)


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
  
  corr_expr_purity_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
    
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
      
      return(data.frame(
        hicds=hicds,
        exprds=exprds,
        entrezID = NA,
        all_corr_gene_purity = NA,
        samp1_corr_gene_purity = NA,
        samp2_corr_gene_purity = NA,
        stringsAsFactors = FALSE
      ))
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
    
    stopifnot(length(purity_values) == ncol(fpkm_dt))
    
    all_puritygene_corr <- apply(fpkm_dt, 1, function(x) cor.test(x=x, y=purity_values, method=corMet, na.rm=TRUE)$estimate)
    if(length(pur_samp1) == 0) {
      samp1_puritygene_corr <- NA
    } else {
      samp1_puritygene_corr <- apply(fpkm_dt[,c(pur_samp1)], 1, function(x) cor.test(x=x, y=purity_values[pur2_samp1], method=corMet, na.rm=TRUE)$estimate)   
    }
    if(length(pur_samp2) == 0) {
      samp2_puritygene_corr <- NA
    } else {
      samp2_puritygene_corr <- apply(fpkm_dt[,c(pur_samp2)], 1, function(x) cor.test(x=x, y=purity_values[pur2_samp2], method=corMet, na.rm=TRUE)$estimate)  
    }
    
    out_genes <- names(all_puritygene_corr)
    stopifnot(out_genes == names(samp1_puritygene_corr))
    stopifnot(out_genes == names(samp2_puritygene_corr))
    
    
    data.frame(
      hicds=hicds,
      exprds=exprds,
      entrezID = out_genes,
      all_corr_gene_purity = as.numeric(all_puritygene_corr[out_genes]),
      samp1_corr_gene_purity = as.numeric(samp1_puritygene_corr[out_genes]),
      samp2_corr_gene_purity = as.numeric(samp2_puritygene_corr[out_genes]),
      stringsAsFactors = FALSE
    )
    
  }
  outFile <- file.path(outFolder, "corr_expr_purity_dt.Rdata")
  save(corr_expr_purity_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {
  outFile <- file.path(outFolder, "corr_expr_purity_dt.Rdata")
  corr_expr_purity_dt <- get(load(outFile))
  
}

plot_dt <- melt(corr_expr_purity_dt, id=c("hicds", "exprds", "entrezID"))
  
plotTit <- paste0("Corr. expression and purity")
subTit <- paste0("corr. method = ", corMet)
  
myxlab <- ""
myylab <- paste0("purity (", pm, ")")

expr_p_boxplot <- ggboxplot(plot_dt, x = "variable",
                    y = "value",
                      xlab = myxlab,
                      # scales='free_y',
                      # y = gene_order,
                      # combine = TRUE,
                      ylab = myylab,
                       palette = "jco") + 
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
  

outFile <- file.path(outFolder, paste0("correlation_expression_purity_boxplot.", plotType))
ggsave(filename = outFile, plot = expr_p_boxplot, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))    


expr_p_density <- ggdensity(plot_dt,
          title = paste0(plotTit),
          subtitle = paste0(subTit),
          x = "value", 
          col = "variable",
          palette="jco") + 
  labs(color="")+
  theme(
    plot.title = element_text(size=16, hjust=0.5, face = "bold"),
    plot.subtitle = element_text(size=14, hjust=0.5),
    strip.text.x =  element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=14)
  )

outFile <- file.path(outFolder, paste0("correlation_expression_purity_densityplot.", plotType))
ggsave(filename = outFile, plot = expr_p_density, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))    



##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
