########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript coexpr_expr_and_purity.R\n"))

script_name <- "cexpr_expr_and_purity.R"

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
  fullcorr_dt$gene1 <- as.character(pmin(fullcorr_dt$Var1, fullcorr_dt$Var2))
  fullcorr_dt$gene2 <- as.character(pmax(fullcorr_dt$Var1, fullcorr_dt$Var2))
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
  partialcorr_dt$gene1 <- as.character(pmin(partialcorr_dt$Var1, partialcorr_dt$Var2))
  partialcorr_dt$gene2 <- as.character(pmax(partialcorr_dt$Var1, partialcorr_dt$Var2))
  colnames(partialcorr_dt)[colnames(partialcorr_dt) == "value"] <- newcol
  partialcorr_dt <- partialcorr_dt[,c("gene1", "gene2", newcol)]
  return(partialcorr_dt)
}

# Rscript coexpr_and_purity.R

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

outFolder <- file.path("COEXPR_AND_PURITY")
dir.create(outFolder, recursive = TRUE)

purity_file <- file.path("tcga_purity_aran2015.csv")
purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
pm <- purity_metrics[1]
# all the ranks are between 1 and 0

all_hicds <- list.files(pipFolder)
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]
# all_ds=all_ds[1]

corMet <- "pearson"

cat(paste0("!!! HARD-CODED !!!\n"))
cat(paste0(">>> corMet\t=\t", corMet, "\n"))
cat(paste0(">>> purity metric\t=\t", pm, "\n"))

if(buildTable) {
  
  coexpr_and_purity_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
    
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
        gene1=NA,
        gene2=NA,
        coexpr=NA,
        coexpr_samp1=NA,
        coexpr_samp2=NA,
        partial_coexpr=NA,
        partial_coexpr_samp1=NA,
        partial_coexpr_samp2=NA,
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
      
   out_dt
  } # end-iterating ds 
    
  
  outFile <- file.path(outFolder, "coexpr_and_purity_dt.Rdata")
  save(coexpr_and_purity_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {
  outFile <- file.path(outFolder, "coexpr_and_purity_dt.Rdata")
  coexpr_and_purity_dt <- get(load(outFile))
  
}

# plot_dt <- melt(corr_expr_purity_dt, id=c("hicds", "exprds", "entrezID"))
# 
# plotTit <- paste0("Corr. expression and purity")
# subTit <- paste0("corr. method = ", corMet)
# 
# myxlab <- ""
# myylab <- paste0("purity (", pm, ")")
# 
# expr_p_boxplot <- ggboxplot(plot_dt, x = "variable",
#                             y = "value",
#                             xlab = myxlab,
#                             # scales='free_y',
#                             # y = gene_order,
#                             # combine = TRUE,
#                             ylab = myylab,
#                             palette = "jco") + 
#   guides(color=FALSE)+
#   # ggtitle(plotTit, sub=plotSubTit)+
#   ggtitle(plotTit)+
#   theme(
#     plot.title = element_text(size=16, hjust=0.5, face = "bold"),
#     plot.subtitle = element_text(size=14, hjust=0.5),
#     strip.text.x =  element_text(size=12),
#     axis.text.x = element_text(size=12),
#     axis.title.y = element_text(size=14)
#   )
# 
# 
# outFile <- file.path(outFolder, paste0("correlation_expression_purity_boxplot.", plotType))
# ggsave(filename = outFile, plot = expr_p_boxplot, height=myHeightGG, width=myWidthGG)
# cat(paste0("... written: ", outFile, "\n"))    
# 
# 
# expr_p_density <- ggdensity(plot_dt,
#                             title = paste0(plotTit),
#                             subtitle = paste0(subTit),
#                             x = "value", 
#                             col = "variable",
#                             palette="jco") + 
#   labs(color="")+
#   theme(
#     plot.title = element_text(size=16, hjust=0.5, face = "bold"),
#     plot.subtitle = element_text(size=14, hjust=0.5),
#     strip.text.x =  element_text(size=12),
#     axis.text.x = element_text(size=12),
#     axis.title.y = element_text(size=14)
#   )
# 
# outFile <- file.path(outFolder, paste0("correlation_expression_purity_densityplot.", plotType))
# ggsave(filename = outFile, plot = expr_p_density, height=myHeightGG, width=myWidthGG)
# cat(paste0("... written: ", outFile, "\n"))    

# 
# x=get(load("COEXPR_AND_PURITY/coexpr_and_purity_dt.Rdata"))
# y=get(load("COEXPR_AND_PURITY_FOREACH//coexpr_and_purity_dt.Rdata"))
# 
# nrow(x)
# nrow(y)
# 
# uu <- merge(x, y, by=c("hicds", "exprds","gene1","gene2"))
# 
# all(zapsmall(uu[,c("coexpr.y")]) == zapsmall(uu[,c("coexpr.x")] ))
# all(zapsmall(uu[,c("partial_coexpr.y")]) == zapsmall(uu[,c("partial_coexpr.x")] ))
# all(zapsmall(uu[,c("partial_coexpr_samp1.y")]) == zapsmall(uu[,c("partial_coexpr_samp1.x")] ))
# all(zapsmall(uu[,c("partial_coexpr_samp2.y")]) == zapsmall(uu[,c("partial_coexpr_samp2.x")] ))
# all(zapsmall(uu[,c("partial_coexpr_samp1.y")]) == zapsmall(uu[,c("partial_coexpr_samp1.x")] ))
# all(zapsmall(uu[,c("partial_coexpr_samp2.y")]) == zapsmall(uu[,c("partial_coexpr_samp2.x")] ))
# 
# uu[,c("gene1", "gene2", "coexpr.x","coexpr.y")]
# 
# # gene1  gene2     coexpr.x     coexpr.y
# # 1   10357 155060 -0.006653467 -0.006347639
# # 2   10357 388795 -0.007596405 -0.009245269
# # 3   10357 390284  0.021291232  0.021014074
# # 
# # all_dt[all_dt$gene1=="10357" & all_dt$gene2 == "155060",]
# 
# uu[,c("gene1", "gene2", "partial_coexpr.x","partial_coexpr.y")]


##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
