
startTime <- Sys.time()
cat(paste0("> Rscript cmp_corr_purity_signif_full_partial.R\n"))

script_name <- "cmp_corr_purity_signif_full_partial.R"


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


# Rscript cmp_corr_purity_signif_full_partial.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- TRUE

plotType <- "png"
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "CMP_CORR_PURITY_SIGNIF_FULL_PARTIAL"
dir.create(outFolder, recursive = TRUE)

full_signif_dt_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(full_signif_dt_file))
full_signif_dt <- get(load(full_signif_dt_file))
full_signif_dt$dataset <- file.path(full_signif_dt$hicds, full_signif_dt$exprds)

partial_signif_dt_file <- file.path("CREATE_FINAL_TABLE_PARTIAL/all_result_dt.Rdata")
stopifnot(file.exists(full_signif_dt_file))
partial_signif_dt <- get(load(partial_signif_dt_file))
partial_signif_dt$dataset <- file.path(partial_signif_dt$hicds, partial_signif_dt$exprds)

stopifnot(setequal(partial_signif_dt$dataset, full_signif_dt$dataset))

all_ds <- unique(partial_signif_dt$dataset)

pvalSignifThresh <- 0.01

full_onlysignif_dt <- full_signif_dt[full_signif_dt$adjPvalComb <= pvalSignifThresh,]
partial_onlysignif_dt <- partial_signif_dt[partial_signif_dt$adjPvalComb <= pvalSignifThresh,]

corr_expr_purity_dt <- get(load("CORR_EXPR_AND_PURITY_EPIC/corr_expr_purity_dt.Rdata"))

ds=all_ds[1]
ds=all_ds[39]
ds=all_ds[8]

all_purityCorr <- foreach(ds = all_ds) %dopar% {
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
  
  g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  g2t_dt <- read.delim(g2t_file, header=FALSE, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
  g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
  
  geneList <- get(load(file.path("PIPELINE/OUTPUT_FOLDER/", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
  stopifnot(geneList %in% g2t_dt$entrezID)
  g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
  
  full_signif_tads <- full_onlysignif_dt$region[full_onlysignif_dt$hicds==hicds & full_onlysignif_dt$exprds==exprds ]
  partial_signif_tads <- partial_onlysignif_dt$region[partial_onlysignif_dt$hicds==hicds & partial_onlysignif_dt$exprds==exprds ]
  
  
  full_signif_genes <- g2t_dt$entrezID[g2t_dt$region %in% full_signif_tads]
  stopifnot(full_signif_genes %in% geneList)
  full_signif_genes_fpkm <- names(pipeline_geneList)[pipeline_geneList %in% full_signif_genes]  
  
  partial_signif_genes <- g2t_dt$entrezID[g2t_dt$region %in% partial_signif_tads]
  stopifnot(partial_signif_genes %in% geneList)  
  partial_signif_genes_fpkm <- names(pipeline_geneList)[pipeline_geneList %in% partial_signif_genes]  
  
  # stopifnot(partial_signif_genes %in% corr_expr_purity_dt$entrezID[corr_expr_purity_dt$hicds == hicds & corr_expr_purity_dt$exprds == exprds])
  stopifnot(partial_signif_genes_fpkm %in% corr_expr_purity_dt$entrezID[corr_expr_purity_dt$hicds == hicds & corr_expr_purity_dt$exprds == exprds])
  # stopifnot(full_signif_genes %in% corr_expr_purity_dt$entrezID[corr_expr_purity_dt$hicds == hicds & corr_expr_purity_dt$exprds == exprds])
  stopifnot(full_signif_genes_fpkm %in% corr_expr_purity_dt$entrezID[corr_expr_purity_dt$hicds == hicds & corr_expr_purity_dt$exprds == exprds])
  
  full_signif_purityCorr <- corr_expr_purity_dt$all_corr_gene_purity[corr_expr_purity_dt$hicds == hicds &
                                                  corr_expr_purity_dt$exprds == exprds &
                                                  corr_expr_purity_dt$entrezID %in% full_signif_genes_fpkm
                                                  ]
  
  partial_signif_purityCorr <- corr_expr_purity_dt$all_corr_gene_purity[corr_expr_purity_dt$hicds == hicds &
                                                                       corr_expr_purity_dt$exprds == exprds &
                                                                       corr_expr_purity_dt$entrezID %in% partial_signif_genes_fpkm
                                                                     ]
  
  
  list(partial_signif_purityCorr=partial_signif_purityCorr,
       full_signif_purityCorr=full_signif_purityCorr)
  
}

names(all_purityCorr) <- all_ds
save(all_purityCorr, file=file.path(outFolder, "all_purityCorr.Rdata"), version=2)

all_partial_values <- unlist(lapply(all_purityCorr, function(x) x[["partial_signif_purityCorr"]]))
all_full_values <- unlist(lapply(all_purityCorr, function(x) x[["full_signif_purityCorr"]]))

plot_dt <- data.frame(
  corrWithPurity = c(all_full_values, all_partial_values),
  corrType = c(rep("full", length(all_full_values)), rep("partial", length(all_partial_values))),
  stringsAsFactors = FALSE
)

plotTit <- "Corr. expression and purity"
subTit <- "(genes from signif. TADs only)"

corrWithPur_signif_density <- ggdensity(plot_dt,
                            title = paste0(plotTit),
                            subtitle = paste0(subTit),
                            x = "corrWithPurity", 
                            col = "corrType",
                            palette="jco") + 
  labs(color="")+
  theme(
    plot.title = element_text(size=16, hjust=0.5, face = "bold"),
    plot.subtitle = element_text(size=14, hjust=0.5),
    strip.text.x =  element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=14)
  )

outFile <- file.path(outFolder, paste0("corrExprPurity_genes_from_signifTADs_full_partial_densityplot.", plotType))
ggsave(filename = outFile, plot = corrWithPur_signif_density, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))    