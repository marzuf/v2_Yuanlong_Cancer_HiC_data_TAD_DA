options(scipen=100)

setDir <- ""

# Rscript cmp_same_conditions.R 

script_name <- "cmp_same_conditions.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(reshape2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4
labCex <- 1.4

tiesMeth <- "min"


pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CMP_SAME_CONDITIONS"
dir.create(outFolder, recursive=TRUE)


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2 | length(args) == 0)
hicds <- args[1]
exprds <- args[2]

if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}


# retrieve the "duplicated" all_exprds
all_exprds_vect <- unlist(all_exprds)
dup_exprds <- unique(all_exprds_vect[duplicated(all_exprds_vect)])

expr_hicds <- sapply(dup_exprds, function(exprds) {
  names(all_exprds)[sapply(all_exprds, function(x) exprds %in% x )]
})

stopifnot(dup_exprds %in% names(expr_hicds))
  
### BUILD SIGNIF ALONG FDR THRESH
foo <- foreach(exprds = dup_exprds, .combine='rbind') %dopar% {
  
  all_hicds <- expr_hicds[[paste0(exprds)]]
  
  all_geneList <- lapply(all_hicds, function(hicds) {
    geneList_file <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneList_file))
    eval(parse(text=load(geneList_file)))
  })
  names(all_geneList) <- all_hicds
  commonGenes <- Reduce(intersect, all_geneList)
  
  stopifnot(!duplicated(commonGenes))
  
  foo <- lapply(all_geneList, function(x) stopifnot(commonGenes %in% x))
  
  # retrieve the gene rank from DE and gene TAD rank from pval comb
  all_gene_DE_Rank <- lapply(all_hicds, function(hicds) {
    de_file <- file.path(pipOutFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")
    stopifnot(file.exists(de_file))
    de_DT <- eval(parse(text = load(de_file)))
    de_DT$genes <- as.character(de_DT$genes)
    g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    geneList_file <- file.path(pipOutFolder, hicds, exprds,script0_name,   "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneList_file))
    geneList <- eval(parse(text = load(geneList_file)))
    stopifnot(names(geneList) %in% de_DT$genes)
    stopifnot(geneList %in% g2t_DT$entrezID)
    de_DT <- de_DT[de_DT$genes %in% names(geneList),]
    stopifnot(nrow(de_DT) == length(geneList))
    de_DT$entrezID <- sapply(de_DT$genes, function(x) as.character(geneList[as.character(x)]))
    stopifnot(!is.na(de_DT$entrezID))
    gene_DE <- setNames(de_DT$adj.P.Val, de_DT$entrezID)
    gene_DE_rank <- rank(gene_DE, ties = tiesMeth) # best rank = smallest value
    stopifnot(gene_DE[names(gene_DE_rank)[gene_DE_rank == 1]] == min(gene_DE))
    gene_rank_DT <- data.frame(entrezID=names(gene_DE_rank), gene_rank = gene_DE_rank, stringsAsFactors = FALSE)
    stopifnot(commonGenes %in% gene_rank_DT$entrezID)
    gene_rank_DT <- gene_rank_DT[gene_rank_DT$entrezID %in% commonGenes,]
    gene_rank_DT <- gene_rank_DT[order(gene_rank_DT$entrezID),]
    stopifnot(nrow(gene_rank_DT) == length(commonGenes))
    gene_rank_DT
  })
  names(all_gene_DE_Rank) <- all_hicds
  
  # RETRIEVE PVAL COMBINED DATA
  
  all_gene_TAD_DE_Rank <- lapply(all_hicds, function(hicds) {
    pvalFile <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
    stopifnot(file.exists(pvalFile))
    combPvals <- eval(parse(text = load(pvalFile)))
    adj_combPvals <- p.adjust(combPvals, method="BH")
    adj_combPvals_rank <- rank(adj_combPvals, ties = tiesMeth) # best rank = smallest value
    stopifnot(adj_combPvals[names(adj_combPvals_rank)[adj_combPvals_rank == 1]] == min(adj_combPvals))
    tad_rank_DT <- data.frame(region=names(adj_combPvals_rank), TAD_rank = adj_combPvals_rank, stringsAsFactors = FALSE)
    geneList_file <- file.path(pipOutFolder, hicds, exprds,script0_name,   "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneList_file))
    geneList <- eval(parse(text = load(geneList_file)))
    g2tFile <- file.path(pipFolder, hicds,  "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    stopifnot(geneList %in% g2t_DT$entrezID)
    g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
    gene_rank_DT <- all_gene_DE_Rank[[paste0(hicds)]]
    g2t_rank_DT <- merge(gene_rank_DT, g2t_DT[,c("entrezID", "region")], by="entrezID", all=TRUE)
    stopifnot(setequal(tad_rank_DT$region, g2t_rank_DT$region))
    gene_tad_rank_DT <- merge(tad_rank_DT, g2t_rank_DT, by ="region")
    stopifnot(commonGenes %in% gene_tad_rank_DT$entrezID)
    gene_tad_rank_DT <- gene_tad_rank_DT[gene_tad_rank_DT$entrezID %in% commonGenes,]
    gene_tad_rank_DT <- gene_tad_rank_DT[order(gene_tad_rank_DT$entrezID),]
    stopifnot(nrow(gene_tad_rank_DT) == length(commonGenes))
    gene_tad_rank_DT
  })
  names(all_gene_TAD_DE_Rank) <- all_hicds
  
  hicds_comb <- combn(all_hicds,2)
  
  for(i in 1:ncol(hicds_comb)) {
    hicds1 <- hicds_comb[1,i]
    hicds2 <- hicds_comb[2,i]
    stopifnot(hicds1 %in% names(all_gene_TAD_DE_Rank))
    stopifnot(hicds2 %in% names(all_gene_TAD_DE_Rank))
    dt1 <- all_gene_TAD_DE_Rank[[paste0(hicds1)]]
    dt2 <- all_gene_TAD_DE_Rank[[paste0(hicds2)]]
    
    stopifnot(nrow(dt1) == length(commonGenes))
    stopifnot(nrow(dt2) == length(commonGenes))
    
    stopifnot(dt1$entrezID == dt2$entrezID)
    
    all_vars <- c("TAD_rank", "gene_rank")
    plot_var = all_vars[1]
    for(plot_var in all_vars) {
      var1 <- dt1[,paste0(plot_var)]
      var2 <- dt2[,paste0(plot_var)]
      outFile <- file.path(outFolder, paste0(exprds, "_", hicds1, "_vs_", hicds2, "_", plot_var, ".", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      densplot(
        x = var1,
        y = var2,
        cex.axis = axisCex,
        cex.lab = labCex,
        xlab = paste0(hicds1, "\n", plot_var), 
        ylab = paste0(hicds2, "\n", plot_var), 
        main = paste0(exprds),
        sub = paste0(hicds2, " vs. ", hicds1)
      )
      mtext(side=3, text = paste0("(n=", nrow(dt1),")"), font=3)
      curve(1*x, add=TRUE, lty=2, col="black", lwd=1.5)
      addCorr(
        x = var1, y = var2,
        legPos = "topleft", bty="n"
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    } # end-for iterating over plot_var (TAD_rank or gene_rank)
  } # end-for iterating over pairs of hic datasets
} # end-for iterating over exprds





##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



