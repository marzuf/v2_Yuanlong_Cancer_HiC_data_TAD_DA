setDir <- "/media/electron"
setDir <- ""

# Rscript revision_expressionData_all_matchedCL.R
script_name="revision_expressionData_all_matchedCL.R"
startTime <- Sys.time()
cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
minGenes <- 3

log10_offset <- 0.01

outFolder <- "REVISION_EXPRESSIONDATA_ALL_MATCHEDCL"
dir.create(outFolder, recursive = TRUE)

#### will not look at mean corr because max number of replicates is 6, with median=2

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


expr_var <- "FPKM"
symb_aggFun <- "mean"  # aggregate expression data to have unique gene symbol
rep_aggFun <- "mean"  # aggregate expression across replicates

tad_aggFun <- "mean" # aggregate expression to TAD level

# all_mrna_data=get(load("sub_expression_allDS.Rdata"))

all_mrna_data <-get(load(file.path(setDir,
                                   "/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/1.Data/exp_CELL_LINEs/exp_info.Rdata")))

source("revision_settings.R")

# all_hicds=all_hicds[1]

hicds=c("Barutcu_MCF-10A_40kb")

hicds="ENCSR489OCU_NCI-H460_40kb"

all_matched_mRNA_dt <- foreach(hicds = names(all_hicds), .combine='rbind') %dopar% {
  
  cat(paste0("> START : ", hicds, "\n"))
  
  stopifnot(hicds %in% names(all_hicds))
  
  ds <- as.character(all_hicds[paste0(hicds)])
  stopifnot(length(ds) > 0)
  stopifnot(!is.na(ds))
  stopifnot(ds %in% names(all_mrna_data))
  
  ds_mrna_data <- all_mrna_data[[ds]]
  
  
  ## laod the gene2tad assignment
  region_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  genes2tad_dt <- read.delim(region_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID",  "chromo", "start", "end", "region"))
  genes2tad_dt$entrezID <- as.character(genes2tad_dt$entrezID)                                                               
  stopifnot(genes2tad_dt$entrezID %in% names(entrez2symb))
  genes2tad_dt$symbol <- entrez2symb[paste0(genes2tad_dt$entrezID)]
  stopifnot(!is.na(genes2tad_dt$symbol))
  # take only TAD regions
  genes2tad_dt <- genes2tad_dt[grepl("_TAD", genes2tad_dt$region),]
  
  expr_values <- lapply(ds_mrna_data, function(x) {
    stopifnot("external_gene_name" %in% colnames(x))
    stopifnot(paste0(expr_var) %in% colnames(x))
    # stopifnot(!duplicated(x[["external_gene_name"]])) ### NOT TRUE
    agg_dt <- aggregate(as.formula(paste0(expr_var, "~external_gene_name")), data=x, FUN=symb_aggFun, na.rm=TRUE)
    stopifnot(!duplicated(agg_dt[["external_gene_name"]])) # now should be
    setNames(agg_dt[[paste0(expr_var)]], agg_dt[["external_gene_name"]])
  })
  # first check that they have the same gene list
  stopifnot(length(unique(lapply(expr_values, names))) == 1)
  # as they are the same, I can sum
  expr_values_agg <- apply( do.call(rbind, expr_values), 2,rep_aggFun)
  stopifnot(names(expr_values_agg) == names(expr_values[[1]]))
  # 5S_rRNA        7SK       A1BG   A1BG-AS1       A1CF        A2M
  # 0.0000000  0.1508333 34.3700000  5.3150000  0.0000000  0.0000000
  # some checks
  tmp <- sort(expr_values_agg, decreasing = TRUE)
  stopifnot(tmp[1] == do.call(rep_aggFun, list(unlist(lapply(expr_values, function(x)x[[names(tmp)[1]]])))))
  stopifnot(tmp[3] == do.call(rep_aggFun, list(unlist(lapply(expr_values, function(x)x[[names(tmp)[3]]])))))  
  
  stopifnot(any(genes2tad_dt$symbol %in% names(expr_values_agg)))
  
  
  genes2tad_dt$gene_expr <- expr_values_agg[paste0(genes2tad_dt$symbol)]
  genes2tad_dt$gene_expr_log10 <- log10(genes2tad_dt$gene_expr+log10_offset)
  
  genes2tad_dt <- genes2tad_dt[!is.na(genes2tad_dt$gene_expr),]
  stopifnot(!is.na(genes2tad_dt))
  
  nGenes_withExpr <- setNames(as.numeric(table(genes2tad_dt$region)), names(table(genes2tad_dt$region)))

  log_dt <- aggregate(gene_expr_log10 ~ region, data=genes2tad_dt, FUN = tad_aggFun)
  colnames(log_dt)[colnames(log_dt) == "gene_expr_log10"] <- paste0(tad_aggFun, "_geneExprLog10")
  stopifnot(log_dt$region %in% names(nGenes_withExpr))
  stopifnot(!duplicated(log_dt$region))
  log_aggValues <- setNames(log_dt[,  paste0(tad_aggFun, "_geneExprLog10")], log_dt$region)
    
    
  agg_dt <- aggregate(gene_expr ~ region, data=genes2tad_dt, FUN = tad_aggFun)
  colnames(agg_dt)[colnames(agg_dt) == "gene_expr"] <- paste0(tad_aggFun, "_geneExpr")
  stopifnot(agg_dt$region %in% names(nGenes_withExpr))
  stopifnot(agg_dt$region %in% names(log_aggValues))
  
  agg_dt[,paste0(tad_aggFun, "_geneExprLog10")] <- log_aggValues[paste0(agg_dt$region)]
  
  agg_dt$nGenes_withExpr <- nGenes_withExpr[paste0(agg_dt$region)]
  agg_dt$hicds <- hicds
  save(agg_dt, file="agg_dt.Rdata", version=2)
  agg_dt[,c("hicds", "region", "nGenes_withExpr", paste0(tad_aggFun, "_geneExpr"), paste0(tad_aggFun, "_geneExprLog10"))]
  
  
}
stopifnot(!duplicated(file.path(all_matched_mRNA_dt$hicds, all_matched_mRNA_dt$region)))
outFile <- file.path(outFolder,paste0("all_matched_mRNA_dt.Rdata"))
save(all_matched_mRNA_dt, file=outFile, version=2)
cat(paste0("written: ", outFile, "\n"))


  

######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

# tmp <- list(c("A"=1,"B"=2,"C"=3), c("A"=3,"B"=6,"C"=9), c("A"=2, "B"=1, "C" =3))
# Reduce(rep_aggFun, tmp)  ### !!! not working
# Reduce(`+`, tmp)/length(tmp)
# colMeans(do.call(rbind, tmp))
# all(colMeans(do.call(rbind, tmp)) == Reduce(`+`, tmp)/length(tmp))
#
# apply( do.call(rbind, tmp), 2,rep_aggFun)
# all(apply( do.call(rbind, tmp), 2,rep_aggFun) == Reduce(`+`, tmp)/length(tmp))

# work only if length >= 2
# all(sapply(lapply(expr_values[2:length(expr_values)], names), FUN = identical, names(expr_values[[1]])))
# all(sapply(lapply(tmp[2:length(tmp)], names), FUN = identical, names(tmp[[1]])))
# > tmp <- list(c("A"=1,"B"=2,"C"=3), c("A"=3,"B"=6,"C"=9))
# > unique(lapply(tmp, names))
# [[1]]
# [1] "A" "B" "C"
#
# > Reduce(`+`, tmp)
# A  B  C
# 4  8 12

# > tmp <- list(c(1,2,3), c(1,2,3))
# > unique(tmp)
# [[1]]
# [1] 1 2 3
#
# > tmp <- list(c(1,2,3), c(1,2,4))
# > unique(tmp)
# [[1]]
# [1] 1 2 3
#
# [[2]]
# [1] 1 2 4
#
# >

# > sub_data2=x[["mega_GSE105318_DLD1"]]
# > str(sub_data2)
# List of 2
# $ rep1:Classes ‘data.table’ and 'data.frame':  54849 obs. of  9 variables:
#   ..$ ensembl_gene_id   : chr [1:54849] "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" "ENSG00000000938" ...
# ..$ TPM               : num [1:54849] 79.87 2.28 17.89 0 0.17 ...
# ..$ FPKM              : num [1:54849] 70.91 2.03 15.88 0 0.15 ...
# ..$ expected_count    : num [1:54849] 499 81.7 256.3 0 2 ...
# ..$ external_gene_name: chr [1:54849] "DPM1" "SCYL3" "C1orf112" "FGR" ...
# ..$ chromosome_name   : num [1:54849] 20 1 1 1 1 6 6 6 1 1 ...
# ..$ start_position    : int [1:54849] 49551404 169818772 169631245 27938575 196621008 143815948 53362139 41040684 24683489 24742284 ...
# ..$ end_position      : int [1:54849] 49575092 169863408 169823221 27961788 196716634 143832827 53481768 41067715 24743424 24799466 ...
# ..$ gene_biotype      : chr [1:54849] "protein_coding" "protein_coding" "protein_coding" "protein_coding" ...
# ..- attr(*, ".internal.selfref")=<externalptr>
#   ..- attr(*, "sorted")= chr "ensembl_gene_id"
# $ rep2:Classes ‘data.table’ and 'data.frame':  54849 obs. of  9 variables:
#   ..$ ensembl_gene_id   : chr [1:54849] "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" "ENSG00000000938" ...
# ..$ TPM               : num [1:54849] 88.64 2.61 16.89 0 0.2 ...
# ..$ FPKM              : num [1:54849] 75.12 2.21 14.32 0 0.17 ...
# ..$ expected_count    : num [1:54849] 627 105 279 0 2 ...
# ..$ external_gene_name: chr [1:54849] "DPM1" "SCYL3" "C1orf112" "FGR" ...
# ..$ chromosome_name   : num [1:54849] 20 1 1 1 1 6 6 6 1 1 ...
# ..$ start_position    : int [1:54849] 49551404 169818772 169631245 27938575 196621008 143815948 53362139 41040684 24683489 24742284 ...
# ..$ end_position      : int [1:54849] 49575092 169863408 169823221 27961788 196716634 143832827 53481768 41067715 24743424 24799466 ...
# ..$ gene_biotype      : chr [1:54849] "protein_coding" "protein_coding" "protein_coding" "protein_coding" ...
# ..- attr(*, ".internal.selfref")=<externalptr>
#   

# -> for each TAD
# -> how many genes available in the  rna expr
# -> mean FPKM
# -> then I can do the same as for norm and tumor
# Danielle prefer TPM, but I use FPKM

