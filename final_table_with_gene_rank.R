# add the rank of the genes in () after gene symbol
 
# Rscript final_table_with_gene_rank.R

require(foreach)
require(doMC)
registerDoMC(40)

outFolder <- "FINAL_TABLE_WITH_GENE_RANK"
dir.create(outFolder, recursive = TRUE)

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)


final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

gene_rank_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))


final_dt_with_rank <- foreach(i = 1:nrow(final_dt), .combine='rbind') %dopar% {
  
  curr_hicds <- final_dt$hicds[i]
  curr_exprds <- final_dt$expr[i]
  curr_region <- final_dt$region[i]
    
  all_geneSymbols <- unlist(strsplit(final_dt$region_genes[i], split = ","))
  stopifnot(all_geneSymbols %in% names(symb2entrez))
  all_entrezID <- symb2entrez[paste0(all_geneSymbols)]
  
  stopifnot(all_entrezID %in% gene_rank_dt$entrezID)
  
  all_geneRanks <- sapply(all_entrezID, function(entrez)  {
      entrez_rank <- gene_rank_dt$gene_rank[gene_rank_dt$hicds == curr_hicds &
                                              gene_rank_dt$exprds == curr_exprds &
                                                gene_rank_dt$entrezID == entrez]
      stopifnot(!is.na(entrez_rank))
      stopifnot(length(entrez_rank) == 1)
      entrez_rank
  })
  stopifnot(names(all_geneRanks) == all_geneSymbols)
  region_genes_with_rank <- paste0(paste0(all_geneSymbols, "(", all_geneRanks, ")"), collapse=",")
  
  out_dt <- final_dt[i,,drop=F]
  out_dt$region_genes <- region_genes_with_rank 
  out_dt
  
}

outFile <- file.path(outFolder, "final_dt_with_rank.Rdata")
save(final_dt_with_rank, file = outFile, version=2) 
cat(paste0("... written: ", outFile, "\n"))
