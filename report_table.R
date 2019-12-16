geneSignifThresh <-  0.05
tadSignifThresh <-  0.01

geneDT <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
tadDT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

tad_dt <- tadDT[,c("hicds", "exprds", "region", "adjPvalComb")]
tad_dt <- unique(tad_dt)

genes_dt <- geneDT[,c("hicds", "exprds", "entrezID", "adj.P.Val")]
genes_dt <- unique(genes_dt)

tad2_dt <- geneDT[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tad2_dt <- unique(tad2_dt)

stopifnot(nrow(tad_dt) == nrow(tad2_dt))


nGenesByDS <- aggregate(entrezID ~hicds + exprds, data=genes_dt, FUN=length)
nSignifGenesByDS <- aggregate(adj.P.Val ~hicds + exprds, data=genes_dt, FUN=function(x) sum(x <= geneSignifThresh))
mean(nGenesByDS$entrezID)
# [1] 8458.31
mean(nSignifGenesByDS$adj.P.Val)
# [1] 3892.741

nTADsByDS <- aggregate(region ~hicds + exprds, data=tad_dt, FUN=length)
nSignifTADsByDS <- aggregate(adjPvalComb ~hicds + exprds, data=tad_dt, FUN=function(x) sum(x <= tadSignifThresh))
mean(nTADsByDS$region)
# [1] 1739.379
mean(nSignifTADsByDS$adjPvalComb)
# [1] 25.05172