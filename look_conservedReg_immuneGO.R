goSignifThresh <- 0.01

tadSignifThresh <- 0.01
geneSignifThresh <- 0.01

inDT <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
inDT$entrezID <- as.character(inDT$entrezID)
inDT$dataset <- file.path(inDT$hicds, inDT$exprds)
tadSignif_inDT <-   inDT[inDT$tad_adjCombPval <= tadSignifThresh,]

go_result_dt <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
nrow(go_result_dt)
go_signif_dt <- go_result_dt[go_result_dt$p.adjust <= goSignifThresh]
nrow(go_signif_dt)

all_gos <- as.character(go_signif_dt$ID)

curr_go = all_gos[1]
foreach(curr_go = all_gos) %dopar% {
  
  curr_entrezID <- as.character(go_signif_dt$geneID[as.character(go_signif_dt$ID) == curr_go])
  stopifnot(length(curr_entrezID) == 1)
  go_entrezID <- unlist(strsplit(x=curr_entrezID, split="/"))
  stopifnot(all(go_entrezID %in% tadSignif_inDT$entrezID))
  
  go_inDT <- tadSignif_inDT[tadSignif_inDT$entrezID %in% go_entrezID,]
  
  occByGenes_dt <- aggregate(dataset ~ entrezID, data = go_inDT, FUN=length)
  
  nSignifByGenes_dt <- aggregate(adj.P.Val ~ entrezID, data = go_inDT, FUN=function(x) sum(x<=geneSignifThresh))
  
  
  occ_signif <- merge(occByGenes_dt, nSignifByGenes_dt, by="entrezID", all=T)
  stopifnot(!is.na(occ_signif))
  occ_signif$signifRatio <- occ_signif$adj.P.Val/occ_signif$dataset
  stopifnot(occ_signif$signifRatio>=0 & occ_signif$signifRatio <= 1)
  
  meanOccByGenes <- mean(occByGenes_dt$dataset)
  meanSignifByGenes <- mean(nSignifByGenes_dt$adj.P.Val)
  meanRatioSiginf <- mean(occ_signif$signifRatio)
  meanFC <- mean(go_inDT$logFC)
   
  
  
}