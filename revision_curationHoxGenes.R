outFolder <- file.path("REVISION_CURATIONHOXGENES")
dir.create(outFolder, recursive = TRUE)

# Rscript revision_curationHoxGenes.R

tadSignifThresh <- 0.01

result_dt <- get(load(file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")))

hox_dt <- result_dt[grepl("HOX", result_dt$region_genes),]


out_dt1 <- hox_dt[,c("hicds","exprds", "region", "region_genes", "meanLogFC","meanCorr" , "adjPvalComb")]
out_dt1 <- out_dt1[order(out_dt1$adjPvalComb, decreasing = F),]
outfile1 <- file.path(outFolder, "all_hoxTADs_all_ds.txt")
write.table(out_dt1, file=outfile1, sep="\t", quote=F, col.names=T, row.names = F)
cat(paste0("... written: ", outfile1, "\n"))

hox_dt$signif_up <- hox_dt$adjPvalComb <= tadSignifThresh & hox_dt$meanLogFC > 0

agg_nHoxTADsByDS_dt <- aggregate(region~hicds+exprds, data=hox_dt, FUN=length)
colnames(agg_nHoxTADsByDS_dt)[colnames(agg_nHoxTADsByDS_dt) == "region"] <- "nHoxTADs"

agg_nUpHoxTADsByDS_dt <- aggregate(meanLogFC~hicds+exprds, data=hox_dt, FUN=function(x)sum(x>0))
colnames(agg_nUpHoxTADsByDS_dt)[colnames(agg_nUpHoxTADsByDS_dt) == "meanLogFC"] <- "nUpHoxTADs"

agg_nSignifHoxTADsByDS_dt <- aggregate(adjPvalComb~hicds+exprds, data=hox_dt, FUN=function(x) sum(x <= tadSignifThresh))
colnames(agg_nSignifHoxTADsByDS_dt)[colnames(agg_nSignifHoxTADsByDS_dt) == "adjPvalComb"] <- "nSignifHoxTADs"

agg_nSignifUpHoxTADsByDS_dt <- aggregate(signif_up~hicds+exprds, data=hox_dt, FUN=sum)
colnames(agg_nSignifUpHoxTADsByDS_dt)[colnames(agg_nSignifUpHoxTADsByDS_dt) == "signif_up"] <- "nSignifUpHoxTADs"


merged_dt <- merge(agg_nSignifUpHoxTADsByDS_dt, 
      merge(agg_nSignifHoxTADsByDS_dt, 
            merge(agg_nHoxTADsByDS_dt, agg_nUpHoxTADsByDS_dt, by=c("hicds", "exprds"), all=T), by=c("hicds", "exprds"), all=T), by=c("hicds", "exprds"), all=T)

stopifnot(!is.na(merged_dt))

merged_dt$ratioSignifTADs <- merged_dt$nSignifHoxTADs/merged_dt$nHoxTADs
stopifnot(merged_dt$ratioSignifTADs >= 0 & merged_dt$ratioSignifTADs <= 1)

out_dt2 <- merged_dt[,c("hicds","exprds","nHoxTADs", "nSignifHoxTADs" , "ratioSignifTADs", "nUpHoxTADs","nSignifUpHoxTADs")]
out_dt2 <- out_dt2[order(out_dt2$ratioSignifTADs, decreasing = TRUE),]

outfile2 <- file.path(outFolder, "sum_and_ratio_hoxTADs_all_ds.txt")
write.table(out_dt2, file=outfile2, sep="\t", quote=F, col.names=T, row.names = F)
cat(paste0("... written: ", outfile2, "\n"))
##### take gene rank tad rank -> for hox genes count how many dataset +  how many times signif + how many times signif up

gene_dt <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

gene_dt$symb <- entrez2symb[paste0(gene_dt$entrezID)]
stopifnot(!is.na(gene_dt$symb))

hoxgenes_dt <- gene_dt[grepl("^HOX", gene_dt$symb),]
hoxgenes_dt$signif_up <- hoxgenes_dt$tad_adjCombPval <= tadSignifThresh & hoxgenes_dt$logFC > 0


nDS_dt <- aggregate(hicds ~ symb, FUN=length, data=hoxgenes_dt)
colnames(nDS_dt)[colnames(nDS_dt) == "hicds"] <- "nDS"

nSignifTAD_dt <- aggregate(tad_adjCombPval ~ symb, data=hoxgenes_dt, FUN=function(x) sum(x<=tadSignifThresh))
colnames(nSignifTAD_dt)[colnames(nSignifTAD_dt) == "tad_adjCombPval"] <- "nInSignifTAD"

nSignifUpTAD_dt <- aggregate(signif_up ~ symb, FUN=sum, data=hoxgenes_dt)
colnames(nSignifUpTAD_dt)[colnames(nSignifUpTAD_dt) == "signif_up"] <- "nUpInSignifTAD"


gene_merged_dt <-  merge(nSignifTAD_dt, 
                         merge(nSignifUpTAD_dt, nDS_dt, by=c("symb"), all=T), by=c("symb"), all=T)
stopifnot(!is.na(gene_merged_dt))
gene_merged_dt$ratioSignif <- gene_merged_dt$nInSignifTAD/gene_merged_dt$nDS
gene_merged_dt$ratioUpSignif <- gene_merged_dt$nUpInSignifTAD/gene_merged_dt$nDS

stopifnot(gene_merged_dt$ratioSignif >= 0 & gene_merged_dt$ratioSignif <= 1)
stopifnot(gene_merged_dt$ratioUpSignif >= 0 & gene_merged_dt$ratioUpSignif <= 1)
stopifnot(gene_merged_dt$nDS <= 58)

out_dt3 <- gene_merged_dt[,c("symb","nDS"  ,"nInSignifTAD","nUpInSignifTAD", "ratioSignif", "ratioUpSignif" )]
out_dt3 <- out_dt3[order(out_dt3$ratioSignif, decreasing = TRUE),]

             
outfile3 <- file.path(outFolder, "sum_and_ratio_hoxGenes_all_ds.txt")
write.table(out_dt3, file=outfile3, sep="\t", quote=F, col.names=T, row.names = F)
cat(paste0("... written: ", outfile3, "\n"))

