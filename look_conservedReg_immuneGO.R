# Rscript look_conservedReg_immuneGO.R

require(foreach)
require(doMC)
require(EPIC)
registerDoMC(40)

goSignifThresh <- 0.01
tadSignifThresh <- 0.01
geneSignifThresh <- 0.01
tieMeth <- "min"
strwdth <- 25

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 9

# setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

Tref_sub <- TRef[["sigGenes"]][TRef[["sigGenes"]] %in% names(symb2entrez)]
Tref_entrez <- symb2entrez[paste0(Tref_sub)]

Bref_sub <- BRef[["sigGenes"]][BRef[["sigGenes"]] %in% names(symb2entrez)]
Bref_entrez <- symb2entrez[paste0(Bref_sub)]

outFolder <- file.path("LOOK_CONSERVEDREG_IMMUNEGO")
dir.create(outFolder, recursive = TRUE)

inDT <- get(load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
inDT$entrezID <- as.character(inDT$entrezID)
inDT$dataset <- file.path(inDT$hicds, inDT$exprds)
tadSignif_inDT <-   inDT[inDT$tad_adjCombPval <= tadSignifThresh,]

go_result_dt <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
nrow(go_result_dt)
go_signif_dt <- go_result_dt[go_result_dt$p.adjust <= goSignifThresh,]
nrow(go_signif_dt)
go_signif_dt$p.adjust_rank <- rank(go_signif_dt$p.adjust, tie=tieMeth)

all_gos <- as.character(go_signif_dt$ID)
curr_go = all_gos[1]

go_stat_dt <- foreach(curr_go = all_gos, .combine='rbind') %dopar% {
  go_rank <- unique(go_signif_dt$p.adjust_rank[as.character(go_signif_dt$ID) == curr_go])
  stopifnot(length(go_rank) == 1)
  curr_entrezID <- as.character(go_signif_dt$geneID[as.character(go_signif_dt$ID) == curr_go])
  stopifnot(length(curr_entrezID) == 1)
  go_entrezID <- unlist(strsplit(x=curr_entrezID, split="/"))
  stopifnot(all(go_entrezID %in% tadSignif_inDT$entrezID))
  go_symbol <- entrez2symb[paste0(go_entrezID)]
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
  median_tadRank <- median(go_inDT$tad_rank)
  median_geneRank <- median(go_inDT$gene_rank)
  data.frame(
    go_id = curr_go,
    go_rank=go_rank,
    median_tadRank=median_tadRank,
    median_geneRank=median_geneRank,
    meanFC=meanFC,
    meanOccByGenes=meanOccByGenes,
    meanSignifByGenes=meanSignifByGenes,
    meanRatioSiginf=meanRatioSiginf,
    go_genes=paste0(go_symbol, collapse=","),
    stringsAsFactors = FALSE
  )
}

outFile <- file.path(outFolder, "go_stat_dt.Rdata")
save(go_stat_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# load("LOOK_CONSERVEDREG_IMMUNEGO/go_stat_dt.Rdata")

go_stat_dt <- go_stat_dt[order(go_stat_dt$go_rank),]
plot_dt <- go_stat_dt
plot_dt$go_id_lab <- unlist(lapply(strwrap(gsub("_", " ", plot_dt$go_id), 
                      width = strwdth, simplify=FALSE), 
              function(x) paste0(x, collapse="\n")))
plot_dt$go_id <- factor(plot_dt$go_id, levels=plot_dt$go_id)
plot_dt$go_id_lab <- factor(plot_dt$go_id_lab, levels=plot_dt$go_id_lab)

subTit <- "agg. mean values"

allcols <- colnames(go_stat_dt)
allcols <- allcols[!allcols %in% c("go_id", "go_rank", "go_genes")]
mycol = allcols[1]
for(mycol in allcols) {
  plotTit <- paste0(mycol, " by ranked GOs")
  bar_p <- ggplot(plot_dt, aes_string(x="go_id_lab", y=paste0(mycol)))+
    ggtitle(plotTit, subtitle=subTit)+
    geom_bar(stat="identity")+
    labs(x="", y =mycol)+
    theme(
      plot.margin = margin(b = 2, unit = "cm"),
      text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_text(size=8, hjust=1, vjust=0.5, angle=90),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    ) 
  outFile <- file.path(outFolder, paste0("agg_", mycol, "_barplot.", plotType))
  ggsave(bar_p, filename=outFile, height = myHeightGG, width=mywidthGG)
  cat(paste0("... written: ", outFile, "\n"))
}
  
  
full_go_stat_dt <- foreach(curr_go = all_gos, .combine='rbind') %dopar% {
  curr_entrezID <- as.character(go_signif_dt$geneID[as.character(go_signif_dt$ID) == curr_go])
  stopifnot(length(curr_entrezID) == 1)
  go_entrezID <- unlist(strsplit(x=curr_entrezID, split="/"))
  stopifnot(all(go_entrezID %in% tadSignif_inDT$entrezID))
  
  go_rank <- unique(go_signif_dt$p.adjust_rank[as.character(go_signif_dt$ID) == curr_go])
  stopifnot(length(go_rank) == 1)
  
  go_symbol <- entrez2symb[paste0(go_entrezID)]
  
  go_inDT <- tadSignif_inDT[tadSignif_inDT$entrezID %in% go_entrezID,]
  
  occByGenes_dt <- aggregate(dataset ~ entrezID, data = go_inDT, FUN=length)
  
  nSignifByGenes_dt <- aggregate(adj.P.Val ~ entrezID, data = go_inDT, FUN=function(x) sum(x<=geneSignifThresh))
  
  go_inDT$go_id <- curr_go
  go_inDT$go_rank <- go_rank
  
  go_inDT
  
}

outFile <- file.path(outFolder, "full_go_stat_dt.Rdata")
save(full_go_stat_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# load("LOOK_CONSERVEDREG_IMMUNEGO/full_go_stat_dt.Rdata")

plot_dt2 <- full_go_stat_dt[,c("go_id", "go_rank", "tad_rank", "gene_rank", "logFC")]
require(reshape2)
m_plot_dt2 <- melt(plot_dt2, id=c("go_id", "go_rank"))
m_plot_dt2 <- m_plot_dt2[order(m_plot_dt2$go_rank),]

m_plot_dt2$go_id_lab <- unlist(lapply(strwrap(gsub("_", " ", m_plot_dt2$go_id), 
                                           width = strwdth, simplify=FALSE), 
                                   function(x) paste0(x, collapse="\n")))

m_plot_dt2$go_id <- factor(m_plot_dt2$go_id, levels = unique(m_plot_dt2$go_id))
m_plot_dt2$go_id_lab <- factor(m_plot_dt2$go_id_lab, levels=unique(m_plot_dt2$go_id_lab))


allvars <- unique(m_plot_dt2$variable)
myvar = allvars[1]
subTit <- "all values"
for(myvar in allvars) {
  sub_dt <- m_plot_dt2[m_plot_dt2$variable == myvar,]
  plotTit <- paste0(myvar, " by ranked GOs")
  box_p <- ggplot(sub_dt, aes(x=go_id_lab, y=value))+
    ggtitle(plotTit, subtitle=subTit)+
    geom_boxplot()+
    labs(x="", y =myvar)+
    theme(
      plot.margin = margin(b = 2, unit = "cm"),
      text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_text(size=8, hjust=1, vjust=0.5, angle=90),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    ) 
  outFile <- file.path(outFolder, paste0("all_", mycol, "_boxplot.", plotType))
  ggsave(box_p, filename=outFile, height = myHeightGG, width=mywidthGG)
  cat(paste0("... written: ", outFile, "\n"))
}



########################################################################################
########################################################################################
########################################################################################

stopifnot(!duplicated(go_signif_dt$ID))


all_go_genes_dt <- do.call(rbind, by(go_signif_dt, go_signif_dt$ID, function(x) {
  
  go <- unique(x$ID)
  stopifnot(length(go) == 1)
  go_rank <- unique(x$p.adjust_rank)
  stopifnot(length(go_rank) == 1)
  
  go_padj <-unique(x$p.adjust)
  stopifnot(length(go_padj) == 1)
  
  go_entrez <-unique(x$geneID)
  stopifnot(length(go_entrez) == 1)
  go_entrezID <- unlist(strsplit(x=go_entrez, split="/"))
  
  data.frame(
    GO_id = go,
    p.adjust=go_padj,
    p.adjust_rank=go_rank,
    nGenes = length(go_entrezID),
    nTRef_genes = sum(Tref_entrez %in% go_entrezID),
    nBRef_genes = sum(Tref_entrez %in% go_entrezID),
    stringsAsFactors = FALSE
  )
  
}))
outFile <- file.path(outFolder, "all_go_genes_dt.Rdata")
save(all_go_genes_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
load("LOOK_CONSERVEDREG_IMMUNEGO/all_go_genes_dt.Rdata")



