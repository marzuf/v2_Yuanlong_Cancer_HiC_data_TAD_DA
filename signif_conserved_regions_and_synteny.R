
# Rscript signif_conserved_regions_and_synteny.R
# Rscript signif_conserved_regions_and_synteny.R norm_vs_tumor
# Rscript signif_conserved_regions_and_synteny.R subtypes
# Rscript signif_conserved_regions_and_synteny.R wt_vs_mut

# Rscript signif_conserved_regions_and_synteny.R <cmpType>

cat("> START ", "signif_conserved_regions_and_synteny.R", "\n")

startTime <- Sys.time()

require(doMC)
require(foreach)
registerDoMC(40)
require(GenomicRanges)
require(ggplot2)
require(ggsci)


ggsci_pal <- "d3"
ggsci_subpal <- ""

plotType <- "png"
source("../FIGURES_V2_YUANLONG/settings.R")

myWidth <- myHeight <- ifelse(plotType=="png", 400, 7)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  cmpType <- ""
  filePrefix <- ""
  cmpTit <- paste0("all")
} else if(length(args) == 1) {
  cmpType <- args[1]  
  filePrefix <- paste0(cmpType, "_")
  cmpTit <- cmpType
}else {
  stop("---error\n")
}

outFolder <- file.path("SIGNIF_CONSERVED_REGIONS_AND_SYNTENY", cmpType)
dir.create(outFolder, recursive = TRUE)

inferCARsBlocks <- read.delim("inferCARs_data/Orthology.Blocks_processed_hg19.txt", header=TRUE, stringsAsFactors = FALSE)

proCARsBlocks <- read.delim("procars_orthology_blocks_processsed.txt", header=TRUE, stringsAsFactors = FALSE)
proCARsBlocks$genome[proCARsBlocks$genome == "homo_sapiens"] <- "hg19"


settingSuffix <- "tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3"

# inFile <- file.path("CREATE_COORD_CONSERVED_SIGNIF_REGIONS", cmpType, paste0(filePrefix, "conserved_signif_tads_coord.Rdata"))
inFile <- file.path("CREATE_COORD_CONSERVED_SIGNIF_REGIONS", cmpType, paste0(filePrefix, "conserved_signif_tads_coord.Rdata"))
coord_conserved_signif_tads <- get(load(inFile))


conserved_regions_dt <- data.frame(
    do.call(cbind, list(
      do.call(rbind, lapply(coord_conserved_signif_tads, function(x) x[["chromo"]])),
      do.call(rbind, lapply(coord_conserved_signif_tads, function(x) x[["min_start"]])),
      do.call(rbind, lapply(coord_conserved_signif_tads, function(x) x[["max_end"]]))
    )), stringsAsFactors = FALSE)
colnames(conserved_regions_dt) <- c("chromo", "start", "end")
conserved_regions_dt$start <- as.numeric(as.character(conserved_regions_dt$start))
conserved_regions_dt$end <- as.numeric(as.character(conserved_regions_dt$end))
stopifnot(is.numeric(conserved_regions_dt$start))
stopifnot(is.numeric(conserved_regions_dt$end))
stopifnot(!is.na(conserved_regions_dt$start))
stopifnot(!is.na(conserved_regions_dt$end))
conserved_regions_dt$region <- rownames(conserved_regions_dt)
stopifnot(conserved_regions_dt$end > conserved_regions_dt$start)

inFile <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmpType, paste0(filePrefix, "conserved_signif_", settingSuffix, ".Rdata"))
cat(inFile,"\n")
stopifnot(file.exists(inFile))
conserved_signif_tads <- get(load(inFile))

inFile <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmpType, paste0("plot_matching_dt_", settingSuffix, ".Rdata"))
cat(inFile,"\n")
stopifnot(file.exists(inFile))
conserved_by_ds_dt <- get(load(inFile))

stopifnot(setequal(colnames(conserved_by_ds_dt), names(coord_conserved_signif_tads)))
stopifnot(setequal(colnames(conserved_by_ds_dt), names(conserved_signif_tads)))


synType = "inferCARs"
for(synType in c("inferCARs", "proCARs")) {
  
  synBlocks_dt <- eval(parse(text=paste0(synType, "Blocks")))
  
  nSynteny_dt <- aggregate(genome~blockID, data=synBlocks_dt, FUN=function(x){length(unique(x))})
  nUniqueSynteny_dt <- aggregate(genome~blockID, data=synBlocks_dt, FUN=function(x){length(x[! x %in% x[duplicated(x)]])})
  
  colnames(nSynteny_dt)[colnames(nSynteny_dt) == "genome"] <- "nSyn"
  colnames(nUniqueSynteny_dt)[colnames(nUniqueSynteny_dt) == "genome"] <- "nUniqueSyn"
  
  nBlocks_dt <- merge(nUniqueSynteny_dt, nSynteny_dt, by="blockID")
  colnames(nBlocks_dt)[colnames(nBlocks_dt) == "blockID"] <- "queryID"
  
  hg19_synBlocks_dt <- synBlocks_dt[synBlocks_dt$genome == "hg19",]
  stopifnot(is.numeric(hg19_synBlocks_dt$start))
  stopifnot(is.numeric(hg19_synBlocks_dt$end))
  
  
  ### PREPARE THE BP OVERLAP
  cat("... preparing bp matching overlap\n")

  region_GRanges <-  GRanges(seqnames=conserved_regions_dt$chromo, 
                           ranges=IRanges(start=conserved_regions_dt$start, end=conserved_regions_dt$end, names=conserved_regions_dt$region))
  
  
  block_GRanges <-  GRanges(seqnames=hg19_synBlocks_dt$chromo, 
                          ranges=IRanges(start=hg19_synBlocks_dt$start, end=hg19_synBlocks_dt$end, names=hg19_synBlocks_dt$blockID))
  
  # determine which features from the query overlap which features in the subject
  overlap_GR <- findOverlaps(query=block_GRanges, subject=region_GRanges)
  
  if(length(overlap_GR) == 0) stop("NULL")
  
  refID <- names(region_GRanges[subjectHits(overlap_GR)])
  queryID <- names(block_GRanges[queryHits(overlap_GR)])
  
  overlapDT_bp <- data.frame(
    refID = refID,
    queryID = queryID,
    overlapBp = width(pintersect(region_GRanges[refID], block_GRanges[queryID])),
    # overlapBpRatio = width(pintersect(region_GRanges[refID], block_GRanges[queryID]))/ref_tad_size[refID],
    stringsAsFactors = FALSE)
  stopifnot(overlapDT_bp$refID %in% conserved_regions_dt$region)
  stopifnot(overlapDT_bp$queryID %in% hg19_synBlocks_dt$blockID)
  stopifnot(!is.na(overlapDT_bp))
  
  overlapDT_bp$entrezID <- gsub("(.+)_chr.+", "\\1", overlapDT_bp$refID)
  overlapDT_bp$regSymbol <- gsub("(.+)_chr.+", "\\1", overlapDT_bp$queryID)
  
  outFile <- file.path(outFolder, paste0(synType, "_overlapDT_bp.Rdata"))
  save(overlapDT_bp, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  nConserv_dt <- data.frame(
    refID = names(colSums(plot_matching_dt)),
    nConserv = as.numeric(colSums(plot_matching_dt)),
    stringsAsFactors = FALSE)
  
  
  conserv_and_match_dt <- nConserv_dt
  conserv_and_match_dt$isInSyn <- conserv_and_match_dt$refID %in% overlapDT_bp$refID
  t1_dt <- aggregate(isInSyn~nConserv , data = conserv_and_match_dt, FUN = function(x) sum(x)/length(x))
  t2_dt <- aggregate(isInSyn~nConserv , data = conserv_and_match_dt, FUN = function(x) sum(!x)/length(x))
  bar_dt <- data.frame(
    nConserv = c(t1_dt$nConserv, t2_dt$nConserv),
    fractRegion = c(t1_dt$isInSyn, t2_dt$isInSyn),
    fractType = c(rep("withSynt", nrow(t1_dt)), rep("noSynt", nrow(t2_dt))),
    stringsAsFacors=FALSE)
  
  xlabs <- as.character(sort((unique(bar_dt$nConserv))))
  bar_dt$nConserv <- factor(bar_dt$nConserv, levels=xlabs)
  
  fract_conserv_synt <- ggplot(bar_dt, aes(x=nConserv, y=fractRegion, fill=fractType, color=fractType)) + 
    geom_bar(position="stack", stat="identity") +
    coord_cartesian(expand = FALSE) +
    ggtitle("Synteny of conserved regions", subtitle = paste0("(", synType, ")"))+
    labs(fill="", color="")+
    eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
    eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
    scale_y_continuous(name=paste0("Fraction of regions from synt. blocks"),
                       breaks = scales::pretty_breaks(n = 10))+
    scale_x_discrete(name="# conservation", labels = xlabs)+
    theme( # Increase size of axis lines
      plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14, family=fontFamily),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"),
      axis.line.x = element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .3, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
      axis.text.x = element_text(color="black", hjust=0.5,vjust = 0.5, size=12, family=fontFamily),
      axis.title.y = element_text(color="black", size=14, family=fontFamily),
      axis.title.x = element_text(color="black", size=14, family=fontFamily),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.background =  element_rect(),
      legend.key = element_blank(),
      legend.title = element_text(face="bold", family=fontFamily)
    )
  
  
  
  outFile <- file.path(outFolder, paste0(filePrefix, synType, "_nbrConserv_and_fractSynteny_barplot.", plotType))
  ggsave(plot = fract_conserv_synt, filename = outFile, height=myHeightGG, width = myWidthGG*2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  tmp1 <- merge(nConserv_dt, overlapDT_bp, by="refID", all=TRUE)
  nConserv_nSyn_dt <- merge(tmp1, nBlocks_dt, by="queryID", all.x=TRUE, all.y=FALSE)
  
  nConserv_nSyn_dt$nUniqueSyn[is.na(nConserv_nSyn_dt$queryID)] <- 0
  nConserv_nSyn_dt$nSyn[is.na(nConserv_nSyn_dt$queryID)] <- 0
  
  stopifnot(!is.na(nConserv_nSyn_dt$nSyn))
  stopifnot(!is.na(nConserv_nSyn_dt$nUniqueSyn))
  
  outFile <- file.path(outFolder, paste0(filePrefix, synType, "_nbrConserv_nbrUniqueSyn.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(nConserv_nSyn_dt$nConserv ~ nConserv_nSyn_dt$nUniqueSyn,
       xlab = "# syn. spec. (unique blocks)",
       ylab = "# dataset conserv.",
       main = paste0("signif. conserv. and synteny - ", cmpTit),
       pch=16,
       cex.lab=plotCex,
       cex.axis=plotCex,
       cex.main=plotCex,
       cex=0.7)
  mtext(side=3, text = paste0(synType))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0(filePrefix, synType, "_nbrConserv_nbrSyn.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot( nConserv_nSyn_dt$nConserv ~nConserv_nSyn_dt$nSyn,
        xlab = "# syn. spec.",
        ylab = "# dataset conserv.",
        main = paste0("signif. conserv. and synteny - ", cmpTit),
        cex.lab=plotCex,
        cex.axis=plotCex,
        cex.main=plotCex,
       pch=16,
       cex=0.7)
  mtext(side=3, text = paste0(synType))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0(filePrefix, synType, "_nbrConserv_nbrSyn_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(nConserv ~ nSyn, data=nConserv_nSyn_dt,
          xlab = "# syn. spec.",
          ylab = "# dataset conserv.",
          main = paste0("signif. conserv. and synteny - ", cmpTit),
          cex.lab=plotCex,
          cex.axis=plotCex,
          cex.main=plotCex
          )
  mtext(side=3, text = paste0(synType))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0(filePrefix, synType, "_nbrConserv_nbrUniqueSyn_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(nConserv ~ nUniqueSyn, data=nConserv_nSyn_dt,
          xlab = "# syn. spec. (unique blocks)",
          ylab = "# dataset conserv.",
          main = paste0("signif. conserv. and synteny - ", cmpTit),
          cex.lab=plotCex,
          cex.axis=plotCex,
          cex.main=plotCex
  )
  mtext(side=3, text = paste0(synType))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}




cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))