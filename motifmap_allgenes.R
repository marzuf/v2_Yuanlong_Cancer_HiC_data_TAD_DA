# Rscript motifmap_allgenes.R

require(GenomicRanges)

outFolder <- "MOTIFMAP_ALLGENES"
dir.create(outFolder, recursive = TRUE)
motifFile <- "HUMAN_hg19_BBLS_1_00_FDR_0_10.bed"

buildData <- TRUE


all_hicds <- list.files(file.path("PIPELINE", "OUTPUT_FOLDER"))


if(buildData) {
  motif_dt <- read.delim(motifFile, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo", "start", "end", "motif", "FDR", "strand"))
  motif_dt$FDR <- NULL
  motif_dt$strand <- NULL
  motif_dt <- unique(motif_dt)
  
  setDir <- "/media/electron"
  setDir <- ""
  entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
  gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
  gff_dt$entrezID <- as.character(gff_dt$entrezID)
  stopifnot(!duplicated(gff_dt$entrezID))
  stopifnot(!duplicated(gff_dt$symbol))
  entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
  symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)
  
  ### PREPARE THE BP OVERLAP
  cat("... preparing bp matching overlap\n")
  gff_dt$entrezIDlong <- paste0(gff_dt$entrezID, "_", gff_dt$chromo, "_", gff_dt$start, "_", gff_dt$end)
  gene_GRanges <-  GRanges(seqnames=gff_dt$chromo, 
                           ranges=IRanges(start=gff_dt$start, end=gff_dt$end, names=gff_dt$entrezIDlong))
  
  motif_dt$motifID <- paste0(motif_dt$motif, "_", motif_dt$chromo, "_", motif_dt$start, "_", motif_dt$end)
  reg_GRanges <-  GRanges(seqnames=motif_dt$chromo, 
                          ranges=IRanges(start=motif_dt$start, end=motif_dt$end, names=motif_dt$motifID))
  
  # determine which features from the query overlap which features in the subject
  overlap_GR <- findOverlaps(query=reg_GRanges, subject=gene_GRanges)
  
  if(length(overlap_GR) == 0) stop("NULL")
  
  refID <- names(gene_GRanges[subjectHits(overlap_GR)])
  queryID <- names(reg_GRanges[queryHits(overlap_GR)])
  
  overlapDT_bp <- data.frame(
    refID = refID,
    queryID = queryID,
    overlapBp = width(pintersect(gene_GRanges[refID], reg_GRanges[queryID])),
    # overlapBpRatio = width(pintersect(gene_GRanges[refID], reg_GRanges[queryID]))/ref_tad_size[refID],
    stringsAsFactors = FALSE)
  stopifnot(overlapDT_bp$refID %in% gff_dt$entrezIDlong)
  stopifnot(overlapDT_bp$queryID %in% motif_dt$motifID)
  stopifnot(!is.na(overlapDT_bp))
  
  overlapDT_bp$entrezID <- gsub("(.+)_chr.+", "\\1", overlapDT_bp$refID)
  overlapDT_bp$regSymbol <- gsub("(.+)_chr.+", "\\1", overlapDT_bp$queryID)
  
  outFile <- file.path(outFolder, "overlapDT_bp.Rdata")
  save(overlapDT_bp, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
} else {
  inFile <- file.path(outFolder, "overlapDT_bp.Rdata")
  overlapDT_bp <- get(load(inFile))
  
}
require(doMC)
require(doMC)
require(dplyr)
registerDoMC(40)

foreach(hicds = all_hicds) %dopar% {
  
  g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
  g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
  
  g2t_reg_dt <- inner_join(g2t_dt, overlapDT_bp, by=c("entrezID"))
  
}






