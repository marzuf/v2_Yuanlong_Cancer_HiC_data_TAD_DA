
# Rscript signif_conserved_regions_and_synteny
# Rscript signif_conserved_regions_and_synteny norm_vs_tumor
# Rscript signif_conserved_regions_and_synteny subtypes
# Rscript signif_conserved_regions_and_synteny wt_vs_mut

# Rscript signif_conserved_regions_and_synteny <cmpType>

cat("> START ", "signif_conserved_regions_and_synteny", "\n")

startTime <- Sys.time()

require(doMC)
require(foreach)
registerDoMC(40)
require(GenomicRanges)

plotType <- "png"
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

inFile <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmpType, paste0(filePrefix, "conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
cat(inFile,"\n")
stopifnot(file.exists(inFile))
conserved_signif_tads <- get(load(inFile))

inFile <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmpType, paste0("plot_matching_dt.Rdata"))
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
          main = paste0("signif. conserv. and synteny - ", cmpTit)
          )
  mtext(side=3, text = paste0(synType))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, paste0(filePrefix, synType, "_nbrConserv_nbrUniqueSyn_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(nConserv ~ nUniqueSyn, data=nConserv_nSyn_dt,
          xlab = "# syn. spec. (unique blocks)",
          ylab = "# dataset conserv.",
          main = paste0("signif. conserv. and synteny - ", cmpTit)
  )
  mtext(side=3, text = paste0(synType))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}




cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))