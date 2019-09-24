# Rscript prep_homer_file.R

require(foreach)
require(doMC)
registerDoMC(40)

outFolder <- file.path("HOMER")

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  inFolder <- "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2"
  args <- ""
  file_prefix <- ""
} else if(length(args) == 1) {
  inFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", args[1])
  args <- args[1]
  file_prefix <- paste0(args, "_")
} else {
  stop("error\n")
}


signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3

mainFolder <- "."

stopifnot(dir.exists(inFolder))
inFile <- file.path(inFolder, paste0(file_prefix,"conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_df.txt"))
stopifnot(file.exists(inFile))
conserved_dt <- read.delim(inFile, header=TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(conserved_dt) > 0)

#    Column1: Unique Peak ID
#    Column2: chromosome
#    Column3: starting position
#    Column4: ending position
#    Column5: Strand (+/- or 0/1, where 0="+", 1="-")


i=1
out_dt <- foreach(i = 1:nrow(conserved_dt), .combine='rbind') %dopar% {
  
  # cat("i = ", i, "\n")
  
  curr_reg <- conserved_dt$conserved_region[i]
  
  curr_tads <- conserved_dt$corresp_tads[i]
    
  all_regs <- unlist(strsplit(x=curr_tads, split=","))
  
  chromo <- unique(unlist(sapply(all_regs, function(reg) gsub("(.+)_.+", "\\1", basename(reg)))))
  stopifnot(length(chromo) == 1)
  
  reg=all_regs[[1]]
  all_starts_ends <- sapply(all_regs, function(reg) {
    hicds <- dirname(dirname(reg))
    tadfile <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    tad <- basename(reg)
    stopifnot(file.exists(tadfile))
    tad_dt <- read.delim(tadfile, stringsAsFactors = FALSE, header=F, col.names=c("chromo", "region", "start", "end"))
    stopifnot(tad %in% tad_dt$region)
    tad_start <- tad_dt$start[tad_dt$region==tad]
    stopifnot(length(tad_start) == 1)
    tad_end <- tad_dt$end[tad_dt$region==tad]
    stopifnot(length(tad_end) == 1)
    stopifnot( tad_dt$chromo[tad_dt$region==tad] == chromo)
    c(tad_start, tad_end)
  })
  stopifnot(ncol(all_starts_ends) == length(all_regs))
  stopifnot(nrow(all_starts_ends) == 2)
  stopifnot(all_starts_ends[1,]<all_starts_ends[2,])
  min_start <- min(all_starts_ends[1,])
  stopifnot(is.numeric(min_start))
  max_end <- max(all_starts_ends[2,])
  stopifnot(is.numeric(max_end))
  
  data.frame(
    region=curr_reg,
    chromo=chromo,
    start=min_start,
    end=max_end,
    stringsAsFactors = FALSE
  )
}



out_dt <- out_dt[,c("region", "chromo", "start", "end")]
out_dt <- out_dt[order(out_dt$chromo, out_dt$start, out_dt$end),]

out_dt$strand <- "+"

outFile <- file.path(outFolder, "conserved_regions", args, paste0(file_prefix, "conserved_regions_plus.txt"))
dir.create(dirname(outFile), recursive = TRUE)
write.table(out_dt, file=outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... written: ", outFile, "\n"))

out_dt$strand <- "-"
outFile <- file.path(outFolder, "conserved_regions", args, paste0(file_prefix, "conserved_regions_minus.txt"))
dir.create(dirname(outFile), recursive = TRUE)
write.table(out_dt, file=outFile, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
cat(paste0("... written: ", outFile, "\n"))










