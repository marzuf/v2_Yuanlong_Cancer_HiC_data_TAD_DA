startTime <- Sys.time()
cat(paste0("> Rscript Rscript conserved_region_genome_plot.R\n"))

# Rscript conserved_region_genome_plot.R
# Rscript conserved_region_genome_plot.R subtypes
# Rscript conserved_region_genome_plot.R norm_vs_tumor
# Rscript conserved_region_genome_plot.R wt_vs_mut

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

library(stringr)
pattern<-"connected_gene=.+?;"


signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3

args="norm_vs_tumor"
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
stopifnot(dir.exists(inFolder))
inFile <- file.path(inFolder, paste0(file_prefix,"conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_df.txt"))
stopifnot(file.exists(inFile))
conserved_dt <- read.delim(inFile, header=TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(conserved_dt) > 0)

outFolder <- file.path("CONSERVED_REGION_GENOME_PLOT", args)
dir.create(outFolder, recursive = TRUE)
# outfile <- file.path(outFolder, "tad_matching_enrich_enhancer_annot_genesonly.txt")
# file.remove(outfile)

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))


genome_nConserved_dt <- foreach(i = 1:nrow(conserved_dt), .combine='rbind') %dopar% {
  
  # cat("i = ", i, "\n")
  
  curr_reg <- conserved_dt$conserved_region[i]
  
  curr_tads <- conserved_dt$corresp_tads[i]
  
  
  curr_entrez <- unlist(strsplit(x=conserved_dt$intersect_genes_entrez[i], split=","))
  curr_symbols <- unlist(strsplit(x=conserved_dt$intersect_genes_symbol[i], split=","))
  
  
    
  all_regs <- unlist(strsplit(x=curr_tads, split=","))
  
  nMatch <- length(all_regs)
  
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
  
  data.frame(chromo=chromo,
             min_start = min_start,
             max_end = max_end,
             nMatch = nMatch, 
             stringsAsFactors = FALSE)
  
}

maxPosChromo <- aggregate(max_end~chromo, data=genome_nConserved_dt, FUN=max)

all_chromo <- paste0("chr", 1:22)

genome_nConserved_dt$chromo <- factor(genome_nConserved_dt$chromo, levels=all_chromo)
genome_nConserved_dt <- genome_nConserved_dt[order(genome_nConserved_dt$chromo),]


# 
# 
# genomePos_nConserv_dt <- foreach(i  = 1:nrow(maxPosChromo), .combine='rbind') %dopar% {
#   
#   chromo <- maxPosChromo$chromo[i]
#   
#   chr_region_dt <- genome_nConserved_dt[genome_nConserved_dt$chromo == chromo,]
#   stopifnot(nrow(chr_region_dt) > 0)
#   
#   chr_region_dt <- chr_region_dt[order(chr_region_dt$nMatch, decreasing = FALSE),] # in case,overwrite with max nMatch
#   
#   maxPos <- maxPosChromo$max_end[i]
#   
#   # maxPos <- 5000
#   
#   all_n_conserved <- rep(0, maxPos)
#   
#   for(j in 1:nrow(chr_region_dt)) {
#     cs <- chr_region_dt$min_start[j]
#     ce <- chr_region_dt$max_end[j]
#     nm <- chr_region_dt$nMatch[j]
#     all_n_conserved[cs:ce] <- nm
#   }
#   
#   data.frame(
#     chromo = chromo,
#     pos  = 1:length(all_n_conserved),
#     nConserv = all_n_conserved
#   )
# }

genomePos_nConserv_dt$genomeIdx <- 1:nrow(genomePos_nConserv_dt)


outFile <- file.path(outFolder, paste0("genomePosNconserv.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(NULL,
  xlim=range(genome_nConserved_dt$min_start),
  ylim = range(genome_nConserved_dt$nMatch)
)
i=0
for(chromo in genome_nConserved_dt$chromo) {
  i=i+1
  # points(x=genome_nConserved_dt$min_start[genome_nConserved_dt$chromo == chromo],
  #        y=genome_nConserved_dt$nMatch[genome_nConserved_dt$chromo == chromo],
  #        col=i)
  segments(x0=genome_nConserved_dt$min_start[genome_nConserved_dt$chromo == chromo],
         x1=genome_nConserved_dt$max_end[genome_nConserved_dt$chromo == chromo],
         y0=genome_nConserved_dt$nMatch[genome_nConserved_dt$chromo == chromo],
         y1=genome_nConserved_dt$nMatch[genome_nConserved_dt$chromo == chromo],
         col=i)
  
}
foo <- dev.off()
cat(paste0("... written: ", outfile, "\n"))

######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

