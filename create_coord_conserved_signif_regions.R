
# Rscript create_coord_conserved_signif_regions
# Rscript create_coord_conserved_signif_regions norm_vs_tumor
# Rscript create_coord_conserved_signif_regions subtypes
# Rscript create_coord_conserved_signif_regions wt_vs_mut

# Rscript create_coord_conserved_signif_regions <cmpType>

require(doMC)
require(foreach)
registerDoMC(40)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  cmpType <- ""
  filePrefix <- ""
} else if(length(args) == 1) {
  cmpType <- args[1]  
  filePrefix <- paste0(cmpType, "_")
}else {
  stop("---error\n")
}


outFolder <- file.path("CREATE_COORD_CONSERVED_SIGNIF_REGIONS", cmpType)
dir.create(outFolder, recursive = TRUE)


inFile <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmpType, paste0(filePrefix, "conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
stopifnot(file.exists(inFile))
conserved_signif_tads <- get(load(inFile))

cr = names(conserved_signif_tads)[1]

conserved_signif_tads_coord <- foreach(cr = names(conserved_signif_tads)) %dopar% {
  all_tads <- conserved_signif_tads[[paste0(cr)]]
  all_starts_ends <- sapply(all_tads, function(x) {
    g2t_dt <- read.delim(file.path(dirname(dirname(x)), "genes2tad", "all_assigned_regions.txt"), header=FALSE, stringsAsFactors=FALSE, col.names=c("chromo", "region", "start", "end"))
    stopifnot(sum(g2t_dt$region == basename(x)) == 1)
    g2t_dt$start[g2t_dt$region == basename(x)]
    c(start=g2t_dt$start[g2t_dt$region == basename(x)], end=g2t_dt$end[g2t_dt$region == basename(x)], chromo = g2t_dt$chromo[g2t_dt$region == basename(x)])
  })
  chromo <- unique(all_starts_ends["chromo",])
  stopifnot(length(chromo) == 1)
  c(min_start=min(all_starts_ends["start",]),
    max_start=max(all_starts_ends["start",]),
    min_end=min(all_starts_ends["end",]),
    max_end=max(all_starts_ends["end",]),
    chromo = chromo
  )
}

names(conserved_signif_tads_coord) <- names(conserved_signif_tads)

outFile <- file.path(outFolder, paste0(filePrefix, "conserved_signif_tads_coord.Rdata"))
save(conserved_signif_tads_coord, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


