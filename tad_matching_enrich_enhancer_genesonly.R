startTime <- Sys.time()
cat(paste0("> Rscript enrich_enhancer_annot.R\n"))

# Rscript tad_matching_enrich_enhancer_genesonly.R
# Rscript tad_matching_enrich_enhancer_genesonly.R subtypes
# Rscript tad_matching_enrich_enhancer_genesonly.R norm_vs_tumor
# Rscript tad_matching_enrich_enhancer_genesonly.R wt_vs_mut

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)


library(stringr)
pattern<-"connected_gene=.+?;"



# enhancer data download from https://www.genecards.org/GeneHancer_version_4-4
enhancerFile <- "GeneHancer_version_4-4.csv"
stopifnot(file.exists(enhancerFile))
enhancerDT <- read.delim(enhancerFile, stringsAsFactors = FALSE, sep="\t", header=TRUE)
stopifnot(nrow(enhancerDT) > 0)


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

outFolder <- file.path("TAD_MATCHING_ENRICH_ENHANCER_ANNOT_GENESONLY", args)
dir.create(outFolder, recursive = TRUE)
outfile <- file.path(outFolder, "tad_matching_enrich_enhancer_annot_genesonly.txt")
file.remove(outfile)

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))

inFolder <- "ENRICH_SAME_DIFF_TAD"
inFile <- file.path(inFolder, "full_enhancerDT.Rdata")
stopifnot(file.exists(inFile))
full_enhancer_DT <- get(load(inFile))
full_enhancer_DT$entrezID <- as.character(full_enhancer_DT$entrezID)
i=2
for(i in 1:nrow(conserved_dt)) {
  
  # cat("i = ", i, "\n")
  
  curr_reg <- conserved_dt$conserved_region[i]
  
  curr_tads <- conserved_dt$corresp_tads[i]
  
  curr_entrez <- unlist(strsplit(x=conserved_dt$intersect_genes_entrez[i], split=","))
  curr_symbols <- unlist(strsplit(x=conserved_dt$intersect_genes_symbol[i], split=","))
  

  matching_gene_enhancer_dt <- full_enhancer_DT[full_enhancer_DT$entrezID %in% curr_entrez | 
                                                  full_enhancer_DT$geneSymbol %in% curr_symbols, ]
    
  
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
  
  matching_gene_enhancer_dt$enhancer_in_tad <- matching_gene_enhancer_dt$enhancer_start >= min_start &  matching_gene_enhancer_dt$enhancer_start <= max_end
  
  if(nrow(matching_gene_enhancer_dt) > 0) {
    
    inTAD_dt <- aggregate(enhancer_in_tad ~ geneSymbol, data=matching_gene_enhancer_dt, FUN=sum)
    colnames(inTAD_dt)[colnames(inTAD_dt) == "enhancer_in_tad"] <- "nEnhancersSameTAD"
    
    outTAD_dt <- aggregate(enhancer_in_tad ~ geneSymbol, data=matching_gene_enhancer_dt, FUN=function(x) sum(!x))
    colnames(outTAD_dt)[colnames(outTAD_dt) == "enhancer_in_tad"] <- "nEnhancersDiffTAD"
    
    stopifnot(sum(inTAD_dt$nEnhancersSameTAD) + sum(outTAD_dt$nEnhancersDiffTAD) == nrow(matching_gene_enhancer_dt))
    
    
    out_dt_tmp <- merge(inTAD_dt, outTAD_dt, by = "geneSymbol")  
  } else {
    out_dt_tmp <-  data.frame(geneSymbol = c(),
                              nEnhancersDiffTAD = c(),
                              nEnhancersSameTAD = c(),
                              stringsAsFactors = FALSE)
  }
  
  
  
  missingSymb <- curr_symbols[ ! (curr_symbols %in% out_dt_tmp$geneSymbol)]
  
  if(length(missingSymb) > 0 ) {
    missingDT <- data.frame(geneSymbol = missingSymb,
                            nEnhancersDiffTAD = NA,
                            nEnhancersSameTAD = NA,
                            stringsAsFactors = FALSE
    )
    out_dt <- rbind(out_dt_tmp, missingDT)
    
  } else {
    out_dt <- out_dt_tmp
  }
  
  out_dt <- out_dt[order(out_dt$nEnhancersSameTAD, out_dt$nEnhancersDiffTAD, decreasing = TRUE),]
  
  
  mycon <- file(outfile,"a")
  writeLines(text=paste0("\n*** ", i, ") ", curr_reg, "(", min_start, "-", max_end, ")"), con=mycon)
  close(mycon)
  
  
  write.table(out_dt, file = outfile, append=T, quote=F,col.names = TRUE,row.names=FALSE,sep="\t")
  
}
  
cat(paste0("... written: ", outfile, "\n"))

######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

