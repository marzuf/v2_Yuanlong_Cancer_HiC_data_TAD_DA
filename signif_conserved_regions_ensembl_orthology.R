
# Rscript signif_conserved_regions_ensembl_orthology.R
# Rscript signif_conserved_regions_ensembl_orthology.R norm_vs_tumor
# Rscript signif_conserved_regions_ensembl_orthology.R subtypes
# Rscript signif_conserved_regions_ensembl_orthology.R wt_vs_mut

# Rscript signif_conserved_regions_ensembl_orthology.R <cmpType>

script_name <- "signif_conserved_regions_ensembl_orthology.R"

cat("> START ", script_name, "\n")

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

myHeight <- ifelse(plotType=="png", 400, 7)
myWidth  <- ifelse(plotType=="png", 600, 9)

buildData <- TRUE

settingSuffix <- "tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3"


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

outFolder <- file.path("SIGNIF_CONSERVED_REGIONS_ENSEMBL_ORTHOLOGY", cmpType)
dir.create(outFolder, recursive = TRUE)


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

inFile <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmpType, paste0("plot_matching_dt_", settingSuffix, ".Rdata"))
cat(inFile,"\n")
stopifnot(file.exists(inFile))
conserved_by_ds_dt <- get(load(inFile))

stopifnot(setequal(colnames(conserved_by_ds_dt), names(coord_conserved_signif_tads)))
stopifnot(setequal(colnames(conserved_by_ds_dt), names(conserved_signif_tads)))


inFile <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmpType, paste0(filePrefix, "conserved_signif_", settingSuffix, ".Rdata"))
cat(inFile,"\n")
stopifnot(file.exists(inFile))
conserved_signif_tads <- get(load(inFile))

inFile <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmpType, paste0(filePrefix, "conserved_regions_with_genes_signif_", settingSuffix, ".Rdata"))
cat(inFile,"\n")
stopifnot(file.exists(inFile))
conserved_with_genes_dt <- get(load(inFile))

orth_dt <- get(load("DATA_ENSEMBL_ORTHOLOGY/ensembl_orthology_prep_data.Rdata"))
orth_dt$entrezID <- as.character(orth_dt$entrezID)

all_orth_sp <- c("Mouse", "Chicken", "Dog", "Chimpanzee")

i=1
if(buildData) {
  conserved_regions_orderScores_dt <- foreach(i = 1:nrow(conserved_with_genes_dt), .combine='rbind') %dopar% {
    
    reg_entrezID <- unlist(strsplit(x=conserved_with_genes_dt$intersect_genes_entrez[i], split=","))
    
    orth_sub_dt <- orth_dt[orth_dt$entrezID %in% reg_entrezID,]
    
    conserv_scores <- sapply(all_orth_sp, function(x) mean(orth_sub_dt[,paste0(x, ".Gene.order.conservation.score")], na.rm=T))
    names(conserv_scores) <- paste0(all_orth_sp, "_geneOrderConservScore")
    
    cbind(conserved_with_genes_dt[i,], t(data.frame(conserv_scores)))
  }
  
  outFile <- file.path(outFolder, paste0(filePrefix, "conserved_regions_orderScores_", settingSuffix, ".Rdata"))
  save(conserved_regions_orderScores_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))  
} else {
  outFile <- file.path(outFolder, paste0(filePrefix, "conserved_regions_orderScores_", settingSuffix, ".Rdata"))
  conserved_regions_orderScores_dt <- get(load(outFile))
}

nMatches <- colSums(plot_matching_dt)

stopifnot(conserved_regions_orderScores_dt$conserved_region %in% names(nMatches))

conserved_regions_orderScores_dt$nConservMatches <- nMatches[paste0(conserved_regions_orderScores_dt$conserved_region)]

sp = all_orth_sp[1]
fo <- foreach(sp = all_orth_sp) %dopar% {
  
  sp_conserv_orth_dt <- aggregate( as.formula(paste0(sp, "_geneOrderConservScore ~ conserved_region + nConservMatches")), data = conserved_regions_orderScores_dt, FUN=mean, na.rm=TRUE)
  stopifnot(!duplicated(sp_conserv_orth_dt$conserved_region))
  
  outFile <- file.path(outFolder, paste0(filePrefix, sp, "_geneOrderConservScore_by_nbrConservMatch_", settingSuffix, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot( as.formula(paste0(sp, "_geneOrderConservScore ~ nConservMatches")), data=sp_conserv_orth_dt, 
           xlab = "# conserved datasets", ylab="Gene order conservation score", main = paste0("Orth. query species: ", sp))
  mtext(side=3, text = paste0("# conserved regions = ", length(unique(sp_conserv_orth_dt$conserved_region))))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
  
}






##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))









