startTime <- Sys.time()

source("revision_settings.R")

require(GenomicRanges)

# tolRad in bp
get_set1_boundaryMatch <- function(set1DT, set2DT, tolRad, start1based=TRUE,
                                   match_type=c("bd_vs_bd", "tad_vs_tad", "tad_vs_bd")) {
  
  stopifnot(c("chromo", "start", "end") %in% colnames(set1DT))
  stopifnot(c("chromo", "start", "end") %in% colnames(set2DT) )
  
  mychr <- unique(set1DT$chromo)
  stopifnot(mychr == set2DT$chromo)
  stopifnot(length(mychr) == 1)
  
  if(start1based){
    set1DT$start <- set1DT$start-1 
    set2DT$start <- set2DT$start-1 
  } 
  
  # 
  # if(start1Based) { 
  #   
  #   all_bd1_DT <- data.frame(chromo=mychr, bdPos=c(set1DT$start-1, set1DT$end))
  #   all_bd2_DT <- data.frame(chromo=mychr, bdPos=c(set2DT$start-1, set2DT$end))
  # } else {
  #   all_bd1_DT <- data.frame(chromo=mychr, bdPos=c(set1DT$start, set1DT$end))
  #   all_bd2_DT <- data.frame(chromo=mychr, bdPos=c(set2DT$start, set2DT$end))
  # }
  
  if(match_type=="bd_vs_bd"){
    all_bd1_DT <- data.frame(chromo=mychr, bdPos=c(set1DT$start-1, set1DT$end))
    all_bd2_DT <- data.frame(chromo=mychr, bdPos=c(set2DT$start-1, set2DT$end))
    
  }
  if(match_type=="tad_vs_tad") {
    # start should match start, end should match end 
  }
  if(match_type="tad_vs_bd") {
    # start should match start or end and end should match start and end
  }
  
  

  
    query1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd1_DT$bdPos, end=bd1_DT$bdPos))
    ref2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd2_DT$bdPos-tolRad, end=bd2_DT$bdPos+tolRad))
    match1 <- query1_GR %over% ref2_GR
    stopifnot(length(match1) == length(query1_GR))
    stopifnot(length(match1) == nrow(bd1_DT))
  
}

get_set1_tadCoverMatch <- function(set1DT, set2DT,coverMatchRatioThresh ){
  stopifnot(c("chromo", "start", "end") %in% colnames(set1DT))
  stopifnot(c("chromo", "start", "end") %in% colnames(set2DT) )
  
  mychr <- unique(set1DT$chromo)
  
  cat(paste0("# of domains in set1\t=\t", nrow(set1DT), "\n"))
  cat(paste0("# of domains in set2\t=\t", nrow(set2DT), "\n"))
  
  
  #' NB: for correct TAD matching using GenomicRanges, start should be "1-based" and end "0-based", otherwise would find an overlap of 1 between the domains here below:
  # chr6 1425000 1600000 set1_chr6_TAD5
  # chr6 775000 1425000 set2_chr6_TAD1
  # (handle 1-based and 0-based start coordinates) # for TAD size calc. # ensure correct 1-based start for GenomicRange overlap calc.
  if(unique(set1DT$start %%10) == 0) set1DT$start <- set1DT$start+1
  if(unique(set2DT$start %%10) == 0) set2DT$start <- set2DT$start+1
  
  set1DT$region <- paste0("set1_", mychr, "_TAD", 1:nrow(set1DT))
  set1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=set1DT$start, end=set1DT$end, names=set1DT$region))
  
  set2DT$region <- paste0("set2_", mychr, "_TAD", 1:nrow(set2DT))
  set2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=set2DT$start, end=set2DT$end, names=set2DT$region))
  
  ref_GR <- set1_GR
  query_GR <- set2_GR
  ref_tad_size_set1 <- setNames(c(set1DT$end-set1DT$start+1), c(set1DT$region))  # +1 since starts are "1-based" for GenomicRanges
  # determine which features from the query overlap which features in the subject
  overlap_GR <- findOverlaps(query=query_GR, subject=ref_GR)
  if(length(overlap_GR) == 0) {  
    set1_match <- 0
    set1DT$set1_match <- FALSE
  } else {
    IDoverlap_hits_all <- findOverlaps(query=query_GR,
                                       subject=query_GR)
    IDoverlap_hits <- IDoverlap_hits_all[queryHits(IDoverlap_hits_all) != subjectHits(IDoverlap_hits_all)]
    
    refID <- names(ref_GR[subjectHits(overlap_GR)])
    queryID <- names(query_GR[queryHits(overlap_GR)])
    
    stopifnot(refID %in% set1DT$region)
    stopifnot(refID %in% names(ref_tad_size_set1))
    
    set1_overlapDT_bp <- data.frame(
      refID = refID,
      queryID = queryID,
      overlapBp = width(pintersect(ref_GR[refID], query_GR[queryID])),
      overlapBpRatio = width(pintersect(ref_GR[refID], query_GR[queryID]))/ref_tad_size_set1[refID],
      stringsAsFactors = FALSE)
    # retain only the matching that pass the threshold
    set1_overlapDT_bp <- set1_overlapDT_bp[set1_overlapDT_bp$overlapBpRatio >= coverMatchRatioThresh,]
    set1_match <- set1DT$region %in% set1_overlapDT_bp$refID
    stopifnot(length(set1_match) == nrow(set1DT))
    set1DT$set1_match <- set1_match
  }
  return(set1DT)
}




# all_pairs <- file.path  ds1 / ds2 / exprds_norm_tumor

# Rscript revision_cellLines_sim_signif.R

outFolder <- file.path("REVISION_CELLLINES_SIM_SIGNIF")

pvalthresh <- 0.01
runFolder <- "."
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
resultData$region_id <- file.path(resultData$dataset, resultData$region)
resultData$signif <- resultData$adjPvalComb <= pvalthresh


require(foreach)
require(doMC)
registerDoMC(50)

onlyPipTADs <- TRUE

runFold <- "."

dspair = all_pairs[1]

all_hicds <- unique(c(basename(dirname(all_pairs)), dirname(dirname(all_pairs))))


all_tads <- foreach(hicds = all_hicds) %dopar% {
  tad_file <- file.path(runFold, hicds, "genes2tad", "all_assigned_regions.txt" )
  genome_dt <- read.delim(tad_file, sep="\t", header=FALSE, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
  tad_dt <- genome_dt[grepl("_TAD", genome_dt$region),]
  
  if(onlyPipTADs) {
    hicds_exprds <- list.files(file.path("PIPELINE", "OUTPUT_FOLDER", hicds))
    hicds_exprds <- hicds_exprds[grepl("_norm_", hicds_exprds)]
    #  if(length(hicds_exprds) == 0) return(NULL) not here, i want norm tumor for all
    stopifnot(length(hicds_exprds) == 1)
    stopifnot(grepl(""))
    pip_tads <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, hicds_exprds, "0_prepGeneData", "pipeline_regionList.Rdata")))
    stopifnot(pip_tads %in% tad_dt$region)
    tad_dt <- tad_dt[tad_dt$region %in% pip_tads,]
  }
  
  stopifnot(nrow(tad_dt) > 0)
  tad_dt
}
names(all_tads) <- all_hicds
outFile <- file.path(outFolder, "all_tads.Rdata")
save(all_tads, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

all_ds_pairs <- combn(x=all_hicds, m=2)
stopifnot(nrow(all_ds_pairs) == 2)
stopifnot(ncol(all_ds_pairs) == length(all_hicds) * (length(all_hicds) -1) * 0.5)

all_pairs_dt <- foreach(dspair = all_pairs, .combine='rbind') %dopar% {
  
  
  hicds1 <- basename(dirname(all_pairs)) 
  hicds2 <- dirname(dirname(all_pairs))
  
  stopifnot(hicds1 %in% names(all_tads))
  stopifnot(hicds2 %in% names(all_tads))
  
  
  hicds1_tads <- all_tads[[paste0(hicds1)]]
  hicds2_tads <- all_tads[[paste0(hicds2)]]
  
  stopifnot(grepl("_TAD", hicds1_tads$region))
  hicds1_tads$region <- NULL
  stopifnot(grepl("_TAD", hicds2_tads$region))
  hicds2_tads$region <- NULL
  
  all_chromo <- unique(c(hicds1_tads$chromo, hicds2_tads$chromo))
  stopifnot(all_chromo %in% hicds1_tads$chromo)
  stopifnot(all_chromo %in% hicds2_tads$chromo)
  
  # all_chromo=all_chromo[1]
  
  chr_dt <- foreach(chr = all_chromo, .combine='rbind') %dopar% {
    
    cat(paste0("... start:\t", hicds1, " vs. ", hicds2, " - ", chr, "\n"))
    
    hicds1_chr_tads <- hicds1_tads[hicds1_tads$chromo == chr,]
    stopifnot(nrow(hicds1_chr_tads) > 0)
    
    hicds2_chr_tads <- hicds2_tads[hicds2_tads$chromo == chr,]
    stopifnot(nrow(hicds2_chr_tads) > 0)
    
    chrsize <- max(c(hicds1_chr_tads$end, hicds2_chr_tads$end))
    stopifnot(is.numeric(chrsize))
    
    cat(paste0("chrsize= ", chrsize, "\n"))
    
    # pour chaque TAD -> est-ce qu'il a un match
    # match d'1 boundary, match des 2 boundaries, cover match
    
    
    # outFile <- file.path(outFolder, "hicds1_chr_tads.Rdata")
    # save(hicds1_chr_tads, file=outFile, version=2)
    # cat(paste0("... written: ", outFile, "\n"))
    # outFile <- file.path(outFolder, "hicds2_chr_tads.Rdata")
    # save(hicds2_chr_tads, file=outFile, version=2)
    # cat(paste0("... written: ", outFile, "\n"))
    
    ds1_ds2_moc <- get_MoC(hicds1_chr_tads, hicds2_chr_tads, chrSize=chrsize)
    ds1_ds2_binJI <- get_bin_JaccardIndex(hicds1_chr_tads, hicds2_chr_tads, binSize=bin_size)
    ds1_ds2_boundJI <- get_boundaries_JaccardIndex(hicds1_chr_tads, hicds2_chr_tads, tolRad=boundariesJI_tolRad,matchFor="set1" )
    ds1_ds2_matchingTADs <- get_ratioMatchingTADs(hicds1_chr_tads, hicds2_chr_tads, coverMatchRatioThresh=coverTADmatchRatio, matchFor="set1" )
    ds1_ds2_VI <- get_variationInformation(hicds1_chr_tads, hicds2_chr_tads)
    
    ds2_ds1_boundJI <- get_boundaries_JaccardIndex(hicds2_chr_tads, hicds1_chr_tads, 
                                                   tolRad=boundariesJI_tolRad,matchFor="set1" )
    ds2_ds1_matchingTADs <- get_ratioMatchingTADs(hicds2_chr_tads, hicds1_chr_tads, 
                                                  coverMatchRatioThresh=coverTADmatchRatio, matchFor="set1" )
    ds2_ds1_VI <- get_variationInformation(hicds2_chr_tads, hicds1_chr_tads)
    
    stopifnot(ds1_ds2_VI == ds2_ds1_VI)
    
    # rbind(
    data.frame(
      ds1 = hicds1,
      ds2 = hicds2, 
      chromo = chr,
      moc =ds1_ds2_moc,
      binJI =ds1_ds2_binJI,
      VI = ds1_ds2_VI,
      ds1_ds2_boundJI =ds1_ds2_boundJI,
      ds1_ds2_matchingTADs = ds1_ds2_matchingTADs,
      ds2_ds1_boundJI =ds2_ds1_boundJI,
      ds2_ds1_matchingTADs = ds2_ds1_matchingTADs,
      # ds2_ds1_VI = ds2_ds1_VI,
      stringsAsFactors = FALSE
    )
    # ,data.frame(
    #   ds1 = hicds2,
    #   ds2 = hicds1, 
    #   moc =ds1_ds2_moc,  # symmetric
    #   binJI =ds1_ds2_binJI, # symmetric
    #   ds1_ds2_boundJI =ds2_ds1_boundJI,
    #   ds1_ds2_matchingTADs = ds2_ds1_matchingTADs,
    #   ds1_ds2_VI = ds2_ds1_VI,
    #   stringsAsFactors = FALSE
    # ))
  } # end iterating over chromo
  chr_dt
} # end iterating over pair

outFile <- file.path(outFolder, "all_pairs_dt.Rdata")
save(all_pairs_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

cat("***Done\n")
cat(paste0(startTime , " - ", Sys.time(), "\n"))