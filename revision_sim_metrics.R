# copied 17.05.21 from /mnt/etemp/marie/tadcomp_shinyApp_feb21/tadcomp/scripts
# copied get_moc from /mnt/etemp/marie/TADcall_protocol/fig2_fig3_fig4_fig5_moc_calc.R


#' Calculate Measure of Concordance (MoC); adapted from formula given in Pfitzner et al. 2009
#'
# 

#' @param file1 Path to file for 1st partition (3-column BED format chromo-start-end; no header)
#' @param file2 Path to file for 2nd partition (3-column BED format chromo-start-end; no header)
#' @param chrSize Chromo size [bp]
#' @param nCpu Nbr of cpu available [use of foreach]
#' @param fillInter Consider inter-TAD gaps as TADs ? [in our analyses: TRUE]
#' @param meanWithInter Compute MoC for the inter-TADs separately and take the mean of TADs and interTADs MoCs [in our analyses: FALSE]
#' @param correctClust Add penalty term for the number of domains ? [in our analyses: FALSE]
#' @param binSize Bin size [bp]
#' @param noMinusOne Do not substract 1 to the MoC ? [in our analyses: FALSE]
#' @param fillInterGapZero Set score to 0 if comparing 1 TAD with interTAD ? [in our analyses: FALSE]

# requirement: doMC, foreach

get_MoC <- function(file1, file2, chrSize, nCpu = 1, fillInter=TRUE, meanWithInter = FALSE,
                    correctClust = FALSE, binSize = NA, noMinusOne = FALSE, fillInterGapZero = FALSE) {

  library(foreach)
  library(doMC)

  if(correctClust)
    if(is.na(binSize))
      stop("should provide binSize for correctClust!\n")

  if(fillInterGapZero)
    fillInter <- TRUE
  
  registerDoMC(cores=nCpu)
  
  if(fillInter & meanWithInter)
    stop("not meaningful")
  
  if(typeof(file1) == "character") {
    if (file.info(as.character(file1))$size == 0) {
      set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
      colnames(set1DT) <- c("chromo", "start", "end")
    }
  } else {
    set1DT <- file1
    stopifnot(ncol(set1DT) == 3)
    colnames(set1DT) <- c("chromo", "start", "end")
    
  }

  if(typeof(file2) == "character") {
    if (file.info(as.character(file2))$size == 0) {
      set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
      colnames(set2DT) <- c("chromo", "start", "end")
    }
    
  } else {
    set2DT <- file2
    stopifnot(ncol(set2DT) == 3)
    colnames(set2DT) <- c("chromo", "start", "end")
    
  }
  
  if(fillInterGapZero){
    set1DT_nofilled <- set1DT
    set2DT_nofilled <- set2DT
  }
  
  if(fillInter) {
    set1DT <- fill_part(set1DT, chrSize)
    set2DT <- fill_part(set2DT, chrSize)
  }
  
  if(fillInterGapZero){
    if(nrow(set1DT) > 0)
    set1DT$partType <- unlist(sapply(1:nrow(set1DT), function(x)
                ifelse(any(set1DT_nofilled$start == set1DT$start[x] & set1DT_nofilled$end == set1DT$end[x]), "domain", "gap")))
    if(nrow(set2DT) > 0)
    set2DT$partType <- unlist(sapply(1:nrow(set2DT), function(x)
      ifelse(any(set2DT_nofilled$start == set2DT$start[x] & set2DT_nofilled$end == set2DT$end[x]), "domain", "gap")))
  }
  
  if(nrow(set1DT) > 0)
    rownames(set1DT) <- paste0("P", 1:nrow(set1DT))
  if(nrow(set2DT) > 0)
    rownames(set2DT) <- paste0("Q", 1:nrow(set2DT))
  
  if(nrow(set1DT) == 0 & nrow(set2DT) > 0){
    MoC_score <- 0
  } else if(nrow(set1DT) > 0 & nrow(set2DT) == 0){
    MoC_score <- 0
  } else if( (nrow(set1DT) == nrow(set2DT))  & ( all(set1DT$start == set2DT$start) & all(set1DT$end == set2DT$end) ) ) {
      MoC_score <- 1
  } else if( (nrow(set1DT) == nrow(set2DT)) & (nrow(set1DT) == 1) ) {
      MoC_score <- 1
  } else {
    all_fragmentDT <- foreach(i = 1:nrow(set1DT), .combine='rbind') %dopar% {
      ref_start <- set1DT$start[i]
      ref_end <- set1DT$end[i]
      # if exact match
      all_matches <- which((set2DT$start == ref_start & set2DT$end == ref_end)  |
                             # nested
                             (set2DT$start >= ref_start & set2DT$end <= ref_end) |
                             # overlap left
                             (set2DT$start <= ref_start & set2DT$end >= ref_start) |
                             # overlap right
                             (set2DT$start <= ref_end & set2DT$end >= ref_end))
      all_matches_2 <- which(set2DT$end >= ref_start & set2DT$start <= ref_end)
      stopifnot(all(all_matches == all_matches_2))
      
      fragmentDT <- foreach(i_match = all_matches, .combine='rbind') %dopar% {
        # "fragment" size
        tmp_range <- c(set2DT$start[i_match] :set2DT$end[i_match])
        tmp_range <- tmp_range[tmp_range >= ref_start & tmp_range <= ref_end]
        frag_overlap <- tmp_range[length(tmp_range)] - tmp_range[1] + 1
        c(rownames(set1DT)[i], rownames(set2DT)[i_match], frag_overlap)
      }
      fragmentDT
    }
  
    # there is no intersect -> 0
    if(is.null(all_fragmentDT)) {
      MoC_score <- 0
    } else {
      # through matrix otherwise drop if nrow=1
      all_fragmentDT <- as.data.frame(matrix(all_fragmentDT, ncol=3), stringsAsFactors=F)
      rownames(all_fragmentDT) <- NULL
      colnames(all_fragmentDT) <- c("set1", "set2", "intersect_size")
      all_fragmentDT$intersect_size <- as.numeric(as.character(all_fragmentDT$intersect_size))
      
      MoC_score <- foreach(i = 1:nrow(all_fragmentDT), .combine = 'sum') %dopar% {
        d1 <- all_fragmentDT$set1[i]
        d2 <- all_fragmentDT$set2[i]
        # get the size of the domain from P  
        Pi <- set1DT[d1, "end"]  - set1DT[d1, "start"] + 1
        # get the size of the domain from Q 
        Qi <- set2DT[d2, "end"]  - set2DT[d2, "start"] + 1
        Fij <- all_fragmentDT$intersect_size[i]
        moc <- ((Fij*Fij)/(Pi * Qi))
        if(fillInterGapZero) { 
           if(set1DT[all_fragmentDT$set1[i], "partType" ] != set2DT[all_fragmentDT$set2[i], "partType" ])
              moc <- 0          
        }
        moc
      }
      # if wanted, correct for the number of clusters
      if(correctClust) {
        # penalize number of clusters / number tot of bins
        penaltyTerm <- nrow(all_fragmentDT) / (chrSize/binSize)
        MoC_score <- MoC_score + (1 - penaltyTerm) 
      }
    }
    
	if(! noMinusOne)
      MoC_score <- MoC_score - 1

    MoC_score <- MoC_score/(sqrt(nrow(set1DT) * nrow(set2DT)) - 1) 
  }
  if(meanWithInter) {
    bd_only_set1DT <- bd_only(set1DT, chrSize = chrSize)
    file1_BD_only <- sub("_final_domains.txt", "_final_domains_BD_only.txt", file1)
    write.table(bd_only_set1DT, file = file1_BD_only, col.names = F, row.names = F, quote=F, sep="\t")
    
    bd_only_set2DT <- bd_only(set2DT, chrSize = chrSize)
    file2_BD_only <- sub("_final_domains.txt", "_final_domains_BD_only.txt", file2)
    write.table(bd_only_set2DT, file = file2_BD_only, col.names = F, row.names = F, quote=F, sep="\t")
    
      MoC_score_inter <- get_MoC(file1_BD_only, file2_BD_only,chrSize, fillInter=fillInter, meanWithInter = FALSE )
    MoC_score <- (MoC_score + MoC_score_inter)/2
  }
  MoC_score
}


################
fill_part <- function(fillDT, chrSize) {
  stopifnot(all(colnames(fillDT) == c("chromo", "start", "end")))
  
  if(nrow(fillDT) == 0) {
    fillDT <- data.frame(chromo = "chr6",
                                  start = 1,
                                  end = chrSize)
    return(fillDT)
  }
  
  # add the intra-TAD also as domains
  fillDT_BD <- foreach(x = 1:nrow(fillDT), .combine='rbind') %dopar% {
    if(x == 1) {
      if(fillDT$start[1] == 1) {
        tmpDT <- data.frame(chromo = fillDT$chromo[1],
                            start = fillDT$start[1],
                            end = fillDT$end[1])
      } else{
        tmpDT <- data.frame(chromo = c(fillDT$chromo[1],fillDT$chromo[1]),
                            start = c(1, fillDT$start[1]),
                            end = c(fillDT$start[1]-1,fillDT$end[1] ))
      }
    } else{
      if(fillDT$start[x] > fillDT$end[x-1] + 1) {
        tmpDT <- data.frame(chromo = c(fillDT$chromo[x], fillDT$chromo[x]),
                            start = c(fillDT$end[x-1]+1, fillDT$start[x]),
                            end = c(fillDT$start[x] -1, fillDT$end[x]))
      } else{
        tmpDT <- data.frame(chromo = fillDT$chromo[x],
                            start = fillDT$start[x],
                            end = fillDT$end[x])
      }
    }
  }
  fillDT_BD <- as.data.frame(fillDT_BD)
  fillDT_BD$start <- as.numeric(as.character(fillDT_BD$start))
  fillDT_BD$end <- as.numeric(as.character(fillDT_BD$end))
  lastrow <- nrow(fillDT_BD)
  if(fillDT_BD$end[lastrow] < chrSize ){
    endDT <- data.frame(chromo = fillDT_BD$chromo[lastrow],
                        start = fillDT_BD$end[lastrow] + 1,
                        end = chrSize)
    fillDT <- rbind(fillDT_BD, endDT)
  } else{
    fillDT <- fillDT_BD
  }
  
  stopifnot(all( fillDT$end > fillDT$start  ))
  stopifnot(!any(duplicated(fillDT$start)))
  stopifnot(!any(duplicated(fillDT$end)))
  return(fillDT)
}



# 
#' Calculate "Jaccard index" based on concordance of each bin
#'
#' @param file1 Path to file for 1st partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param file2 Path to file for 2nd partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param binSize Bin size [bp]
#' @param chrSize Chromo size [bp]
#' @param nCpu Nbr of cpu available [use of foreach]

# options(repos = BiocInstaller::biocinstallRepos()) # for shiny - version 2019
# getOption("repos") # for shiny  - version 2019
options(repos = BiocManager::repositories())  # - version 2020

# requirement: foreach, doMC (if nCpu > 1)

get_bin_JaccardIndex <- function(file1, file2, binSize, chrSize=NULL, nCpu = 1){
  
  library(foreach)
  if(nCpu > 1) {
    library(doMC)
    registerDoMC(cores=nCpu)
  }
  
  if(typeof(file1) == "character") {
    if (file.info(file1)$size == 0) {
      set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
      colnames(set1DT) <- c("chromo", "start", "end")
    }
  } else {
    set1DT <- file1
    stopifnot(ncol(set1DT) == 3)
    colnames(set1DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set1DT$start))
  stopifnot(is.numeric(set1DT$end))
  
  if(typeof(file2) == "character") {
    if (file.info(file2)$size == 0) {
      set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
      colnames(set2DT) <- c("chromo", "start", "end")
    }
  } else {
    set2DT <- file2
    stopifnot(ncol(set2DT) == 3)
    colnames(set2DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set2DT$start))
  stopifnot(is.numeric(set2DT$end))

  stopifnot(length(unique(set1DT$chromo)) == 1)
  stopifnot(length(unique(set2DT$chromo)) == 1)
  stopifnot(unique(set2DT$chromo) == unique(set1DT$chromo))
  
  cat(paste0("# of domains in set1\t=\t", nrow(set1DT), "\n"))
  cat(paste0("# of domains in set2\t=\t", nrow(set2DT), "\n"))
  
  
  if(is.null(chrSize)) {
    chrSize <- max(c(set1DT$end, set2DT$end))
  }
  
  # check if the mid-position of the bin is within domain or not
  # "to" should be round up to bin size
  # e.g. if binSize=10 and chrSize=101 -> should go up to 110
  all_bins_mid <- seq(from=binSize/2, to=ceiling(chrSize/binSize)*binSize, by=binSize)
  
  stopifnot(!all_bins_mid %in% set1DT$start)
  stopifnot(!all_bins_mid %in% set1DT$end)
  stopifnot(!all_bins_mid %in% set2DT$start)
  stopifnot(!all_bins_mid %in% set2DT$end)
  pos=all_bins_mid[1]
  matchingVect <- foreach(pos=all_bins_mid, .combine='c') %dopar% {
    as.numeric(any(set1DT$start < pos & set1DT$end > pos) == any(set2DT$start < pos & set2DT$end > pos))
  }
  return(mean(matchingVect))
}
##################################################################################################################################################################################################  
##################################################################################################################################################################################################
##################################################################################################################################################################################################
  
#' Calculate "Jaccard index" for the boundaries 
#'
#' @param file1 Path to file for 1st partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param file2 Path to file for 2nd partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param tolRad Tolerance radius for boundary matching [bp]
#' @param matchFor Matching for set1 ("set1") or set2 ("set2") only or for both ("all"; default)
#' @param nCpu Nbr of cpu available [use of foreach]

# requirement: GenomicRanges

get_boundaries_JaccardIndex <- function(file1, file2, tolRad, matchFor="all", nCpu = 1){

  stopifnot(length(tolRad) == 1)
  stopifnot(is.numeric(tolRad))
  
  stopifnot(matchFor %in% c("set1", "set2", "all"))
  
  library(GenomicRanges)
  
  if(typeof(file1) == "character") {
    if (file.info(file1)$size == 0) {
      set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
      colnames(set1DT) <- c("chromo", "start", "end")
    }
  } else {
    set1DT <- file1
    stopifnot(ncol(set1DT) == 3)
    colnames(set1DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set1DT$start))
  stopifnot(is.numeric(set1DT$end))
  
  if(typeof(file2) == "character") {
    if (file.info(file2)$size == 0) {
      set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
      colnames(set2DT) <- c("chromo", "start", "end")
    }
  } else {
    set2DT <- file2
    stopifnot(ncol(set2DT) == 3)
    colnames(set2DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set2DT$start))
  stopifnot(is.numeric(set2DT$end))
  
  stopifnot(length(unique(set1DT$chromo)) == 1)
  stopifnot(length(unique(set2DT$chromo)) == 1)
  stopifnot(unique(set2DT$chromo) == unique(set1DT$chromo))
  
  mychr <- unique(set1DT$chromo)
  
  cat(paste0("# of domains in set1\t=\t", nrow(set1DT), "\n"))
  cat(paste0("# of domains in set2\t=\t", nrow(set2DT), "\n"))

  # boundaries matching -> treat starts and ends equally:
  # (handle 1-based and 0-based start coordinates)
  if(unique(set1DT$start %%10) == 1) { 
    bd1_DT <- data.frame(chromo=mychr, bdPos=c(set1DT$start-1, set1DT$end))
  } else {
    bd1_DT <- data.frame(chromo=mychr, bdPos=c(set1DT$start, set1DT$end))
  }
  bd1_DT <- unique(bd1_DT)
  bd1_DT <- bd1_DT[order(bd1_DT$bdPos),]
  if(unique(set2DT$start %%10) == 1) { 
    bd2_DT <- data.frame(chromo=mychr, bdPos=c(set2DT$start-1, set2DT$end))
  } else {
    bd2_DT <- data.frame(chromo=mychr, bdPos=c(set2DT$start, set2DT$end))
  }
  bd2_DT <- unique(bd2_DT)
  bd2_DT <- bd2_DT[order(bd2_DT$bdPos),]
  cat(paste0("# of boundaries in set1\t=\t", nrow(bd1_DT), "\n"))
  cat(paste0("# of boundaries in set2\t=\t", nrow(bd2_DT), "\n"))

  
  if(matchFor %in% c("set1", "all")) {
    query1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd1_DT$bdPos, end=bd1_DT$bdPos))
    ref2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd2_DT$bdPos-tolRad, end=bd2_DT$bdPos+tolRad))
    match1 <- query1_GR %over% ref2_GR
    stopifnot(length(match1) == length(query1_GR))
    stopifnot(length(match1) == nrow(bd1_DT))
    if(matchFor == "set1") return(mean(match1))
  }
  if(matchFor %in% c("set2", "all")) {
    query2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd2_DT$bdPos, end=bd2_DT$bdPos))
    ref1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bd1_DT$bdPos-tolRad, end=bd1_DT$bdPos+tolRad)) 
    match2 <- query2_GR %over% ref1_GR
    stopifnot(length(match2) == length(query2_GR))
    stopifnot(length(match2) == nrow(bd2_DT))
    if(matchFor == "set2") return(mean(match2))
  }
  # if come to here -> matchFor=="all"
  return(mean(c(match1, match2)))
}

##################################################################################################################################################################################################  
##################################################################################################################################################################################################
##################################################################################################################################################################################################

#' Calculate ratio of matching TAD
#'
#' @param file1 Path to file for 1st partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param file2 Path to file for 2nd partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param coverMatchRatioThresh Ratio to be covered to be considered matching [0-1]
#' @param matchFor Matching for set1 ("set1") or set2 ("set2") only or for both ("all"; default)
#' @param nCpu Nbr of cpu available [use of foreach]
#' 
#' WARNING: coverMatchRatioThresh=0 -> will not return 1 -> will return ratio of TADs that have a match, without considering a threshold of bp matching !

# requirement: GenomicRanges

get_ratioMatchingTADs <- function(file1, file2, coverMatchRatioThresh, matchFor="all", nCpu = 1){
  
  stopifnot(matchFor %in% c("set1", "set2", "all"))
  
  stopifnot(coverMatchRatioThresh >= 0 & coverMatchRatioThresh <= 1)
  
  library(GenomicRanges)
  
  if(typeof(file1) == "character") {
    if (file.info(file1)$size == 0) {
      set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
      colnames(set1DT) <- c("chromo", "start", "end")
    }
  } else {
    set1DT <- file1
    stopifnot(ncol(set1DT) == 3)
    colnames(set1DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set1DT$start))
  stopifnot(is.numeric(set1DT$end))
  
  if(typeof(file2) == "character") {
    if (file.info(file2)$size == 0) {
      set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
      colnames(set2DT) <- c("chromo", "start", "end")
    }
  } else {
    set2DT <- file2
    stopifnot(ncol(set2DT) == 3)
    colnames(set2DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set2DT$start))
  stopifnot(is.numeric(set2DT$end))
  
  stopifnot(length(unique(set1DT$chromo)) == 1)
  stopifnot(length(unique(set2DT$chromo)) == 1)
  stopifnot(unique(set2DT$chromo) == unique(set1DT$chromo))
  
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
  
  if(matchFor %in% c("set1", "all")) {
    ref_GR <- set1_GR
    query_GR <- set2_GR
    ref_tad_size_set1 <- setNames(c(set1DT$end-set1DT$start+1), c(set1DT$region))  # +1 since starts are "1-based" for GenomicRanges
    # determine which features from the query overlap which features in the subject
    overlap_GR <- findOverlaps(query=query_GR, subject=ref_GR)
    if(length(overlap_GR) == 0) {  
      set1_match <- 0
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
    }
    if(matchFor == "set1") return(mean(set1_match))
  }
  
  if(matchFor %in% c("set2", "all")) {
    ref_GR <- set2_GR
    query_GR <- set1_GR
    ref_tad_size_set2 <- setNames(c(set2DT$end-set2DT$start+1), c(set2DT$region)) # +1 since starts are "1-based" for GenomicRanges
    # determine which features from the query overlap which features in the subject
    overlap_GR <- findOverlaps(query=query_GR, subject=ref_GR)
    if(length(overlap_GR) == 0) {  
      set2_match <- 0
    } else {
      IDoverlap_hits_all <- findOverlaps(query=query_GR,
                                         subject=query_GR)
      IDoverlap_hits <- IDoverlap_hits_all[queryHits(IDoverlap_hits_all) != subjectHits(IDoverlap_hits_all)]
      
      refID <- names(ref_GR[subjectHits(overlap_GR)])
      queryID <- names(query_GR[queryHits(overlap_GR)])
      
      stopifnot(refID %in% set2DT$region)
      stopifnot(refID %in% names(ref_tad_size_set2))
      
      set2_overlapDT_bp <- data.frame(
        refID = refID,
        queryID = queryID,
        overlapBp = width(pintersect(ref_GR[refID], query_GR[queryID])),
        overlapBpRatio = width(pintersect(ref_GR[refID], query_GR[queryID]))/ref_tad_size_set2[refID],
        stringsAsFactors = FALSE)
      # retain only the matching that pass the threshold
      set2_overlapDT_bp <- set2_overlapDT_bp[set2_overlapDT_bp$overlapBpRatio >= coverMatchRatioThresh,]
      set2_match <- set2DT$region %in% set2_overlapDT_bp$refID
      stopifnot(length(set2_match) == nrow(set2DT))
    }
    if(matchFor == "set2") return(mean(set2_match))
  }
  # if come to here -> matchFor=="all"
  return(mean(c(set1_match, set2_match)))
  
}

##################################################################################################################################################################################################  
##################################################################################################################################################################################################
##################################################################################################################################################################################################

#' Calculate variation of information
#'
#' @param file1 Path to file for 1st partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param file2 Path to file for 2nd partition (3-column BED format chromo-start-end; no header) or 3-column dataframe object
#' @param chrSize Chromo size [bp]
#' @param gapsAsDomains Handle gaps as domains
#' 
#' # requirement: GenomicRanges, get_MoC.R (for fill_part function)

get_variationInformation <- function(file1, file2, chrSize=NULL, gapsAsDomains=FALSE){
  
  library(GenomicRanges)

  if(typeof(file1) == "character") {
    if (file.info(file1)$size == 0) {
      set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
      colnames(set1DT) <- c("chromo", "start", "end")
    }
  } else {
    set1DT <- file1
    stopifnot(ncol(set1DT) == 3)
    colnames(set1DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set1DT$start))
  stopifnot(is.numeric(set1DT$end))
  
  if(typeof(file2) == "character") {
    if (file.info(file2)$size == 0) {
      set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
    } else {
      set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
      colnames(set2DT) <- c("chromo", "start", "end")
    }
  } else {
    set2DT <- file2
    stopifnot(ncol(set2DT) == 3)
    colnames(set2DT) <- c("chromo", "start", "end")
  }
  stopifnot(is.numeric(set2DT$start))
  stopifnot(is.numeric(set2DT$end))
  
  stopifnot(length(unique(set1DT$chromo)) == 1)
  stopifnot(length(unique(set2DT$chromo)) == 1)
  stopifnot(unique(set2DT$chromo) == unique(set1DT$chromo))
  
  mychr <- unique(set1DT$chromo)
  
  cat(paste0("# of domains in set1\t=\t", nrow(set1DT), "\n"))
  cat(paste0("# of domains in set2\t=\t", nrow(set2DT), "\n"))
  
  #' NB: for correct TAD matching using GenomicRanges, start should be "1-based" and end "0-based", otherwise would find an overlap of 1 between the domains here below:
  # chr6 1425000 1600000 set1_chr6_TAD5
  # chr6 775000 1425000 set2_chr6_TAD1
  # (handle 1-based and 0-based start coordinates) # for TAD size calc. # ensure correct 1-based start for GenomicRange overlap calc.
  if(unique(set1DT$start %%10) == 0) set1DT$start <- set1DT$start+1
  if(unique(set2DT$start %%10) == 0) set2DT$start <- set2DT$start+1
  
  
  if(is.null(chrSize)) {
    chrSize <- max(c(set1DT$end, set2DT$end))
  }
  
  if(gapsAsDomains){
    source("get_MoC.R")
    set1DT <- fill_part(set1DT, chrSize)
    set2DT <- fill_part(set2DT, chrSize)  
  }
  

  
  set1DT$size <- set1DT$end-set1DT$start+1
  set2DT$size <- set2DT$end-set2DT$start+1
  

  
  
  
  set1DT$region <- paste0("set1_", mychr, "_TAD", 1:nrow(set1DT))
  set1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=set1DT$start, end=set1DT$end, names=set1DT$region))
  
  set2DT$region <- paste0("set2_", mychr, "_TAD", 1:nrow(set2DT))
  set2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=set2DT$start, end=set2DT$end, names=set2DT$region))
  

  set1_ref_GR <- set1_GR
  set1_query_GR <- set2_GR
  ref_tad_size_set1 <- setNames(c(set1DT$size), c(set1DT$region)) # +1 since starts are "1-based" for GenomicRanges
  # determine which features from the query overlap which features in the subject
  set1_overlap_GR <- findOverlaps(query=set1_query_GR, subject=set1_ref_GR)
  if(length(set1_overlap_GR) == 0) {  
    set1_overlapDT_bp <- data.frame(
      refID = character(0),
      queryID = character(0),
      overlapBp = numeric(0),
      overlapBpRatio = numeric(0),
      stringsAsFactors = FALSE)
    
  } else {
    set1_IDoverlap_hits_all <- findOverlaps(query=set1_query_GR,
                                       subject=set1_query_GR)
    set1_IDoverlap_hits <- set1_IDoverlap_hits_all[queryHits(set1_IDoverlap_hits_all) != subjectHits(set1_IDoverlap_hits_all)]
    
    set1_refID <- names(set1_ref_GR[subjectHits(set1_overlap_GR)])
    set1_queryID <- names(set1_query_GR[queryHits(set1_overlap_GR)])
    
    stopifnot(set1_refID %in% set1DT$region)
    stopifnot(set1_refID %in% names(ref_tad_size_set1))
    
    set1_overlapDT_bp <- data.frame(
      refID = set1_refID,
      queryID = set1_queryID,
      overlapBp = width(pintersect(set1_ref_GR[set1_refID], set1_query_GR[set1_queryID])),
      overlapBpRatio = width(pintersect(set1_ref_GR[set1_refID], set1_query_GR[set1_queryID]))/ref_tad_size_set1[set1_refID],
      stringsAsFactors = FALSE)
  }
  
  
  
    set2_ref_GR <- set2_GR
    set2_query_GR <- set1_GR
    ref_tad_size_set2 <- setNames(c(set2DT$size), c(set2DT$region)) # +1 since starts are "1-based" for GenomicRanges
    # determine which features from the query overlap which features in the subject
    set2_overlap_GR <- findOverlaps(query=set2_query_GR, subject=set2_ref_GR)
    if(length(set2_overlap_GR) == 0) {  
      set2_overlapDT_bp <- data.frame(
        refID = character(0),
        queryID = character(0),
        overlapBp = numeric(0),
        overlapBpRatio = numeric(0),
        stringsAsFactors = FALSE)
    } else {
      set2_IDoverlap_hits_all <- findOverlaps(query=set2_query_GR,
                                              subject=set2_query_GR)
      set2_IDoverlap_hits <- set2_IDoverlap_hits_all[queryHits(set2_IDoverlap_hits_all) != subjectHits(set2_IDoverlap_hits_all)]
      
      set2_refID <- names(set2_ref_GR[subjectHits(set2_overlap_GR)])
      set2_queryID <- names(set2_query_GR[queryHits(set2_overlap_GR)])
      
      stopifnot(set2_refID %in% set2DT$region)
      stopifnot(set2_refID %in% names(ref_tad_size_set2))
      
      set2_overlapDT_bp <- data.frame(
        refID = set2_refID,
        queryID = set2_queryID,
        overlapBp = width(pintersect(set2_ref_GR[set2_refID], set2_query_GR[set2_queryID])),
        overlapBpRatio = width(pintersect(set2_ref_GR[set2_refID], set2_query_GR[set2_queryID]))/ref_tad_size_set2[set2_refID],
        stringsAsFactors = FALSE)
      
    }
    # VI definition from Filippova et al. 2014
    # https://almob.biomedcentral.com/articles/10.1186/1748-7188-9-14 
    # calculate entropy for set1
    # set1DT$domain_p <- set1DT$size/chrSize          #### ???????????????????????????
    set1DT$domain_p <- set1DT$size/sum(set1DT$size)
    set1_H <- sum(set1DT$domain_p * log10(set1DT$domain_p))
    
    # calculate entropy for set2
    # set2DT$domain_p <- set2DT$size/chrSize         #### ???????????????????????????
    set2DT$domain_p <- set2DT$size/sum(set2DT$size)
    set2_H <- sum(set2DT$domain_p * log10(set2DT$domain_p))
    
    # calculate mutual information
    set1_p <- setNames(set1DT$domain_p, set1DT$region)
    set2_p <- setNames(set2DT$domain_p, set2DT$region)
    set1_set2_p <- c(set1_p, set2_p)

    stopifnot(nrow(set1_overlapDT_bp) == nrow(set2_overlapDT_bp))
    stopifnot(sum(set1_overlapDT_bp$overlapBp)==sum(set2_overlapDT_bp$overlapBp)    )
    
    # all_overlap_dt <- rbind(set1_overlapDT_bp, set2_overlapDT_bp) # they are the same !!! reciprocal matching !!!
    all_overlap_dt <-  set2_overlapDT_bp
    stopifnot(all_overlap_dt$refID %in% names(set1_set2_p))
    stopifnot(all_overlap_dt$queryID %in% names(set1_set2_p))
    
    # all_overlap_dt$domain_pij <- all_overlap_dt$overlapBp/chrSize
    all_overlap_dt$domain_pij <- all_overlap_dt$overlapBp/max(c(sum(set1DT$size), sum(set2DT$size)))
    all_overlap_dt$domain_pi <- set1_set2_p[paste0(all_overlap_dt$refID)]
    all_overlap_dt$domain_pj <- set1_set2_p[paste0(all_overlap_dt$queryID)]
    
    set1_set2_MI <- sum(all_overlap_dt$domain_pij * log10(all_overlap_dt$domain_pij/(all_overlap_dt$domain_pi*all_overlap_dt$domain_pj)))
    
    VI <- set1_H+set2_H-2*set1_set2_MI
    
    return(VI)
}





######################################################################## trial toy stuff
# 
# tolRad = 10*2
# 
# bddt1 = data.frame(chr="chr6",
#                  bdPos = c(0, 30, 50,100,140,180))
# 
# bddt2  = data.frame(chr="chr6",
#                    bdPos = c(10, 80, 140, 190, 220))
# 
# 
# query1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bddt1$bdPos, end=bddt1$bdPos))
# ref1_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bddt1$bdPos-tolRad, end=bddt1$bdPos+tolRad))
# 
# query2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bddt2$bdPos, end=bddt2$bdPos))
# ref2_GR <- GRanges(seqnames=mychr, ranges=IRanges(start=bddt2$bdPos-tolRad, end=bddt2$bdPos+tolRad))
# 
# 
# query1_GR %over% ref1_GR
# query2_GR %over% ref2_GR
# 
# 
# query1_GR %over% ref2_GR
# query2_GR %over% ref1_GR








                   
