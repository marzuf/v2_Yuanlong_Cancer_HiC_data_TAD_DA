startTime <- Sys.time()
require(ggplot2)

source("revision_settings.R")

require(GenomicRanges)

binSize <- 40* 10^3
match_tolRad <- 2*binSize
match_coverRatio <- 0.8

source("revision_sim_metrics.R")


# all_pairs <- file.path  ds1 / ds2 / exprds_norm_tumor

# Rscript revision_cellLines_sim_signif.R

outFolder <- file.path("REVISION_CELLLINES_SIM_SIGNIF")
dir.create(outFolder, recursive=TRUE )



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

onlyPipTADs <- FALSE
if(onlyPipTADs) stop("--not correctly implemented")

runFold <- "."

# rev_pairs <- file.path(  basename(dirname(all_pairs)), dirname(dirname(all_pairs)), basename(all_pairs))

dspair = all_pairs[1]

all_hicds <- unique(c(basename(dirname(all_pairs)), dirname(dirname(all_pairs))))


# all_pairs=all_pairs[1]

all_pairs_dt <- foreach(dspair = c(all_pairs), .combine='rbind') %dopar% {
  
  
  hicds1 <- basename(dirname(dspair)) 
  hicds2 <- dirname(dirname(dspair))
  exprds <- basename(dspair)
  
  tad2_file <- file.path(runFold, hicds2, "genes2tad", "all_assigned_regions.txt" )
  genome2_dt <- read.delim(tad2_file, sep="\t", header=FALSE, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
  tad2_dt <- genome2_dt[grepl("_TAD", genome2_dt$region),]
  if(onlyPipTADs) {
    pip_tads <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds2, exprds, "0_prepGeneData", "pipeline_regionList.Rdata")))
    stopifnot(pip_tads %in% tad2_dt$region)
    tad2_dt <- tad2_dt[tad2_dt$region %in% pip_tads,]
  }
  stopifnot(nrow(tad2_dt) > 0)

    
  tad1_file <- file.path(runFold, hicds1, "genes2tad", "all_assigned_regions.txt" )
  genome1_dt <- read.delim(tad1_file, sep="\t", header=FALSE, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
  tad1_dt <- genome1_dt[grepl("_TAD", genome1_dt$region),]
  if(onlyPipTADs) {
    pip_tads <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds1, exprds, "0_prepGeneData", "pipeline_regionList.Rdata")))
    stopifnot(pip_tads %in% tad1_dt$region)
    tad1_dt <- tad1_dt[tad1_dt$region %in% pip_tads,]
  }
  stopifnot(nrow(tad1_dt) > 0)
  
  
  hicds1_tads <- tad1_dt
  hicds2_tads <- tad2_dt
  
  
  stopifnot(grepl("_TAD", hicds1_tads$region))
  # hicds1_tads$region <- NULL
  stopifnot(grepl("_TAD", hicds2_tads$region))
  # hicds2_tads$region <- NULL
  
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
   
    
    dt_set1_boundaryMatch <- get_set1_boundaryMatch(set1DT=hicds1_chr_tads, 
                                               set2DT=hicds2_chr_tads, 
                                               tolRad=match_tolRad)
    
    head(dt_set1_boundaryMatch)
    
    # because they are converted
    dt_set1_boundaryMatch$start <- dt_set1_boundaryMatch$start + 1
    
    cat(paste0("nrow dt_set1_boundaryMatch = ", nrow(dt_set1_boundaryMatch), "\n"))
    
      dt_set1_tadCoverMatch <- get_set1_tadCoverMatch(set1DT=hicds1_chr_tads, 
                                                       set2DT=hicds2_chr_tads,
                                                       coverMatchRatioThresh= match_coverRatio) 
    
      head(dt_set1_tadCoverMatch)
        
      cat(paste0("nrow dt_set1_tadCoverMatch = ", nrow(dt_set1_tadCoverMatch), "\n"))
      
      # save(dt_set1_boundaryMatch, file="dt_set1_boundaryMatch", version=2)
      # save(dt_set1_tadCoverMatch, file="dt_set1_tadCoverMatch", version=2)
      match_dt1 <- merge(dt_set1_boundaryMatch, dt_set1_tadCoverMatch, by=c("chromo", "start", "end", "region"))
      # match_dt1 <- full_join(dt_set1_boundaryMatch, dt_set1_tadCoverMatch, by=c("chromo", "start", "end", "region"))
      
      cat(paste0("nrow match_dt1 = ", nrow(match_dt1), "\n"))
      
      
      stopifnot(nrow(match_dt1) == nrow(dt_set1_tadCoverMatch))
      stopifnot(nrow(match_dt1) == nrow(dt_set1_boundaryMatch))
      stopifnot(nrow(match_dt1) == nrow(hicds1_chr_tads))
      
      match_dt1$ref_ds <- hicds1
      match_dt1$match_ds <- hicds2
    
      ################## do the same but ds2 vs ds1
      
      
      
      dt_set2_boundaryMatch <- get_set1_boundaryMatch(set1DT=hicds2_chr_tads, 
                                                      set2DT=hicds1_chr_tads, 
                                                      tolRad=match_tolRad)
      
      head(dt_set2_boundaryMatch)
      
      # because they are converted
      dt_set2_boundaryMatch$start <- dt_set2_boundaryMatch$start + 1
      
      cat(paste0("nrow dt_set2_boundaryMatch = ", nrow(dt_set2_boundaryMatch), "\n"))
      
      dt_set2_tadCoverMatch <- get_set1_tadCoverMatch(set1DT=hicds2_chr_tads, 
                                                      set2DT=hicds1_chr_tads,
                                                      coverMatchRatioThresh= match_coverRatio) 
      
      head(dt_set2_tadCoverMatch)
      
      cat(paste0("nrow dt_set2_tadCoverMatch = ", nrow(dt_set2_tadCoverMatch), "\n"))
      
      # save(dt_set2_boundaryMatch, file="dt_set2_boundaryMatch", version=2)
      # save(dt_set2_tadCoverMatch, file="dt_set2_tadCoverMatch", version=2)
      match_dt2 <- merge(dt_set2_boundaryMatch, dt_set2_tadCoverMatch, by=c("chromo", "start", "end", "region"))
      # match_dt2 <- full_join(dt_set2_boundaryMatch, dt_set2_tadCoverMatch, by=c("chromo", "start", "end", "region"))
      
      cat(paste0("nrow match_dt2 = ", nrow(match_dt2), "\n"))
      
      
      stopifnot(nrow(match_dt2) == nrow(dt_set2_tadCoverMatch))
      stopifnot(nrow(match_dt2) == nrow(dt_set2_boundaryMatch))
      stopifnot(nrow(match_dt2) == nrow(hicds2_chr_tads))
      
      match_dt2$ref_ds <- hicds2
      match_dt2$match_ds <- hicds1
      
      
      match_dt <- rbind(match_dt1, match_dt2)
      match_dt$exprds <- exprds
      match_dt
      
    
  } # end iterating over chromo
  chr_dt
} # end iterating over pair

outFile <- file.path(outFolder, "all_pairs_dt.Rdata")
save(all_pairs_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

stopifnot(!duplicated(file.path(all_pairs_dt$chromo, all_pairs_dt$start, all_pairs_dt$end, 
                                all_pairs_dt$region, all_pairs_dt$ref_ds, all_pairs_dt$match_ds, all_pairs_dt$exprds)))
stopifnot(which(all_pairs_dt$end1_vs_end2_match) %in% which(all_pairs_dt$end1_vs_all2_match))   
stopifnot(all(which(all_pairs_dt$start1_vs_start2_match) %in% which(all_pairs_dt$start1_vs_all2_match))           )

cat("***Done\n")
cat(paste0(startTime , " - ", Sys.time(), "\n"))