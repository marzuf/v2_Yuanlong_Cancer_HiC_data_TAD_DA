startTime <- Sys.time()

source("revision_settings.R")

# Rscript revision_cellLines_sim_signif.R

outFolder <- file.path("REVISION_CELLLINES_SIM_SIGNIF")

require(foreach)
require(doMC)
registerDoMC(50)

all_hicds <- names(hicds_names)

runFold <- "."

bin_size <- 40*10^3
coverTADmatchRatio <-  0.8
boundariesJI_tolRad <- 2*bin_size


onlyPipTADs <- FALSE  # only implemented for normal vs. tumor !!!
if(onlyPipTADs) outFolder <- file.path("REVISION_CELLLINES_SIM_ONLYPIPTADS")
dir.create(outFolder, recursive = TRUE)


all_hicds <- all_hicds[1:3] 

hicds = all_hicds[1]

all_tads <- foreach(hicds = all_hicds) %dopar% {
  tad_file <- file.path(runFold, hicds, "genes2tad", "all_assigned_regions.txt" )
  genome_dt <- read.delim(tad_file, sep="\t", header=FALSE, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
  tad_dt <- genome_dt[grepl("_TAD", genome_dt$region),]
  
  if(onlyPipTADs) {
    hicds_exprds <- list.files(file.path("PIPELINE", "OUTPUT_FOLDER", hicds))
    hicds_exprds <- hicds_exprds[grepl("_norm_", hicds_exprds)]
    if(length(hicds_exprds) == 0) return(NULL)
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

all_pairs_dt <- foreach(i = 1:ncol(all_ds_pairs), .combine='rbind') %dopar% {
  
  
  hicds1 <- all_ds_pairs[1,i]
  hicds2 <- all_ds_pairs[2,i]
  
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