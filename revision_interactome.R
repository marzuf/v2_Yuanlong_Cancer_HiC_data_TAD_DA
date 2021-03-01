# load the interactome data
# assign the signif. interactions to the TADs

binSize <- 40*10^3

all_chrs <- paste0("chr", 1:22)

all_interactMatch_dt <- foreach(hicds = all_hicds, .combine= 'rbind') %dopar% {
  
  
  interactome_dt <- paste0(">>> CHANGE HERE")
  stopifnot(interactome_dt$binA %% 1 == 0)
  stopifnot(interactome_dt$binB %% 1 == 0)
  
  # ensure only intra
  stopifnot(interactome_dt$chromoA == interactome_dt$chromoB)
  
  
  all_domains_dt <- read.delim(file.path(hicds, "genes2tad/all_assigned_regions.txt"), col.names=c("chromo", "region", "start", "end"), 
                               header=FALSE, stringsAsFactors = FALSE)
  ### keep only the TADs !!!
  all_tads_dt <- all_domains_dt[grepl("_TAD", all_domains_dt$region),]
  stopifnot(nrow(all_tads_dt) > 0)
  
  # convert to 0-based bin
  all_tads_dt$startBin <- (all_tads_dt$start-1)/binSize
  stopifnot(all_tads_dt$startBin %% 1 == 0)
  all_tads_dt$endBin <- (all_tads_dt$end)/binSize-1
  stopifnot(all_tads_dt$endBin %% 1 == 0)
  
  stopifnot(all_chrs %in% all_tads_dt$chromo)
  
  interactome_dt$hicds <- hicds
  interactome_dt$binA_tadMatch <- NA
  interactome_dt$binB_tadMatch <- NA
  
  hicds_interactMatch_dt  <- foreach(i = 1:nrow(interactome_dt), .combine='rbind') %dopar% {
    
    # find the region of the 1st bin
    binA_tadMatch <- all_tads_dt$region[all_tads_dt$chromo <= interactome_dt$chromoA[i] & 
                                          all_tads_dt$startBin <= interactome_dt$binA[i] & 
                                          all_tads_dt$endBin >= interactome_dt$binA[i]]
    
    binB_tadMatch <- all_tads_dt$region[all_tads_dt$chromo <= interactome_dt$chromoB[i] & 
                                          all_tads_dt$startBin <= interactome_dt$binB[i] & 
                                          all_tads_dt$endBin >= interactome_dt$binB[i]]
    
    if(length(binA_tadMatch) == 0) binA_tadMatch <- NA
    if(length(binB_tadMatch) == 0) binB_tadMatch <- NA
    
    interactome_dt$binA_tadMatch <- binA_tadMatch
    interactome_dt$binB_tadMatch <- binB_tadMatch
    
  }
  
  
  
}
