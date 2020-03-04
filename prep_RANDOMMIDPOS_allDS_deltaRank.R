














####### => NEED TO BE IMPLEMENTED !!!!!!!!!!!!!!
####### => NEED TO BE IMPLEMENTED !!!!!!!!!!!!!!
####### => NEED TO BE IMPLEMENTED !!!!!!!!!!!!!!





















# RANDOMMIDPOS -> new set of TAD
# start and end mid position of adjacent TADs

# prepare for each chromo the file that looks like 
# ../ENCSR489OCU_NCI-H460_40kb/FINAL_DOMAINS/ENCSR489OCU_NCI-H460_chr1_YL_40kb_final_domains.txt
# (3-column BED format, only domain coordinates, separate file for each chromo, no header)

# -> then I can run ./3_assign_genes.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb

# Rscript prep_RANDOMMIDPOS_allDS_deltaRank.R

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

require(doMC)
require(foreach)
registerDoMC(40)

runFolder <- "."
binSize <- 40*10^3


all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds)]
all_hicds <- all_hicds[!grepl("PERMUT", all_hicds)]

hicds = "Barutcu_MCF-10A_40kb"

foo <- foreach(hicds = all_hicds) %dopar% {
  
    
  newFolder <- gsub("_40kb", "_RANDOMMIDPOS_40kb", hicds)
  
  outFolder <- file.path(newFolder,  "ASSIGNED_REGIONS_DELTARANK")
  dir.create(outFolder, recursive = TRUE)
  
  rd_hicds <- gsub("_40kb", "", newFolder)
  
  tad_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt")
  tad_dt <- read.delim(tad_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo","region", "start", "end"))
  
  onlyTAD_dt <- tad_dt[grep("_TAD", tad_dt$region),]
  
  
  onlyTAD_dt$midPos <- 0.5*(onlyTAD_dt$end + onlyTAD_dt$start)
  ### => get closest bin to midPos
  onlyTAD_dt$midPos_bin_ceiling <- abs(ceiling( onlyTAD_dt$midPos/binSize) * binSize - onlyTAD_dt$midPos)
  onlyTAD_dt$midPos_bin_floor<- abs(floor( onlyTAD_dt$midPos/binSize) * binSize - onlyTAD_dt$midPos)
  onlyTAD_dt$midPos_bin_rd <- abs(round( onlyTAD_dt$midPos/binSize) * binSize - onlyTAD_dt$midPos)
  onlyTAD_dt$midPos_check <- pmin(onlyTAD_dt$midPos_bin_ceiling, onlyTAD_dt$midPos_bin_floor)
  stopifnot(onlyTAD_dt$midPos_bin_rd == onlyTAD_dt$midPos_check)
  
  onlyTAD_dt$midPos_bin <- round( onlyTAD_dt$midPos/binSize) * binSize 
  
  
  chr = "chr1"
  
  all_chr_dt <- foreach(chr = unique(tad_dt$chromo), .combine='rbind') %dopar% {
    
    
    tadScore_dt_file <- file.path(runFolder, hicds, "FINAL_DOMAINS_WITH_SCORES", paste0(gsub("_40kb", "", hicds), "_", chr, "_YL_40kb_final_domains_with_scores.txt"))
    stopifnot(file.exists(tadScore_dt_file))
    tadScore_dt <- read.delim(tadScore_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo","start", "end", "rank"))
    
  
    
    sub_dt <- onlyTAD_dt[onlyTAD_dt$chromo == chr,c("chromo", "start", "end", "midPos_bin")]
    
    stopifnot(sub_dt$end > sub_dt$start)
    stopifnot(sub_dt$end %% binSize == 0)
    stopifnot(sub_dt$start %% binSize == 1)
    
    sub_dt_withScores <- merge(sub_dt, tadScore_dt, by=c("chromo", "start", "end"), all.x=TRUE, all.y=FALSE )
    stopifnot(!is.na(sub_dt_withScores))
    stopifnot(nrow(sub_dt_withScores) == nrow(sub_dt))
    
    stopifnot(is.numeric(sub_dt_withScores$midPos_bin))
    
    sub_dt_withScores <- sub_dt_withScores[order(sub_dt_withScores$midPos_bin),]
    
    new_tad_dt <- data.frame(
      chromo = chr,
      start = sub_dt_withScores$midPos_bin[1:(nrow(sub_dt_withScores)-1)]+1,
      end = sub_dt_withScores$midPos_bin[2:nrow(sub_dt_withScores)],
      abs_delta_rank = abs(diff(sub_dt_withScores$rank)),
      stringsAsFactors = FALSE
    )
    
    stopifnot(new_tad_dt$end > new_tad_dt$start)
    stopifnot(new_tad_dt$end %% binSize == 0)
    stopifnot(new_tad_dt$start %% binSize == 1)
    
    new_tad_dt
  }
  
  # outFile <- file.path(outFolder, "all_chr_dt.Rdata")
  # save(all_chr_dt, file=outFile, version=2)
  # cat(paste0("... written:" , outFile, "\n"))
  
  
  tad_noRank_file <- file.path(newFolder, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(tad_noRank_file))
  tad_noRank_dt <- read.delim(tad_noRank_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo","region", "start", "end"))
  
  assigned_withRank_dt <- merge(tad_noRank_dt, all_chr_dt, by=c("chromo", "start", "end"), all.x=TRUE, all.y=FALSE)
  stopifnot(nrow(assigned_withRank_dt) == nrow(tad_noRank_dt))
  stopifnot(!is.na(assigned_withRank_dt$abs_delta_rank[grepl("_TAD", assigned_withRank_dt$region)]))
  
  outFile <- file.path(outFolder, "all_assigned_regions_withDeltaRank.Rdata")
  save(assigned_withRank_dt, file=outFile, version=2)
  cat(paste0("... written:" , outFile, "\n"))
  
  
}



