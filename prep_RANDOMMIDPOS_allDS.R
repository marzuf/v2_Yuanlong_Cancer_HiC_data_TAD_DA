# RANDOMMIDPOS -> new set of TAD
# start and end mid position of adjacent TADs

# prepare for each chromo the file that looks like 
# ../ENCSR489OCU_NCI-H460_40kb/FINAL_DOMAINS/ENCSR489OCU_NCI-H460_chr1_YL_40kb_final_domains.txt
# (3-column BED format, only domain coordinates, separate file for each chromo, no header)

# -> then I can run ./3_assign_genes.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb

# Rscript prep_RANDOMMIDPOS_allDS.R

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

require(doMC)
require(foreach)
registerDoMC(40)

runFolder <- "."
binSize <- 40*10^3


all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("NCI-H460", all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds)]
all_hicds <- all_hicds[!grepl("PERMUT", all_hicds)]

hicds = "Barutcu_MCF-10A_40kb"
foo <- foreach(hicds = all_hicds) %dopar% {
  
    
  newFolder <- gsub("_40kb", "_RANDOMMIDPOS_40kb", hicds)
  stopifnot(!dir.exists(newFolder))
  outFolder <- file.path(newFolder,  "FINAL_DOMAINS")
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
  
  foo <- foreach(chr = unique(tad_dt$chromo)) %dopar% {
    
    sub_dt <- onlyTAD_dt[onlyTAD_dt$chromo == chr,c("chromo", "start", "end", "midPos_bin")]
    
    stopifnot(sub_dt$end > sub_dt$start)
    stopifnot(sub_dt$end %% binSize == 0)
    stopifnot(sub_dt$start %% binSize == 1)
    
    new_tad_dt <- data.frame(
      chromo = chr,
      start = sub_dt$midPos_bin[1:(nrow(sub_dt)-1)]+1,
      end = sub_dt$midPos_bin[2:nrow(sub_dt)],
      stringsAsFactors = FALSE
    )
    
    stopifnot(new_tad_dt$end > new_tad_dt$start)
    stopifnot(new_tad_dt$end %% binSize == 0)
    stopifnot(new_tad_dt$start %% binSize == 1)
    
    outFile <- file.path(outFolder, paste0(rd_hicds, "_", chr, "_YL_", binSize/1000, "kb_final_domains.txt"))
    write.table(new_tad_dt, file = outFile, sep="\t", col.names=F, row.names=F, quote=F, append=F )
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
  
  
  
  # call assign_genes
  mycmd <- paste0("./3_assign_genes.sh ", rd_hicds)
  cat(paste0("> ", mycmd, "\n"))
  system(mycmd)
}



