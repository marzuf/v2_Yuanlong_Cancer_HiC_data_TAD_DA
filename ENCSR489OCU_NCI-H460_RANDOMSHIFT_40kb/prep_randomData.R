# SHIFT ALL START AND END POSITIONS WITH HALF TAD SIZE

# prepare for each chromo the file that looks like 
# ../ENCSR489OCU_NCI-H460_40kb/FINAL_DOMAINS/ENCSR489OCU_NCI-H460_chr1_YL_40kb_final_domains.txt
# (3-column BED format, only domain coordinates, separate file for each chromo, no header)

# -> then I can run ./3_assign_genes.sh ENCSR489OCU_NCI-H460_40kb_RANDOMSHIFT

# Rscript prep_randomData.R

require(doMC)
require(foreach)
registerDoMC(40)

runFolder <- ".."
hicds <- "ENCSR489OCU_NCI-H460_40kb"
binSize <- 40*10^3

outFolder <- "FINAL_DOMAINS"
dir.create(outFolder, recursive = TRUE)

rd_hicds <- gsub("_40kb", "", basename(getwd()))

tad_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt")
tad_dt <- read.delim(tad_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo","region", "start", "end"))

onlyTAD_dt <- tad_dt[grep("_TAD", tad_dt$region),]

onlyTAD_dt$size <- onlyTAD_dt$end - onlyTAD_dt$start + 1

meanTADsize <- mean(onlyTAD_dt$size)
#meanTADsize*0.5 = 128642.2

# so that start and end positions of domains still concordant with genomic bins

shiftBp <- ceiling( (meanTADsize*0.5)/binSize) * binSize
#160000

onlyTAD_dt$start <- onlyTAD_dt$start + shiftBp
onlyTAD_dt$end <- onlyTAD_dt$end + shiftBp

foo <- foreach(chr = unique(tad_dt$chromo)) %dopar% {
  
  sub_dt <- onlyTAD_dt[onlyTAD_dt$chromo == chr,c("chromo", "start", "end")]
  
  stopifnot(sub_dt$end > sub_dt$start)
  stopifnot(sub_dt$end %% binSize == 0)
  stopifnot(sub_dt$start %% binSize == 1)
  
  outFile <- file.path(outFolder, paste0(rd_hicds, "_", chr, "_YL_", binSize/1000, "kb_final_domains.txt"))
  write.table(sub_dt, file = outFile, sep="\t", col.names=F, row.names=F, quote=F, append=F )
  cat(paste0("... written: ", outFile, "\n"))
  
  
}

