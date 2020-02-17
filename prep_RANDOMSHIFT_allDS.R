# SHIFT ALL START AND END POSITIONS WITH HALF TAD SIZE

# prepare for each chromo the file that looks like 
# ../ENCSR489OCU_NCI-H460_40kb/FINAL_DOMAINS/ENCSR489OCU_NCI-H460_chr1_YL_40kb_final_domains.txt
# (3-column BED format, only domain coordinates, separate file for each chromo, no header)

# -> then I can run ./3_assign_genes.sh ENCSR489OCU_NCI-H460_40kb_RANDOMSHIFT

# Rscript prep_RANDOMSHIFT_allDS.R

require(doMC)
require(foreach)
registerDoMC(40)

runFolder <- "."

binSize <- 40*10^3


setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)

gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$start, gff_dt$end),]

gff_dt$start_bin <- round( gff_dt$start/binSize) * binSize



all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("NCI-H460", all_hicds)]

hicds = "Barutcu_MCF-10A_40kb"
foo <- foreach(hicds = all_hicds) %dopar% {
  
  newFolder <- gsub("_40kb", "_RANDOMSHIFT_40kb", hicds)
  stopifnot(!dir.exists(newFolder))
  outFolder <- file.path(newFolder,  "FINAL_DOMAINS")
  dir.create(outFolder, recursive = TRUE)
  
  rd_hicds <- gsub("_40kb", "", newFolder)
  
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
  
  
  
  
  # call assign_genes
  mycmd <- paste0("./3_assign_genes.sh ", rd_hicds)
  cat(paste0("> ", mycmd, "\n"))
  system(mycmd)
  
  
  
}




