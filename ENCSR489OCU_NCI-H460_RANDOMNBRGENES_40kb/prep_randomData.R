
# CREATE TADS BY TAKING NUMBER OF GENES 

# prepare for each chromo the file that looks like 
# ../ENCSR489OCU_NCI-H460_40kb/FINAL_DOMAINS/ENCSR489OCU_NCI-H460_chr1_YL_40kb_final_domains.txt
# (3-column BED format, only domain coordinates, separate file for each chromo, no header)

# -> then I can run ./3_assign_genes.sh ENCSR489OCU_NCI-H460_RANDOMNBRGENES

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


g2t_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
tad_g2t_dt <- g2t_dt[grep("_TAD", g2t_dt$region),]


meanNbrGenes <- round(mean(as.numeric(table(tad_g2t_dt$region))))

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)

stopifnot(setequal(gff_dt$entrezID[!gff_dt$chromo %in% c("chrX", "chrY")], g2t_dt$entrezID))

gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$start, gff_dt$end),]

gff_dt$start_bin <- round( gff_dt$start/binSize) * binSize


foo <- foreach(chr = unique(tad_g2t_dt$chromo)) %dopar% {
  
  sub_dt <- gff_dt[gff_dt$chromo == chr,]
  
  end_chromo_bin <- ceiling( max(sub_dt$end)/binSize) * binSize
  

  newStarts <- sub_dt$start_bin[seq(from=1, to=nrow(sub_dt), by=meanNbrGenes)]
  stopifnot(diff(newStarts) >= 0)
  newStarts <- newStarts[!duplicated(newStarts)]
  stopifnot(diff(newStarts) > 0)
  
  if(all(end_chromo_bin > newStarts)) {
    stopifnot(end_chromo_bin > newStarts)
    
    new_dt <- data.frame(
      chromo = chr,
      start = newStarts+1,
      end = c(newStarts[2:length(newStarts)], end_chromo_bin),
      stringsAsFactors = FALSE
    )
  } else {
    stopifnot(end_chromo_bin > newStarts[1:(length(newStarts)-1)])
    
    new_dt <- data.frame(
      chromo = chr,
      start = newStarts[1:(length(newStarts)-1)]+1,
      end = c(newStarts[2:length(newStarts)]),
      stringsAsFactors = FALSE
    )
  }

  
  stopifnot(new_dt$end > new_dt$start)
  stopifnot(new_dt$end %% binSize == 0)
  stopifnot(new_dt$start %% binSize == 1)
  
  outFile <- file.path(outFolder, paste0(rd_hicds, "_", chr, "_YL_", binSize/1000, "kb_final_domains.txt"))
  write.table(new_dt, file = outFile, sep="\t", col.names=F, row.names=F, quote=F, append=F )
  cat(paste0("... written: ", outFile, "\n"))
  
}



