tad_dt_file <- file.path("ENCSR489OCU_NCI-H460_RANDOMNBRGENES_40kb", "genes2tad", "all_assigned_regions.txt")
tad_dt <- read.delim(tad_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo", "region", "start", "end"))



g2t_dt_file <- file.path("ENCSR489OCU_NCI-H460_RANDOMNBRGENES_40kb", "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)

nGenes_dt <- data.frame(
  nGenes = as.numeric(table(g2t_dt$region)),
  region = names(table(g2t_dt$region)),
  stringsAsFactors = FALSE
)

nGenes_dt <- nGenes_dt[order(nGenes_dt$nGenes, decreasing = T),]
# > head(nGenes_dt)
# nGenes       region
# 1542     35  chr14_TAD45
# 1472     28 chr14_TAD231
# 1469     24 chr14_TAD229
# 1785     23  chr15_TAD30

# > tad_dt[tad_dt$region == "chr14_TAD45",]
# chromo      region    start      end
# 1373  chr14 chr14_TAD45 22960001 23000000

hicds="ENCSR489OCU_NCI-H460_40kb"

g2t_dt_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
tad_g2t_dt <- g2t_dt[grep("_TAD", g2t_dt$region),]

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


# > sub_dt[235:250,]
# entrezID chromo    start      end   assembly strand symbol start_bin
# 10691    28694  chr14 22944306 22944365 GRCh37.p13      + TRAJ61  22960000
# 10692    28695  chr14 22945296 22945352 GRCh37.p13      + TRAJ60  22960000
# 10693    28696  chr14 22945543 22945596 GRCh37.p13      + TRAJ59  22960000
# 10694    28697  chr14 22946696 22946758 GRCh37.p13      + TRAJ58  22960000
# 10695    28698  chr14 22947861 22947923 GRCh37.p13      + TRAJ57  22960000
# 10696    28699  chr14 22948510 22948571 GRCh37.p13      + TRAJ56  22960000
# 10697    28700  chr14 22950686 22950742 GRCh37.p13      + TRAJ55  22960000
# 10698    28701  chr14 22951276 22951335 GRCh37.p13      + TRAJ54  22960000
# 10699    28702  chr14 22951993 22952058 GRCh37.p13      + TRAJ53  22960000
# 10700    28703  chr14 22955216 22955284 GRCh37.p13      + TRAJ52  22960000
# 10701    28704  chr14 22956171 22956233 GRCh37.p13      + TRAJ51  22960000
# 10702    28705  chr14 22957581 22957640 GRCh37.p13      + TRAJ50  22960000
# 10703    28706  chr14 22958476 22958531 GRCh37.p13      + TRAJ49  22960000
# 10704    28707  chr14 22959479 22959541 GRCh37.p13      + TRAJ48  22960000
# 10705    28708  chr14 22961838 22961894 GRCh37.p13      + TRAJ47  22960000
# 10706    28709  chr14 22962389 22962451 GRCh37.p13      + TRAJ46  22960000

