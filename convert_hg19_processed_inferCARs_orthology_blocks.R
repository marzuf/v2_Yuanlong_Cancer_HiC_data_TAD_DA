
# Rscript convert_hg19_processed_inferCARs_orthology_blocks.R

require(doMC)
require(foreach)

hg18_dt_raw <- read.delim(file.path("inferCARs_data/Orthology.Blocks_hg18.txt.bed"), stringsAsFactors = FALSE, sep=" ", header=FALSE)
hg18_dt <- foreach(i=1:nrow(hg18_dt_raw), .combine='rbind') %dopar% {
  rowCoord <- hg18_dt_raw[i,]
  data.frame(
    chromo = gsub("(chr.+):(.+)-(.+)", "\\1", rowCoord),
    start = gsub("(chr.+):(.+)-(.+)", "\\2", rowCoord),
    end = gsub("(chr.+):(.+)-(.+)", "\\3", rowCoord),
    stringsAsFactors = FALSE
  )
}


hg19_dt_raw <- read.delim(file.path("inferCARs_data/Orthology.Blocks_hg18.txt.bed_hg19conv.bed"), stringsAsFactors = FALSE, sep=" ", header=FALSE)
hg19_dt <- foreach(i=1:nrow(hg19_dt_raw), .combine='rbind') %dopar% {
  rowCoord <- hg19_dt_raw[i,]
  data.frame(
    chromo = gsub("(chr.+):(.+)-(.+)", "\\1", rowCoord),
    start = gsub("(chr.+):(.+)-(.+)", "\\2", rowCoord),
    end = gsub("(chr.+):(.+)-(.+)", "\\3", rowCoord),
    stringsAsFactors = FALSE
  )
}

failed_dt_raw <- read.delim(file.path("inferCARs_data/convert_failure_bed.txt"), stringsAsFactors = FALSE, sep=" ", header=FALSE)
failed_dt <- foreach(i=1:nrow(failed_dt_raw), .combine='rbind') %dopar% {
  rowCoord <- failed_dt_raw[i,]
  data.frame(
    chromo = gsub("(chr.+):(.+)-(.+)", "\\1", rowCoord),
    start = gsub("(chr.+):(.+)-(.+)", "\\2", rowCoord),
    end = gsub("(chr.+):(.+)-(.+)", "\\3", rowCoord),
    stringsAsFactors = FALSE
  )
}


hg18_dt$id <- file.path( hg18_dt$chromo, hg18_dt$start, hg18_dt$end)
hg19_dt$id <- file.path(hg19_dt$chromo, hg19_dt$start, hg19_dt$end)
failed_dt$id <- file.path(failed_dt$chromo, failed_dt$start, failed_dt$end)

stopifnot(failed_dt$id %in% hg18_dt$id)

converted_hg18_dt <- hg18_dt[!hg18_dt$id %in% failed_dt$id,]

stopifnot(nrow(converted_hg18_dt) == nrow(hg19_dt))

colnames(converted_hg18_dt) <- paste0(colnames(converted_hg18_dt), "_hg18")
colnames(hg19_dt) <- paste0(colnames(hg19_dt), "_hg19")

hg18_hg19_match_dt <- cbind(converted_hg18_dt, hg19_dt)

hg18_hg19_match <- setNames(hg18_hg19_match_dt$id_hg19, hg18_hg19_match_dt$id_hg18)


hg18_orthoBlocks_dt <- read.delim(file.path("inferCARs_data/Orthology.Blocks_processed.txt"), stringsAsFactors = FALSE, header=TRUE)
hg18_orthoBlocks_dt$id <- file.path(hg18_orthoBlocks_dt$chromo, hg18_orthoBlocks_dt$start, hg18_orthoBlocks_dt$end)
hg18_orthoBlocks_dt$matchID <- hg18_hg19_match[paste0(hg18_orthoBlocks_dt$id)]
hg18_orthoBlocks_dt$matchID[hg18_orthoBlocks_dt$genome != "hg18"] <- NA
hg18_orthoBlocks_dt$matchStart <- basename(dirname(hg18_orthoBlocks_dt$matchID))
hg18_orthoBlocks_dt$matchEnd <- basename(hg18_orthoBlocks_dt$matchID)
hg18_orthoBlocks_dt$matchGenome <- ifelse(is.na(hg18_orthoBlocks_dt$matchID), NA, "hg19")


hg19_orthoBlocks_dt <- hg18_orthoBlocks_dt


hg19_orthoBlocks_dt$genome <- ifelse(is.na(hg19_orthoBlocks_dt$matchGenome), hg19_orthoBlocks_dt$genome, hg19_orthoBlocks_dt$matchGenome)
hg19_orthoBlocks_dt$start <- ifelse(is.na(hg19_orthoBlocks_dt$matchStart), hg19_orthoBlocks_dt$start, hg19_orthoBlocks_dt$matchStart)
hg19_orthoBlocks_dt$end <- ifelse(is.na(hg19_orthoBlocks_dt$matchEnd), hg19_orthoBlocks_dt$end, hg19_orthoBlocks_dt$matchEnd)

hg19_orthoBlocks_dt$newID <- file.path(hg19_orthoBlocks_dt$chromo, hg19_orthoBlocks_dt$start, hg19_orthoBlocks_dt$end)

stopifnot(hg19_dt$id_hg19 %in% hg19_orthoBlocks_dt$newID[hg19_orthoBlocks_dt$genome == "hg19"])
stopifnot(hg19_orthoBlocks_dt$id[hg19_orthoBlocks_dt$genome == "hg18"] %in% failed_dt$id)


hg19_orthoBlocks_dt <- hg19_orthoBlocks_dt[hg19_orthoBlocks_dt$genome != "hg18",]

newBlocks <- aggregate(genome ~ blockID, data=hg19_orthoBlocks_dt, FUN=function(x)sum(x=="hg19"))
lostBlocks <- newBlocks$blockID[newBlocks$genome == 0]

stopifnot(length(lostBlocks) == nrow(failed_dt))

stopifnot(lostBlocks %in% hg19_orthoBlocks_dt$blockID)

out_dt <- hg19_orthoBlocks_dt[! hg19_orthoBlocks_dt$blockID %in% lostBlocks,]

stopifnot(setequal(newBlocks$blockID[newBlocks$genome == 1], out_dt$blockID))

outFile <- file.path("inferCARs_data/Orthology.Blocks_processed_hg19.txt")
write.table(out_dt[,c("blockID", "genome", "chromo", "start", "end")], file = outFile, col.names = TRUE, row.names=FALSE, sep="\t", append=F, quote=F)
cat(paste0("... written: ", outFile, "\n"))






















