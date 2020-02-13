
require(foreach)
require(doMC)
registerDoMC(40)

# Rscript tad_size_desc.R

plotType <- "png"
myHeight <- myWidth <- 400

outFolder <- "TAD_SIZE_DESC"
dir.create(outFolder, recursive = TRUE)
  
all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
hicds = all_hicds[1]

all_sizes_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  tad_dt <- read.delim(file.path(hicds, "genes2tad", "all_assigned_regions.txt"), header=FALSE, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
  stopifnot(!duplicated(tad_dt$region))
  tad_dt$hicds <- hicds
  tad_dt
}

all_sizes_dt$size <- all_sizes_dt$end - all_sizes_dt$start + 1
all_sizes_dt$size_kb <- all_sizes_dt$size/1000

outFile <- file.path(outFolder, paste0("allDS_tad_size_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot(density(
  all_sizes_dt$size_kb
), main = paste0("all hicds - n=", length(all_hicds), " - TAD size [kb]"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

nTADs_dt <- aggregate(region~hicds, data=all_sizes_dt, FUN=length)

outFile <- file.path(outFolder, paste0("allDS_tad_nbr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot(density(
  nTADs_dt$region
), main = paste0("all hicds - n=", nrow(nTADs_dt), " - # TADs"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
