
require(foreach)
require(doMC)
registerDoMC(40)

# Rscript tad_size_desc_permut_allDS.R

plotType <- "png"
myHeight <- myWidth <- 400
plotType <- "svg"
myHeight <- myWidth <- 7

outFolder <- "TAD_SIZE_DESC_PERMUT_ALLDS"
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
#all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]

all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))

# 
rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T" )

buildData <- TRUE

hicds=all_hicds[1]


all_sizes_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  
  tad_dt <- read.delim(file.path(hicds, "genes2tad", "all_assigned_regions.txt"), header=FALSE, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
  stopifnot(!duplicated(tad_dt$region))
  tad_dt$hicds <- hicds

  tad_dt
}

save(all_sizes_dt, file ="all_sizes_dt_allDS.Rdata", version=2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_sizes_dt$size <- all_sizes_dt$end - all_sizes_dt$start + 1
all_sizes_dt$size_kb <- all_sizes_dt$size/1000
all_sizes_dt$size_log10 <- log10(all_sizes_dt$size)

all_sizes_dt$hicds_lab <- gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_sizes_dt$hicds)
all_sizes_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", all_sizes_dt$hicds_lab)
all_sizes_dt$hicds_lab <- gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES", all_sizes_dt$hicds_lab)
all_sizes_dt$hicds_lab <- gsub(".+PERMUTG2T_40kb", "PERMUTG2T", all_sizes_dt$hicds_lab)
all_sizes_dt$hicds_lab[! all_sizes_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"


outFile <- file.path(outFolder, paste0("allDS_tad_sizeKb_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")

plot_multiDens(
  split(all_sizes_dt$size_kb, all_sizes_dt$hicds_lab),
  plotTit = paste0("allDS - obs+permut data - n=", length(all_hicds), " - TAD size [kb]"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("allDS_tad_sizeLog10_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")

plot_multiDens(
  split(all_sizes_dt$size_log10, all_sizes_dt$hicds_lab),
  plotTit = paste0("all DS - obs+permut data - n=", length(all_hicds), " - TAD size [log10]"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

nTADs_dt <- aggregate(region~hicds, data=all_sizes_dt, FUN=length)

# outFile <- file.path(outFolder, paste0("allDS_tad_nbr_density.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
# par(bty="L")
# plot_multiDens(
#   split(nTADs_dt$region, nTADs_dt$hicds),
#   plotTit = paste0("obs+permut data - n=", nrow(nTADs_dt), " - # TADs"))
# 
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))



all_g2t_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  g2t_dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), 
                       header=FALSE, col.names=c("entrezID", "chromo",  "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
  data.frame(
    region = names(table(g2t_dt$region)),
    nGenes = as.numeric(table(g2t_dt$region)),
    hicds = hicds,
    stringsAsFactors = FALSE) 
}
all_g2t_dt$nGenes_log10 <- log10(all_g2t_dt$nGenes)

all_g2t_dt$hicds_lab <- gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_g2t_dt$hicds)
all_g2t_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", all_g2t_dt$hicds_lab)
all_g2t_dt$hicds_lab <- gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES", all_g2t_dt$hicds_lab)
all_g2t_dt$hicds_lab <- gsub(".+PERMUTG2T_40kb", "PERMUTG2T", all_g2t_dt$hicds_lab)
all_g2t_dt$hicds_lab[! all_g2t_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"


save(all_g2t_dt, file ="all_g2t_dt.Rdata", version=2)


outFile <- file.path(outFolder, paste0("allDS_nGenesByTad_nbr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot_multiDens(
  split(all_g2t_dt$nGenes,  all_g2t_dt$hicds_lab), 
  plotTit = paste0("all DS - n=", nrow(nTADs_dt), " - # genes/TAD"), legPos="topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("allDS_nGenesByTad_nbr_log10_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot_multiDens(
  split(all_g2t_dt$nGenes_log10, all_g2t_dt$hicds_lab), 
  plotTit = paste0("all DS - n=", nrow(nTADs_dt), " - # genes/TAD") ,legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))


all_g2t_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
  g2t_dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), 
                       header=FALSE, col.names=c("entrezID", "chromo",  "start", "end", "region"), stringsAsFactors = FALSE)
  
  geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
  
  g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region) & g2t_dt$entrezID %in% geneList,]
  data.frame(
    region = names(table(g2t_dt$region)),
    nGenes = as.numeric(table(g2t_dt$region)),
    hicds = hicds,
    exprds = exprds,
    stringsAsFactors = FALSE) 
  }
}
all_g2t_dt$nGenes_log10 <- log10(all_g2t_dt$nGenes)

all_g2t_dt$hicds_lab <- gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_g2t_dt$hicds)
all_g2t_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", all_g2t_dt$hicds_lab)
all_g2t_dt$hicds_lab <- gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES", all_g2t_dt$hicds_lab)
all_g2t_dt$hicds_lab <- gsub(".+PERMUTG2T_40kb", "PERMUTG2T", all_g2t_dt$hicds_lab)
all_g2t_dt$hicds_lab[! all_g2t_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"


save(all_g2t_dt, file ="all_g2t_dt.Rdata", version=2)


outFile <- file.path(outFolder, paste0("allDS_nGenesByTad_nbr_density_pipGenes.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_g2t_dt$nGenes,  all_g2t_dt$hicds_lab), 
  plotTit = paste0("all DS - n=", nrow(nTADs_dt), " - # genes/TAD"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("allDS_nGenesByTad_nbr_log10_density_pipGenes.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_g2t_dt$nGenes_log10,  all_g2t_dt$hicds_lab), 
  plotTit = paste0("all DS - n=", nrow(nTADs_dt), " - # genes/TAD"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


