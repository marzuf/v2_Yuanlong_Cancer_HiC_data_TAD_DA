
require(foreach)
require(doMC)
registerDoMC(40)

# Rscript tad_size_desc_permut.R

plotType <- "png"
myHeight <- myWidth <- 400
plotType <- "svg"
myHeight <- myWidth <- 7

outFolder <- "TAD_SIZE_DESC_PERMUT"
dir.create(outFolder, recursive = TRUE)
  
all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
hicds = all_hicds[1]
all_hicds <- all_hicds[grep("NCI-H460", all_hicds)]

myHicds <- "ENCSR489OCU_NCI-H460_"

all_sizes_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  tad_dt <- read.delim(file.path(hicds, "genes2tad", "all_assigned_regions.txt"), header=FALSE, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
  tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
  stopifnot(!duplicated(tad_dt$region))
  tad_dt$hicds <- hicds
  tad_dt
}

save(all_sizes_dt, file ="all_sizes_dt.Rdata", version=2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_sizes_dt$size <- all_sizes_dt$end - all_sizes_dt$start + 1
all_sizes_dt$size_kb <- all_sizes_dt$size/1000
all_sizes_dt$size_log10 <- log10(all_sizes_dt$size)

all_sizes_dt$hicds_lab <- gsub(myHicds, "", all_sizes_dt$hicds)


outFile <- file.path(outFolder, paste0("allDS_tad_sizeKb_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")

plot_multiDens(
  split(all_sizes_dt$size_kb, all_sizes_dt$hicds_lab),
  plotTit = paste0(myHicds, " - obs+permut data - n=", length(all_hicds), " - TAD size [kb]"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("allDS_tad_sizeLog10_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")

plot_multiDens(
  split(all_sizes_dt$size_log10, all_sizes_dt$hicds),
  plotTit = paste0(myHicds, " - obs+permut data - n=", length(all_hicds), " - TAD size [log10]"), legPos = "topright")

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



all_g2t_dt <- foreach(hicds = all_hicds) %dopar% {
  g2t_dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), 
                       header=FALSE, col.names=c("entrezID", "chromo",  "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
  setNames(
    as.numeric(table(g2t_dt$region)),
    names(table(g2t_dt$region))
  )
}
# names(all_g2t_dt) <- all_hicds
names(all_g2t_dt) <- gsub(myHicds, "", all_hicds)

save(all_g2t_dt, file ="all_g2t_dt.Rdata", version=2)


outFile <- file.path(outFolder, paste0("allDS_nGenesByTad_nbr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot_multiDens(
  all_g2t_dt, 
  plotTit = paste0(myHicds, " - n=", nrow(nTADs_dt), " - # genes/TAD"), legPos="topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("allDS_nGenesByTad_nbr_log10_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L")
plot_multiDens(
  lapply(all_g2t_dt, log10),
  plotTit = paste0(myHicds, " - n=", nrow(nTADs_dt), " - # genes/TAD") ,legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



all_g2t_dt <- foreach(hicds = all_hicds) %dopar% {
  g2t_dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), 
                       header=FALSE, col.names=c("entrezID", "chromo",  "start", "end", "region"), stringsAsFactors = FALSE)
  
  geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, "TCGAluad_norm_luad", "0_prepGeneData", "pipeline_geneList.Rdata")))
  
  g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region) & g2t_dt$entrezID %in% geneList,]
  setNames(
    as.numeric(table(g2t_dt$region)),
    names(table(g2t_dt$region))
  )
}
# names(all_g2t_dt) <- all_hicds
names(all_g2t_dt) <- gsub(myHicds, "", all_hicds)

save(all_g2t_dt, file ="all_g2t_dt.Rdata", version=2)


outFile <- file.path(outFolder, paste0("allDS_nGenesByTad_nbr_density_pipGenes.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  all_g2t_dt, 
  plotTit = paste0(myHicds, " - n=", nrow(nTADs_dt), " - # genes/TAD"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("allDS_nGenesByTad_nbr_log10_density_pipGenes.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  lapply(all_g2t_dt, log10),
  plotTit = paste0(myHicds, " - n=", nrow(nTADs_dt), " - # genes/TAD"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


