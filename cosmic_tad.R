cancer_genes_file <- "cancer_gene_census.csv"
cg_dt <- read.delim(cancer_genes_file, header = TRUE, stringsAsFactors = FALSE, sep=",")
nrow(cg_dt)
onco_dt <- cg_dt[grepl("oncogene", cg_dt$Role.in.Cancer),]
nrow(onco_dt)
tsg_dt <- cg_dt[grepl("TSG", cg_dt$Role.in.Cancer),]
nrow(tsg_dt)


# setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
colnames(gff_dt)[colnames(gff_dt) == "symbol"] <- "Gene.Symbol"
onco_entrez_dt <- merge(onco_dt[, "Gene.Symbol", drop=F], gff_dt[, c("Gene.Symbol", "entrezID")], by="Gene.Symbol", all.x=FALSE, all.y=FALSE)
nrow(onco_entrez_dt)
tsg_entrez_dt <- merge(tsg_dt[, "Gene.Symbol", drop=F], gff_dt[, c("Gene.Symbol", "entrezID")], by="Gene.Symbol", all.x=FALSE, all.y=FALSE)
nrow(tsg_entrez_dt)

require(foreach)
require(doMC)
registerDoMC(40)

# Rscript cosmic_tad.R

outFolder <- "COSMIC_TAD"
dir.create(outFolder)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <-  "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4

qt_dt<- read.delim("hicds_sparsity.csv", sep="\t", header = TRUE)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
hicds = all_hicds[1]
all_hicds <-all_hicds[ grepl("H460", all_hicds)]

all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  g2t_dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"), header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
  g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
  g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
  
  tot_dt <- data.frame(region = as.character(names(table(g2t_dt$region))),
                       nbrGenes = as.numeric(table(g2t_dt$region)),
                       stringsAsFactors = FALSE
                         )
  onco_agg_dt <- aggregate(entrezID ~ region, data=g2t_dt, function(x) sum(x %in% onco_entrez_dt$entrezID))
  colnames(onco_agg_dt)[colnames(onco_agg_dt) == "entrezID"] <- "nbrOncogenes"
  tsg_agg_dt <- aggregate(entrezID ~ region, data=g2t_dt, function(x) sum(x %in% tsg_entrez_dt$entrezID))
  colnames(tsg_agg_dt)[colnames(tsg_agg_dt) == "entrezID"] <- "nbrTSG"

  merged_dt <- merge(merge(tot_dt, onco_agg_dt, by="region", all=TRUE), tsg_agg_dt, by="region", all=TRUE)
  merged_dt$hicds <- hicds
  merged_dt
}

all_dt$hicds_lab <- gsub(".+_(.+?_40kb)","\\1",  all_dt$hicds)
all_dt$hicds_lab <- ifelse(grepl("RANDOM",all_dt$hicds_lab) | grepl("PERMUT", all_dt$hicds_lab),all_dt$hicds_lab,  "OBSERVED" )

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFile <- file.path(outFolder, paste0("nbrOncogenes_densityplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  split(all_dt$nbrOncogenes, all_dt$hicds_lab),
  plotTit = paste0("# oncogenes")
)
foo <- dev.off()



outFile <- file.path(outFolder, paste0("nbrOncogenes_noZero_densityplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  split(all_dt$nbrOncogenes[all_dt$nbrOncogenes > 0], all_dt$hicds_lab[all_dt$nbrOncogenes > 0]),
  plotTit = paste0("# oncogenes/TAD (no 0)")
)
foo <- dev.off()

outFile <- file.path(outFolder, paste0("nbrOncogenes_densityplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  split(all_dt$nbrOncogenes, all_dt$hicds_lab),
  plotTit = paste0("# oncogenes/TAD")
)
foo <- dev.off()


outFile <- file.path(outFolder, paste0("nbrTSG_densityplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  split(all_dt$nbrTSG, all_dt$hicds_lab),
  plotTit = paste0("# TSGs/TAD")
)
foo <- dev.off()

outFile <- file.path(outFolder, paste0("nbrTSG_noZero_densityplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  split(all_dt$nbrTSG[all_dt$nbrTSG > 0], all_dt$hicds_lab[all_dt$nbrTSG > 0]),
  plotTit = paste0("# TSGs/TAD (no 0)")
)
foo <- dev.off()

outFile <- file.path(outFolder, paste0("nbrTSG_noZero_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(nbrTSG~hicds_lab, data = all_dt[all_dt$nbrTSG > 0,], main = "# TSG/TAD (no 0)")
foo <- dev.off()

outFile <- file.path(outFolder, paste0("nbrOncogenes_noZero_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(nbrOncogenes~hicds_lab, data = all_dt[all_dt$nbrOncogenes > 0,], main = "# oncogenes/TAD (no 0)")
foo <- dev.off()




