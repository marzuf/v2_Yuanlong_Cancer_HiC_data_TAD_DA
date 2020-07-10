
# INPUT:
#   - table with signif. TADs:
#   dataset / region / 
#   
#   - coordinates of the TADs:
#   dataset / region / chromo / start / end 
# 
#  - gene-to-TAD dt [should be processed to contain only the pip genes ! ]
#     dataset / entrezID / chromo / start / end / region / symbol

# column order does not matter

# Rscript FUNC_tad_matching_signif_across_hicds_allMatch_v2_vFunc.R  

require(foreach)
require(doMC)
require(GenomicRanges)

nCpu = 40

logFile=""

### > Filter1: retain only matches with >= *minOverlapBpRatio* (80%) bp overlap
minOverlapBpRatio = 0.8


### > Filter2: retain only "conserved regions" with >= *minIntersectGenes*(3) genes at the intersect
minIntersectGenes = 3


### !!! >>> v2 update here 
### Merge conserved regions that have >= *gene_matching_fuse_threshold* % (80%) gene overlap
gene_matching_fuse_threshold = 0.8

registerDoMC(nCpu)

script0_name <- "0_prepGeneData"

outFolder <- "TMP_FUNC_CONS_vFUNC"
dir.create(outFolder)

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
entrez2symb <- setNames(gff_dt$symbol,gff_dt$entrezID)

# for testing
mainFolder <- file.path(".")
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("PERMUT", all_hicds) & !grepl("RANDOM", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

# for testing - prep the signif. tables
final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))
final_dt$dataset <- file.path(final_dt$hicds,final_dt$exprds) 
signif_dt <- final_dt[final_dt$adjPvalComb <= 0.01,c("dataset", "region")]
save(signif_dt, file=file.path(outFolder, "signif_dt.Rdata"))
load(file.path(outFolder, "signif_dt.Rdata"))

# for testing, prep the coordinates of TADs
all_tad_pos_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
  hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(hicds_file))
  tadpos_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
  tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]
  stopifnot(nrow(tadpos_dt) > 0 )
  tadpos_dt$dataset <- file.path(hicds, exprds)
  tadpos_dt <- tadpos_dt[, c("dataset", "chromo", "region", "start", "end")]
  # take only the signif ones
  signif_tads <- signif_dt$region[signif_dt$dataset == file.path(hicds, exprds)]
  tadpos_dt[tadpos_dt$region %in% signif_tads,]
  }
  exprds_dt
}
save(all_tad_pos_dt, file=file.path(outFolder, "all_tad_pos_dt.Rdata"))
load(file.path(outFolder, "all_tad_pos_dt.Rdata"))

# for testing, prep the g-2-t assignment TADs
# for testing, prep the coordinates of TADs
all_g2t_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
    stopifnot(geneList %in% g2t_dt$entrezID)
    g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
    # take only the signif ones
    signif_tads <- signif_dt$region[signif_dt$dataset == file.path(hicds, exprds)]
    g2t_dt$dataset <- file.path(hicds, exprds)
    g2t_dt <- g2t_dt[, c("dataset", "entrezID", "chromo", "start", "end", "region")]
    g2t_dt[g2t_dt$region %in% signif_tads,]
  }
  exprds_dt
}
all_g2t_dt$symbol <- entrez2symb[all_g2t_dt$entrezID]
save(all_g2t_dt, file=file.path(outFolder, "all_g2t_dt.Rdata"))
load(file.path(outFolder, "all_g2t_dt.Rdata"))


source("../MANUSCRIPT_FIGURES/COCODATA/R/conserv_func.R")
source("../MANUSCRIPT_FIGURES/COCODATA/R/utils_func.R")
conserv_minIntersectGenes <- 3
conserv_geneMatching <- 0.8
conserv_minOverlapBp <- 0.8
signif_conserv_data <- get_conservedRegion(
  signif_dt=signif_dt, all_tad_pos_dt=all_tad_pos_dt, all_g2t_dt=all_g2t_dt, 
  minOverlapBpRatio=conserv_minOverlapBp,
  minIntersectGenes=conserv_minIntersectGenes,
  gene_matching_fuse_threshold = conserv_geneMatching,
  verbose=FALSE
)

save(signif_conserv_data, file=file.path(outFolder, "signif_conserv_data.Rdata"))


