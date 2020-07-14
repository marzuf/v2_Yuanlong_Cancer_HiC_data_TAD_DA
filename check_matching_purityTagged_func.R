
# Rscript check_matching_purityTagged_func.R


source("conserv_func.R")
source("../MANUSCRIPT_FIGURES/COCODATA/R/utils_func.R")

options(scipen=100)


# Rscript check_matching_purityFilter_func.R

require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
  purity_plot_name <- "EPIC"
} else{
  purity_ds <- ""
  purity_plot_name <- "aran"
}

### HARD-CODED - MAIN SETTINGS
corMet <- "pearson"
transfExpr <- "log10"
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)

purity_file <- file.path("ALLTADS_AND_PURITY", purity_ds, transfExpr, "all_ds_corrPurity_dt.Rdata")
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

result_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

merge_dt$region_id <- file.path(merge_dt$dataset, merge_dt$region)
taggedTADs <- merge_dt$region_id[merge_dt$purityCorr <=  purityCorrThresh]

# for testing - prep the signif. tables
final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))
final_dt$dataset <- file.path(final_dt$hicds,final_dt$exprds) 
final_dt$regID <- file.path(final_dt$dataset, final_dt$region)
stopifnot(taggedTADs %in% final_dt$regID)
tagged_dt <- final_dt[final_dt$regID %in% taggedTADs,]

# for testing
mainFolder <- file.path(".")
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("PERMUT", all_hicds) & !grepl("RANDOM", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

script0_name <- "0_prepGeneData"

# for testing, prep the coordinates of TADs
all_tad_pos_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(hicds_file))
    tadpos_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
    tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]
    stopifnot(nrow(tadpos_dt) > 0 )
    tadpos_dt$dataset <- file.path(hicds, exprds)
    tadpos_dt$regID <- file.path(tadpos_dt$dataset, tadpos_dt$region)
    tadpos_dt <- tadpos_dt[, c("dataset", "chromo", "region", "start", "end", "regID")]
    # take only the tagged ones
    tadpos_dt[tadpos_dt$regID %in% taggedTADs,]
  }
  exprds_dt
}

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
    # take only the tagged ones
    g2t_dt$dataset <- file.path(hicds, exprds)
    g2t_dt$regID <- file.path(g2t_dt$dataset, g2t_dt$region)
    g2t_dt[g2t_dt$regID %in% taggedTADs,]
  }
  exprds_dt
}

all_conserv_data <- get_conservedRegion(signif_dt = tagged_dt, 
                                        all_g2t_dt = all_g2t_dt, 
                                        all_tad_pos_dt = all_tad_pos_dt)

save(all_conserv_data, file="all_conserv_data.Rdata", version=2)



vFunc_conserved_tads <- all_conserv_data[[paste0("conserved_signif_tads")]]
vScript_conserved_tads <- get(load(file.path("TAD_MATCHING_PURITYTAGGED_ACROSS_HICDS_ALLMATCH_v2", purity_plot_name, transfExpr, "conserved_purityTagged_tads_qt0.05_minBpRatio0.8_minInterGenes3.Rdata")))

stopifnot(length(vFunc_conserved_tads) == length(vScript_conserved_tads))

names(vScript_conserved_tads) <- names(vFunc_conserved_tads) <- NULL

vFunc <- sort(unlist(lapply(vFunc_conserved_tads, function(x) paste0(sort(x), collapse=","))))
vScript <- sort(unlist(lapply(vScript_conserved_tads, function(x) paste0(sort(x), collapse=","))))
stopifnot(all.equal(vFunc, vScript))

stop("ok\n")



















