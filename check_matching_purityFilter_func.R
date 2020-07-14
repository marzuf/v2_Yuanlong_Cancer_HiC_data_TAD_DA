
# Rscript check_matching_purityFilter_func.R

source("../MANUSCRIPT_FIGURES/COCODATA/R/conserv_func.R")
source("../MANUSCRIPT_FIGURES/COCODATA/R/utils_func.R")

options(scipen=100)

set.seed(20052020)

allSignif_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/data/allSignif_dt.RData"))
allSignif_dt$regID <- file.path(allSignif_dt$dataset, allSignif_dt$region)
allSignif_g2t_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/data/allSignif_g2t_dt.RData"))
allSignif_g2t_dt$regID <- file.path(allSignif_g2t_dt$dataset, allSignif_g2t_dt$region)
allSignif_tad_pos_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/data/allSignif_tad_pos_dt.RData"))
allSignif_tad_pos_dt$regID <- file.path(allSignif_tad_pos_dt$dataset, allSignif_tad_pos_dt$region)
# > head(allSignif_dt)
# dataset       region
# 1514 Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr6_TAD633
# 121  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr1_TAD551
# 932  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas chr19_TAD244
# > head(allSignif_g2t_dt)
# dataset entrezID chromo    start      end       region  symbol
# 147  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas    51182  chr10 14880159 14913740  chr10_TAD67  HSPA14
# 149  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas    79723  chr10 14920782 14946314  chr10_TAD67 SUV39H2
# > head(allSignif_tad_pos_dt)
# dataset chromo       region     start       end
# 68   Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr10  chr10_TAD67  14880001  15040000
# 821  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr11 chr11_TAD238  61920001  62160000

# Rscript check_matching_purityFilter_func.R

require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)



### HARD-CODED - MAIN SETTINGS
# purity_ds <- "EPIC"
# purity_plot_name <- "EPIC"
purity_ds <- ""
purity_plot_name <- "aran"
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
signifTADs_to_discard <- merge_dt$region_id[merge_dt$purityCorr <=  purityCorrThresh & merge_dt$signif  ]

stopifnot(signifTADs_to_discard %in% allSignif_dt$regID)
allSignif_dt <- allSignif_dt[! allSignif_dt$regID %in% signifTADs_to_discard,]
stopifnot(signifTADs_to_discard %in% allSignif_tad_pos_dt$regID)
allSignif_tad_pos_dt <- allSignif_tad_pos_dt[! allSignif_tad_pos_dt$regID %in% signifTADs_to_discard,]
stopifnot(signifTADs_to_discard %in% allSignif_g2t_dt$regID)
allSignif_g2t_dt <- allSignif_g2t_dt[! allSignif_g2t_dt$regID %in% signifTADs_to_discard,]

all_conserv_data <- get_conservedRegion(signif_dt = allSignif_dt, all_g2t_dt = allSignif_g2t_dt, all_tad_pos_dt = allSignif_tad_pos_dt)

vFunc_conserved_tads <- all_conserv_data[[paste0("conserved_signif_tads")]]
vScript_conserved_tads <- get(load(file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", purity_ds, transfExpr, "conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))

stopifnot(length(vFunc_conserved_tads) == length(vScript_conserved_tads))

names(vScript_conserved_tads) <- names(vFunc_conserved_tads) <- NULL

vFunc <- sort(unlist(lapply(vFunc_conserved_tads, function(x) paste0(sort(x), collapse=","))))
vScript <- sort(unlist(lapply(vScript_conserved_tads, function(x) paste0(sort(x), collapse=","))))
stopifnot(all.equal(vFunc, vScript))

stop("ok\n")



















