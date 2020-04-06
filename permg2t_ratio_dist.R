# Rscript permg2t_ratio_dist.R

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"
fccFile <- file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown/prodSignedRatio_permDT.Rdata")
ratioDownFile <- file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown/ratioDown_permDT.Rdata") 
ratioFcFile <- file.path("PERMG2T_TAD_FC_RATIO", hicds, exprds, "ratioFC_1000Permut_permDT.Rdata")
ratioFc_dt <- get(load(ratioFcFile))
ratioDown_dt <- get(load(ratioDownFile))
ratioDown_dt <- ratioDown_dt[,1:1000]
fcc_dt <- get(load(fccFile))
fcc_dt <- fcc_dt[,1:1000]

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- "PERMG2T_RATIO_DIST"
dir.create(outFolder)

outFile <- file.path(outFolder, "cmp_dist_ratio.png")
do.call("png", list(outFile, height=400, width=600))
plot_multiDens(
  list(
    ratioFc = as.numeric(ratioFc_dt),
    ratioDown = as.numeric(ratioDown_dt),
    FCC = as.numeric(fcc_dt)
  )
)
foo <- dev.off()
