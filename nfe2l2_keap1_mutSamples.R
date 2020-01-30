
# Rscript nfe2l2_keap1_mutSamples.R

outFolder <- "NFE2L2_KEAP1_MUTSAMPLES"
dir.create(outFolder, recursive = TRUE)

setDir = "/media/electron"
setDir = ""
gamFile <- file.path(setDir, "/mnt/ed4/marco/cancerAlterations/output/pancanAtlas32_gam_mc3_and_cnas_v5.RData")
gam_dt <- get(load(gamFile))

keap1_col <- colnames(gam_dt)[grep("KEAP1", colnames(gam_dt))]
# "DEL.KEAP1" "MUT.KEAP1"
nfe2l2_col <- colnames(gam_dt)[grep("NFE2L2", colnames(gam_dt))]
# "MUT.NFE2L2"
sub_gam_dt <- gam_dt[,c(nfe2l2_col,keap1_col)]

mut_samples <- names(rowSums(sub_gam_dt))[rowSums(sub_gam_dt) > 0]

selected_dt <- gam_dt[mut_samples, c(nfe2l2_col, keap1_col)]
stopifnot(rowSums(selected_dt) > 0)

outFile <- file.path(outFolder, "mut_samples.Rdata")
save(mut_samples, file = outFile, version=2)
cat(paste0("...written: ", outFile, "\n"))