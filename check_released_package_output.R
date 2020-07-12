
#### CHECK STEP 10 - ok
vpip1_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/doc/package_emp_pval_combined.RData"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/11sameNbr_runEmpPvalCombined/emp_pval_combined.Rdata"))
stopifnot(all.equal(vpip1_dt, v0_dt))

#### CHECK STEP 4 # => ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/doc/package_meanCorr.RData"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))

#### CHECK STEP 3  # => ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/doc/package_obs_meanFC.RData"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP conserv
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/doc/package_signif_conserv_data.Rdata"))
v0_dt <- get(load("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2/conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))

x <- sort(unlist(lapply(vpip_dt[["conserved_signif_tads"]], function(x) paste0(sort(x), collapse=","))))
names(x) <- NULL
y <- sort(unlist(lapply(v0_dt, function(x) paste0(sort(x), collapse=","))))
names(y) <- NULL

stopifnot(all.equal(x,y))

#### CHECK STEP 11 - ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/doc/package_obs_ratioDown.RData"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/8cOnlyRatioDownFastSave_runAllDown/all_obs_ratioDown.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))

#### CHECK STEP 11 - ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/doc/package_permut_dt.RData"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/8cOnlyRatioDownFastSave_runAllDown/all_obs_ratioDown.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))

#### CHECK STEP 5corr
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/5sameNbr_runPermutationsCorr/sample_around_TADs_sameNbr.Rdata"))
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/COCODATA/doc/package_sampAcrossBD_data.RData"))
stopifnot(setequal(names(v0_dt) ,names(vpip_dt)))
x <- lapply(v0_dt, function(x)x[["genes"]])
y <- lapply(vpip_dt, function(x)x[["genes"]])
stopifnot(setequal(names(x) ,names(y)))
x <- x[names(y)]
stopifnot(setequal(x,y))


