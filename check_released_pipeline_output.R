#### CHECK STEP 0/1

vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/1_prepGeneData/pipeline_geneList.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/0_prepGeneData/pipeline_geneList.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 1/2
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/2_runGeneDE/DE_topTable.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/1_runGeneDE/DE_topTable.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))

#### CHECK STEP 3
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 4
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 5fc
# (check for 1000 on electron)
vpip1_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/5fc_runPermutationsMedian/1000_v1_permutationsDT.Rdata"))
vpip2_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/5fc_runPermutationsMedian/1000_v2_permutationsDT.Rdata"))
v0_dt <- get(load("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/5fc_runPermutationsMedian/permutationsDT.Rdata"))                                                                                                
v0_cut <- v0_dt[,1:1000]
stopifnot(rownames(vpip1_dt) == rownames(v0_dt)) 
# cannot check the g2t assignment -> was not run with seed (not sure)
stopifnot(all.equal(vpip1_dt, vpip2_dt))f -> ok, reproducible

#### CHECK STEP 5corr
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))

#### CHECK STEP 6
EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/7_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))

#### CHECK STEP 7
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/7_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/7sameNbr_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 8

EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/7_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata

vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))

#### CHECK STEP 9
EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/7_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 10

vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))



#### CHECK STEP 11
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/11_runFCC/all_FCC_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))



#### CHECK STEP 12
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/11_runFCC/all_FCC_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))



#### CHECK STEP 13
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/11_runFCC/all_FCC_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))