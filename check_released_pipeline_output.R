#### CHECK STEP 0/1  # => ok

vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000//1_prepGeneData/pipeline_geneList.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/0_prepGeneData/pipeline_geneList.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 1/2  # => ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000//2_runGeneDE/DE_topTable.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/1_runGeneDE/DE_topTable.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))

#### CHECK STEP 3  # => ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000//3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 4 # => ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000//4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 5fc  # => ok-not testable
# # (check for 1000 on electron)
# vpip1_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/5fc_runPermutationsMedian/1000_v1_permutationsDT.Rdata"))
# vpip2_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/5fc_runPermutationsMedian/1000_v2_permutationsDT.Rdata"))
# v0_dt <- get(load("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/5fc_runPermutationsMedian/permutationsDT.Rdata"))                                                                                                
# v0_cut <- v0_dt[,1:1000]
# stopifnot(rownames(vpip1_dt) == rownames(v0_dt)) 
# # cannot check the g2t assignment -> was not run with seed (not sure)
# stopifnot(all.equal(vpip1_dt, vpip2_dt))f -> ok, reproducible
# was not reproducible -> so copy paste the permutDT to check steps 6 and 8

#### CHECK STEP 5corr # => ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000/5corr_runPermutationsCorr/sample_around_TADs_sameNbr.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/5sameNbr_runPermutationsCorr/sample_around_TADs_sameNbr.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 6
# not testable for FULL100000 
#
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000/5corr_runPermutationsCorr/sample_around_TADs_sameNbr.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/5sameNbr_runPermutationsCorr/sample_around_TADs_sameNbr.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 7 # => ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000/7_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/7sameNbr_runPermutationsMeanTADCorr/meanCorr_sample_around_TADs_sameNbr.Rdata"))
stopifnot(all.equal(vpip_dt, v0_dt))


#### CHECK STEP 8 # => ok
vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000/8_runEmpPvalMeanTADLogFC/emp_pval_meanLogFC.Rdata"))
vpip_dt <- sort(vpip_dt)
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/9_runEmpPvalMeanTADLogFC/emp_pval_meanLogFC.Rdata"))
v0_dt <- sort(v0_dt)
stopifnot(setequal(names(vpip_dt), names(v0_dt)))
v0_dt <- v0_dt[names(vpip_dt)]
cor(v0_dt, vpip_dt) # 0.999
cor(v0_dt, vpip_dt, method="spearman") # 0.999
v0_dt <- sort(v0_dt)
vpip_dt <- sort(vpip_dt)
length(intersect(names(vpip_dt)[1:10], names(v0_dt)[1:10])) # -> 10 not the same but they are the same top 10
names(vpip_dt)[1:10] == names(v0_dt)[1:10]
length(intersect(names(vpip_dt)[1:100], names(v0_dt)[1:100])) # -> 99


#### CHECK STEP 9  # => ok
vpip1_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000/9_runEmpPvalMeanTADCorr/fromFile_emp_pval_meanCorr.Rdata"))
vpip2_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000/9_runEmpPvalMeanTADCorr/fromFolder_emp_pval_meanCorr.Rdata"))
stopifnot(all.equal(vpip2_dt, vpip1_dt))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/10sameNbr_runEmpPvalMeanTADCorr/emp_pval_meanCorr.Rdata"))
stopifnot(all.equal(vpip1_dt, v0_dt))


#### CHECK STEP 10
vpip1_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000/10_runEmpPvalCombined/fromFilePermCorr_emp_pval_combined.Rdata"))
vpip1a_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000/10_runEmpPvalCombined/fromFilePermCorr_adj_emp_pval_combined.Rdata"))
vpip2_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad_FULL100000/10_runEmpPvalCombined/fromFilePermCorr_emp_pval_combined.Rdata"))
stopifnot(all.equal(vpip2_dt, vpip1_dt))
stopifnot(p.adjust(vpip1_dt, method="BH") == vpip1a_dt)
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/11sameNbr_runEmpPvalCombined/emp_pval_combined.Rdata"))
stopifnot(setequal(names(vpip1_dt), names(v0_dt)))
stopifnot(setequal(names(vpip2_dt), names(v0_dt)))

vpip1_dt <- vpip1_dt[names(v0_dt)]
cor(vpip1_dt, v0_dt) # 0.999 pearson, 0.999 spearman

vpip1_dt <- sort(vpip1_dt)
v0_dt <- sort(v0_dt)
length(intersect(names(vpip1_dt)[1:10], names(v0_dt)[1:10])) # -> 9 
names(vpip1_dt)[1:10] == names(v0_dt)[1:10]
length(intersect(names(vpip1_dt)[1:100], names(v0_dt)[1:100])) # -> 99

vpip1_dt <-p.adjust(vpip1_dt, method="BH")
v0_dt <-p.adjust(v0_dt, method="BH")
vpip1_dt <- vpip1_dt[names(v0_dt)]
cor(vpip1_dt, v0_dt) # 0.999 pearson, 0.999 spearman








vpip_dt <- get(load("../MANUSCRIPT_FIGURES/code/EXAMPLE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))
v0_dt <- get(load("PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/4_runMeanTADCorr/all_meanCorr_TAD.Rdata"))




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