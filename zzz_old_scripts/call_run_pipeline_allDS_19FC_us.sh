#!/usr/bin/bash

# => run all on electron
######
# usual datasets
######
#>>> kidney
./run_pipeline_19fc.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich  # 5,6,7 - el ok

#>>> skin
./run_pipeline_19fc.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf  # 5,6,7 - el ok
./run_pipeline_19fc.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF  # 5,6,7 - pos ok
./run_pipeline_19fc.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1  # 5,6,7 - pos ok

./run_pipeline_56.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf  # neutrino - pip1 5,6 ok
./run_pipeline_56.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF   # neutrino - pip2 5,6 ok
./run_pipeline_19fc.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1   # neutrino - ok


# prostate
./run_pipeline_19fc.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad   # neutrino - pip1 => ok

./run_pipeline_19fc.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad  # el => ok

#./run_pipeline_19fc.sh GSE73782_PC3_40kb TCGAprad_norm_prad  # discard low quality

#>>> kidney
./run_pipeline_19fc.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich  # pos - ok

#>>> lung
./run_pipeline_19fc.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR  # pos - ok
./run_pipeline_19fc.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker  #  pos ok
./run_pipeline_19fc.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad  # pos   - ok                 # => from datasets below, run with 5fastSavePermut and 6fastSave
./run_pipeline_19fc.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS  # neutrino pip1
./run_pipeline_19fc.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc  # el # ok
./run_pipeline_19fc.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  # neutrino pip2 ok
./run_pipeline_19fc.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker  # pos # ok
./run_pipeline_19fc.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad  # pos ok
 # => from dataset below, added set.seed in TAD_DE_utils_fasterPermut
./run_pipeline_19fc.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS  # neutrino pip 2 - ok
./run_pipeline_19fc.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc  # neutrino pip 1- ok

#>>> breast
./run_pipeline_19fc.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas  # el ok

#./run_pipeline_19fc.sh GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas  # update name
./run_pipeline_19fc.sh Barutcu_MCF-7_40kb TCGAbrca_lum_bas # pos ok

#>>> cerebellum and spinal cord

./run_pipeline_19fc.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal  # => run multi on electron! (savePermut but not all (?) set.seed) only the 1st without set.seed ?! ok
./run_pipeline_19fc.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural # => run multi on electron!(savePermut but not all (?) set.seed) only the 1st without set.seed ?! ok
./run_pipeline_19fc.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural  # => run multi on electron! (savePermut but not all (?) set.seed) only the 1st without set.seed ?! ok
./run_pipeline_19fc.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc  # => run multi on electron!(savePermut but not all (?) set.seed) only the 1st without set.seed ?! ok
./run_pipeline_19fc.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal  # => run multi on electron!(savePermut but not all (?) set.seed) only the 1st without set.seed ?! ok
./run_pipeline_19fc.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural  # => run multi on electron!(savePermut but not all (?) set.seed) only the 1st without set.seed ?! ok
./run_pipeline_19fc.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural  # => run multi on electron!(savePermut but not all (?) set.seed) only the 1st without set.seed ?! ok
./run_pipeline_19fc.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc  # => run multi on electron!(savePermut but not all (?) set.seed) only the 1st without set.seed ?! ok

#>>> colon
./run_pipeline_19fc.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss  # pos ok

#>>> liver
./run_pipeline_19fc.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc  # pos - ok
./run_pipeline_19fc.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  # pos - ok
#./run_pipeline_19fc.sh GSE58752_liver_40kb TCGAlihc_norm_lihc    # update name
#./run_pipeline_19fc.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1   # update name  
./run_pipeline_19fc.sh LI_40kb TCGAlihc_norm_lihc  # el ok
./run_pipeline_19fc.sh LI_40kb TCGAlihc_wt_mutCTNNB1  # el ok


#>>> leukemia
./run_pipeline_19fc.sh K562_40kb TCGAlaml_wt_mutFLT3  # neutrino pip2 ok

#>>> pancreas
./run_pipeline_19fc.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS # neu pip 1 ok



