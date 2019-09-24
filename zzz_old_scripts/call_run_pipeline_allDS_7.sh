#!/usr/bin/bash

######
# usual datasets
######
#>>> kidney
./run_pipeline_7.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich  # ok

#>>> skin
./run_pipeline_7.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf  # ok
./run_pipeline_7.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF  # ok
./run_pipeline_7.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1  #ok

./run_pipeline_7.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf  # pip1
./run_pipeline_7.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF   #- pip1
./run_pipeline_7.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1   #-pip1


# prostate
./run_pipeline_7.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad   # pos

./run_pipeline_7.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad  # pos

#./run_pipeline_7.sh GSE73782_PC3_40kb TCGAprad_norm_prad  # discard low quality

#>>> kidney
./run_pipeline_7.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich  # pos

#>>> lung
./run_pipeline_7.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR  # pos
./run_pipeline_7.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker  #  pos
./run_pipeline_7.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad   # pos
./run_pipeline_7.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS  # pos
./run_pipeline_7.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc  # #pos
./run_pipeline_7.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  # el
./run_pipeline_7.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker  #el
./run_pipeline_7.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad #el
 # => from dataset below, added set.seed in TAD_DE_utils_fasterPermut
./run_pipeline_7.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS  # el
./run_pipeline_7.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc  #el

#>>> breast
./run_pipeline_7.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas  # pip1

#./run_pipeline_7.sh GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas  # update name
./run_pipeline_7.sh Barutcu_MCF-7_40kb TCGAbrca_lum_bas # pip1

#>>> cerebellum and spinal cord

./run_pipeline_7.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal  ## pip1
./run_pipeline_7.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural # # pip1
./run_pipeline_7.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural  ## pip1
./run_pipeline_7.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc  # # pip1
./run_pipeline_7.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal  # # pip1
./run_pipeline_7.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural  # # pip1
./run_pipeline_7.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural  # # pip1
./run_pipeline_7.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc  # # pip1

#>>> colon
./run_pipeline_7.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss # pip2

#>>> liver
./run_pipeline_7.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc  # pip2-
./run_pipeline_7.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  # pip2-
#./run_pipeline_7.sh GSE58752_liver_40kb TCGAlihc_norm_lihc    # update name
#./run_pipeline_7.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1   # update name  
./run_pipeline_7.sh LI_40kb TCGAlihc_norm_lihc  ## pip2
./run_pipeline_7.sh LI_40kb TCGAlihc_wt_mutCTNNB1  ## pip2


#>>> leukemia
./run_pipeline_7.sh K562_40kb TCGAlaml_wt_mutFLT3  # pip2

#>>> pancreas
./run_pipeline_7.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS ## pip2


######
# additional datasets (already run 0-4)
######

#>>> lung
./run_pipeline_7.sh LG1_40kb TCGAluad_mutKRAS_mutEGFR # el 
./run_pipeline_7.sh LG1_40kb TCGAluad_nonsmoker_smoker  # el
./run_pipeline_7.sh LG1_40kb TCGAluad_norm_luad    # el 
./run_pipeline_7.sh LG1_40kb TCGAluad_wt_mutKRAS    # el 
./run_pipeline_7.sh LG1_40kb TCGAlusc_norm_lusc    # el

./run_pipeline_7.sh LG2_40kb TCGAluad_mutKRAS_mutEGFR # pos1
./run_pipeline_7.sh LG2_40kb TCGAluad_nonsmoker_smoker   # pos1
./run_pipeline_7.sh LG2_40kb TCGAluad_norm_luad   # pos1
./run_pipeline_7.sh LG2_40kb TCGAluad_wt_mutKRAS    # pos1
./run_pipeline_7.sh LG2_40kb TCGAlusc_norm_lusc   # pos1
 
#>>> pancreas
./run_pipeline_7.sh PA2_40kb TCGApaad_wt_mutKRAS # pos2
./run_pipeline_7.sh PA3_40kb TCGApaad_wt_mutKRAS # pos2
./run_pipeline_7.sh GSE118588_Panc_beta_40kb TCGApaad_wt_mutKRAS # pos2


#>>> breast
./run_pipeline_7.sh Barutcu_MCF-10A_40kb TCGAbrca_lum_bas # el2
./run_pipeline_7.sh HMEC_40kb TCGAbrca_lum_bas  # el2
./run_pipeline_7.sh GSE109229_BT474_40kb TCGAbrca_lum_bas  #  el2
./run_pipeline_7.sh GSE109229_SKBR3_40kb TCGAbrca_lum_bas  # el2



#>>> colon
./run_pipeline_7.sh Rao_HCT-116_2017_40kb TCGAcoad_msi_mss # ops1

./run_pipeline_7.sh ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss # pos2


#>>> kidney
./run_pipeline_7.sh GSE99051_786_O_40kb TCGAkich_norm_kich # el1

#>>> prostate 
./run_pipeline_7.sh GSE118514_22Rv1_40kb TCGAprad_norm_prad # el2




