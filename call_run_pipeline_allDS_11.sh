#!/usr/bin/bash

######
# usual datasets
######
#>>> kidney
./run_pipeline_11.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich  # el

#>>> skin
./run_pipeline_11.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf  # pip1
./run_pipeline_11.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF  # pip1
./run_pipeline_11.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1  #pip1

./run_pipeline_11.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf  #pip1
./run_pipeline_11.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF   #pip1
./run_pipeline_11.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1   #pip1


# prostate
./run_pipeline_11.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad   #pip1

./run_pipeline_11.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad  #pip1

#./run_pipeline_11.sh GSE73782_PC3_40kb TCGAprad_norm_prad  # discard low quality

#>>> kidney
./run_pipeline_11.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich  # pip1

#>>> lung
./run_pipeline_11.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR # pip2
./run_pipeline_11.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker  # pip2
./run_pipeline_11.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad   # pip2
./run_pipeline_11.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS  # pip2
./run_pipeline_11.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc  # ## pip2
./run_pipeline_11.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  # pip2
./run_pipeline_11.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker # pip2
./run_pipeline_11.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad # pip2
 # => from dataset below, added set.seed in TAD_DE_utils_fasterPermut
./run_pipeline_11.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS  # pip2
./run_pipeline_11.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc # pip2

#>>> breast
./run_pipeline_11.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas  #pos1

#./run_pipeline_11.sh GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas  # update name
./run_pipeline_11.sh Barutcu_MCF-7_40kb TCGAbrca_lum_bas #pos1

#>>> cerebellum and spinal cord

./run_pipeline_11.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal  ##pos1
./run_pipeline_11.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural # #pos1
./run_pipeline_11.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural  ##pos1
./run_pipeline_11.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc  # #pos1
./run_pipeline_11.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal  # #pos1
./run_pipeline_11.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural  # #pos1
./run_pipeline_11.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural  # #pos1
./run_pipeline_11.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc  # #pos1

#>>> colon
./run_pipeline_11.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss # #pos1

#>>> liver
./run_pipeline_11.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc  # pos2
./run_pipeline_11.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  # # pos2
#./run_pipeline_11.sh GSE58752_liver_40kb TCGAlihc_norm_lihc    # update name
#./run_pipeline_11.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1   # update name  
./run_pipeline_11.sh LI_40kb TCGAlihc_norm_lihc  ## # pos2
./run_pipeline_11.sh LI_40kb TCGAlihc_wt_mutCTNNB1  ## # pos2


#>>> leukemia
./run_pipeline_11.sh K562_40kb TCGAlaml_wt_mutFLT3  # # pos2

#>>> pancreas
./run_pipeline_11.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS ### pos2 


######
# additional datasets (already run 0-4)
######

#>>> lung
./run_pipeline_11.sh LG1_40kb TCGAluad_mutKRAS_mutEGFR  # pos2
./run_pipeline_11.sh LG1_40kb TCGAluad_nonsmoker_smoker  # pos2
./run_pipeline_11.sh LG1_40kb TCGAluad_norm_luad     # pos2
./run_pipeline_11.sh LG1_40kb TCGAluad_wt_mutKRAS     # pos2
./run_pipeline_11.sh LG1_40kb TCGAlusc_norm_lusc    # pos2

./run_pipeline_11.sh LG2_40kb TCGAluad_mutKRAS_mutEGFR # el1
./run_pipeline_11.sh LG2_40kb TCGAluad_nonsmoker_smoker   # el1
./run_pipeline_11.sh LG2_40kb TCGAluad_norm_luad   # el1
./run_pipeline_11.sh LG2_40kb TCGAluad_wt_mutKRAS # el1
./run_pipeline_11.sh LG2_40kb TCGAlusc_norm_lusc   # el1
 
#>>> pancreas
./run_pipeline_11.sh PA2_40kb TCGApaad_wt_mutKRAS # el1
./run_pipeline_11.sh PA3_40kb TCGApaad_wt_mutKRAS # el1
./run_pipeline_11.sh GSE118588_Panc_beta_40kb TCGApaad_wt_mutKRAS # el1


#>>> breast
./run_pipeline_11.sh Barutcu_MCF-10A_40kb TCGAbrca_lum_bas # el1
./run_pipeline_11.sh HMEC_40kb TCGAbrca_lum_bas  # el2
./run_pipeline_11.sh GSE109229_BT474_40kb TCGAbrca_lum_bas  ## el2
./run_pipeline_11.sh GSE109229_SKBR3_40kb TCGAbrca_lum_bas  # el2



#>>> colon
./run_pipeline_11.sh Rao_HCT-116_2017_40kb TCGAcoad_msi_mss # # el2

./run_pipeline_11.sh ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss # el2


#>>> kidney
./run_pipeline_11.sh GSE99051_786_O_40kb TCGAkich_norm_kich # el2

#>>> prostate 
./run_pipeline_11.sh GSE118514_22Rv1_40kb TCGAprad_norm_prad # el2




