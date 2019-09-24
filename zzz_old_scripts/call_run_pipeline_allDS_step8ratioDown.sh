#!/usr/bin/bash

######
# usual datasets
######
#>>> kidney
./run_pipeline_8ratioDown.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich  # pos1
#>>> skin
./run_pipeline_8ratioDown.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf  #  pos1
./run_pipeline_8ratioDown.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF  #  pos1
./run_pipeline_8ratioDown.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1  #  pos1

./run_pipeline_8ratioDown.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf  # pos1
./run_pipeline_8ratioDown.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF   #  pos1
./run_pipeline_8ratioDown.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1   #  pos1


# prostate
./run_pipeline_8ratioDown.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad   # neu2

./run_pipeline_8ratioDown.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad  # neu2

#./run_pipeline_8ratioDown.sh GSE73782_PC3_40kb TCGAprad_norm_prad  # 

#>>> kidney
./run_pipeline_8ratioDown.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich  # neu2

#>>> lung
./run_pipeline_8ratioDown.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR  # pos0
./run_pipeline_8ratioDown.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker  #   pos0
./run_pipeline_8ratioDown.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad  # pos0
./run_pipeline_8ratioDown.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS  #  pos0
./run_pipeline_8ratioDown.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc  # el #  pos0
./run_pipeline_8ratioDown.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  #   pos0
./run_pipeline_8ratioDown.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker  #  pos0
./run_pipeline_8ratioDown.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad  #  pos0
 # => from dataset below, added set.seed in TAD_DE_utils_fasterPermut 
./run_pipeline_8ratioDown.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS  #  pos0
./run_pipeline_8ratioDown.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc  #  pos0

#>>> breast
./run_pipeline_8ratioDown.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas  # neu2

#./run_pipeline_8ratioDown.sh GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas  # 
./run_pipeline_8ratioDown.sh Barutcu_MCF-7_40kb TCGAbrca_lum_bas # neu2

#>>> cerebellum and spinal cord

./run_pipeline_8ratioDown.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal # el 4
./run_pipeline_8ratioDown.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural # el 4
./run_pipeline_8ratioDown.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural  # el 4
./run_pipeline_8ratioDown.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc  # el 4
./run_pipeline_8ratioDown.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal  # el 4
./run_pipeline_8ratioDown.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural  # el 4
./run_pipeline_8ratioDown.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural  # el 4
./run_pipeline_8ratioDown.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc  # el 4

#>>> colon
./run_pipeline_8ratioDown.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss  # neupip1

#>>> liver
./run_pipeline_8ratioDown.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc  #  neupip1
./run_pipeline_8ratioDown.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  #  neupip1
#./run_pipeline_8ratioDown.sh GSE58752_liver_40kb TCGAlihc_norm_lihc    # update name
#./run_pipeline_8ratioDown.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1   # update name  
./run_pipeline_8ratioDown.sh LI_40kb TCGAlihc_norm_lihc  # neupip1
./run_pipeline_8ratioDown.sh LI_40kb TCGAlihc_wt_mutCTNNB1  # neupip1


#>>> leukemia
./run_pipeline_8ratioDown.sh K562_40kb TCGAlaml_wt_mutFLT3  #  neupip1

#>>> pancreas
./run_pipeline_8ratioDown.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS # neupip1


######
# additional datasets (already run 0-4)
######

#>>> lung
./run_pipeline_8ratioDown.sh LG1_40kb TCGAluad_mutKRAS_mutEGFR  # pip1
./run_pipeline_8ratioDown.sh LG1_40kb TCGAluad_nonsmoker_smoker  # pip1
./run_pipeline_8ratioDown.sh LG1_40kb TCGAluad_norm_luad    # pip1
./run_pipeline_8ratioDown.sh LG1_40kb TCGAluad_wt_mutKRAS    # pip1
./run_pipeline_8ratioDown.sh LG1_40kb TCGAlusc_norm_lusc    # pip1

./run_pipeline_8ratioDown.sh LG2_40kb TCGAluad_mutKRAS_mutEGFR  # pip1
./run_pipeline_8ratioDown.sh LG2_40kb TCGAluad_nonsmoker_smoker   # pip1
./run_pipeline_8ratioDown.sh LG2_40kb TCGAluad_norm_luad   # pip1
./run_pipeline_8ratioDown.sh LG2_40kb TCGAluad_wt_mutKRAS    # pip1
./run_pipeline_8ratioDown.sh LG2_40kb TCGAlusc_norm_lusc   # pip1
 
#>>> pancreas
./run_pipeline_8ratioDown.sh PA2_40kb TCGApaad_wt_mutKRAS # pip1

./run_pipeline_8ratioDown.sh PA3_40kb TCGApaad_wt_mutKRAS # pip1

./run_pipeline_8ratioDown.sh GSE118588_Panc_beta_40kb TCGApaad_wt_mutKRAS # pip1


#>>> breast
./run_pipeline_8ratioDown.sh Barutcu_MCF-10A_40kb TCGAbrca_lum_bas #  pip2
./run_pipeline_8ratioDown.sh HMEC_40kb TCGAbrca_lum_bas  #  pip2
./run_pipeline_8ratioDown.sh GSE109229_BT474_40kb TCGAbrca_lum_bas  #   pip2
./run_pipeline_8ratioDown.sh GSE109229_SKBR3_40kb TCGAbrca_lum_bas  #  pip2



#>>> colon
./run_pipeline_8ratioDown.sh Rao_HCT-116_2017_40kb TCGAcoad_msi_mss # pip2

./run_pipeline_8ratioDown.sh ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss # pip2


#>>> kidney
./run_pipeline_8ratioDown.sh GSE99051_786_O_40kb TCGAkich_norm_kich # pip2

#>>> prostate 
./run_pipeline_8ratioDown.sh GSE118514_22Rv1_40kb TCGAprad_norm_prad # pip2




