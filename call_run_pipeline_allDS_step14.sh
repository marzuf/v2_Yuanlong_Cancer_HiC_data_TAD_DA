#!/usr/bin/bash


######
# usual datasets #  ALL !
######
#>>> kidney
#./run_pipeline_14f2Folders.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich  # el

#>>> skin
./run_pipeline_14f2Folders.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf  #  el
./run_pipeline_14f2Folders.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF  #  el
./run_pipeline_14f2Folders.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1  #  el

./run_pipeline_14f2Folders.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf  # el
./run_pipeline_14f2Folders.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF   #  el
./run_pipeline_14f2Folders.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1   #  el


# prostate
./run_pipeline_14f2Folders.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad   #el

./run_pipeline_14f2Folders.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad  # el

#./run_pipeline_14f2Folders.sh GSE73782_PC3_40kb TCGAprad_norm_prad  # 

#>>> kidney
./run_pipeline_14f2Folders.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich  # el

#>>> lung
./run_pipeline_14f2Folders.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR  # el2
./run_pipeline_14f2Folders.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker  #   el2
./run_pipeline_14f2Folders.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad  # el2
./run_pipeline_14f2Folders.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS  #  el2
./run_pipeline_14f2Folders.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc   #  el2
./run_pipeline_14f2Folders.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  #el2   
./run_pipeline_14f2Folders.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker  #el2  
./run_pipeline_14f2Folders.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad  #  el2
 # => from dataset below, added set.seed in TAD_DE_utils_fasterPermut 
./run_pipeline_14f2Folders.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS  #  el2
./run_pipeline_14f2Folders.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc  #  el2

#>>> breast
./run_pipeline_14f2Folders.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas  # 

#./run_pipeline_14f2Folders.sh GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas  # 
./run_pipeline_14f2Folders.sh Barutcu_MCF-7_40kb TCGAbrca_lum_bas # 

#>>> cerebellum and spinal cord

./run_pipeline_14f2Folders.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal  
./run_pipeline_14f2Folders.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural 
./run_pipeline_14f2Folders.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural  
./run_pipeline_14f2Folders.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc  
./run_pipeline_14f2Folders.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal  
./run_pipeline_14f2Folders.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural  
./run_pipeline_14f2Folders.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural  
./run_pipeline_14f2Folders.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc  

#>>> colon
./run_pipeline_14f2Folders.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss  # 

#>>> liver
./run_pipeline_14f2Folders.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc  #  
./run_pipeline_14f2Folders.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  #  
#./run_pipeline_14f2Folders.sh GSE58752_liver_40kb TCGAlihc_norm_lihc    # update name
#./run_pipeline_14f2Folders.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1   # update name  
./run_pipeline_14f2Folders.sh LI_40kb TCGAlihc_norm_lihc  # 
./run_pipeline_14f2Folders.sh LI_40kb TCGAlihc_wt_mutCTNNB1  # 


#>>> leukemia
./run_pipeline_14f2Folders.sh K562_40kb TCGAlaml_wt_mutFLT3  #  

#>>> pancreas
./run_pipeline_14f2Folders.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS # 



######
# additional datasets (already run 0-4)
######

#>>> lung
./run_pipeline_14f2Folders.sh LG1_40kb TCGAluad_mutKRAS_mutEGFR  # 
./run_pipeline_14f2Folders.sh LG1_40kb TCGAluad_nonsmoker_smoker  # 
./run_pipeline_14f2Folders.sh LG1_40kb TCGAluad_norm_luad    # 
./run_pipeline_14f2Folders.sh LG1_40kb TCGAluad_wt_mutKRAS    # 
./run_pipeline_14f2Folders.sh LG1_40kb TCGAlusc_norm_lusc    # 

./run_pipeline_14f2Folders.sh LG2_40kb TCGAluad_mutKRAS_mutEGFR  # 
./run_pipeline_14f2Folders.sh LG2_40kb TCGAluad_nonsmoker_smoker   # 
./run_pipeline_14f2Folders.sh LG2_40kb TCGAluad_norm_luad   # 
./run_pipeline_14f2Folders.sh LG2_40kb TCGAluad_wt_mutKRAS    # 
./run_pipeline_14f2Folders.sh LG2_40kb TCGAlusc_norm_lusc   # 
 
#>>> pancreas
./run_pipeline_14f2Folders.sh PA2_40kb TCGApaad_wt_mutKRAS # 

./run_pipeline_14f2Folders.sh PA3_40kb TCGApaad_wt_mutKRAS # 

./run_pipeline_14f2Folders.sh GSE118588_Panc_beta_40kb TCGApaad_wt_mutKRAS # 


#>>> breast
./run_pipeline_14f2Folders.sh Barutcu_MCF-10A_40kb TCGAbrca_lum_bas #  
./run_pipeline_14f2Folders.sh HMEC_40kb TCGAbrca_lum_bas  #  
./run_pipeline_14f2Folders.sh GSE109229_BT474_40kb TCGAbrca_lum_bas  #   
./run_pipeline_14f2Folders.sh GSE109229_SKBR3_40kb TCGAbrca_lum_bas  #  



#>>> colon
./run_pipeline_14f2Folders.sh Rao_HCT-116_2017_40kb TCGAcoad_msi_mss # 

./run_pipeline_14f2Folders.sh ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss # 


#>>> kidney
./run_pipeline_14f2Folders.sh GSE99051_786_O_40kb TCGAkich_norm_kich # 

#>>> prostate 
./run_pipeline_14f2Folders.sh GSE118514_22Rv1_40kb TCGAprad_norm_prad # 





