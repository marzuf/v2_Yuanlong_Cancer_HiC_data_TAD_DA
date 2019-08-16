#!/usr/bin/bash

######
# usual datasets
######
#>>> kidney
./run_pipeline.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich  # 5,6,7 - el ok

#>>> skin
./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf  # 5,6,7 - el ok
./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF  # 5,6,7 - pos ok
./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1  # 5,6,7 - pos ok

./run_pipeline_56.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf  # neutrino - pip1 5,6 ok
./run_pipeline_56.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF   # neutrino - pip2 5,6 ok
./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1   # neutrino - pip1 


# prostate
./run_pipeline.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad   # neutrino - pip1 

./run_pipeline.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad  # el 

#./run_pipeline.sh GSE73782_PC3_40kb TCGAprad_norm_prad  # discard low quality

#>>> kidney
./run_pipeline.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich  # pos - ok

#>>> lung
./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR  # pos - ok
./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker  # el
./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad  # pos                   # => from datasets below, run with 5fastSavePermut and 6fastSave
./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS  
./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc  
./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  
./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker  
./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad  
./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS  
./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc  

#>>> breast
./run_pipeline.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas  

#./run_pipeline.sh GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas  # update name
./run_pipeline.sh Barutcu_MCF-7_40kb TCGAbrca_lum_bas 

#>>> cerebellum and spinal cord

./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal  
./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural 
./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural  
./run_pipeline.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc  
./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal  
./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural  
./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural  
./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc  

#>>> colon
./run_pipeline.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss  

#>>> liver
./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc  
./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  
#./run_pipeline.sh GSE58752_liver_40kb TCGAlihc_norm_lihc    # update name
#./run_pipeline.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1   # update name  
./run_pipeline.sh LI_40kb TCGAlihc_norm_lihc  
./run_pipeline.sh LI_40kb TCGAlihc_wt_mutCTNNB1  


#>>> leukemia
./run_pipeline.sh K562_40kb TCGAlaml_wt_mutFLT3  

#>>> pancreas
./run_pipeline.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS 


######
# additional datasets
######

#>>> lung
./run_pipeline.sh LG1_40kb TCGAluad_mutKRAS_mutEGFR  
./run_pipeline.sh LG1_40kb TCGAluad_nonsmoker_smoker  
./run_pipeline.sh LG1_40kb TCGAluad_norm_luad  
./run_pipeline.sh LG1_40kb TCGAluad_wt_mutKRAS  
./run_pipeline.sh LG1_40kb TCGAlusc_norm_lusc  

./run_pipeline.sh LG2_40kb TCGAluad_mutKRAS_mutEGFR  
./run_pipeline.sh LG2_40kb TCGAluad_nonsmoker_smoker  
./run_pipeline.sh LG2_40kb TCGAluad_norm_luad  
./run_pipeline.sh LG2_40kb TCGAluad_wt_mutKRAS  
./run_pipeline.sh LG2_40kb TCGAlusc_norm_lusc  

#>>> pancreas
./run_pipeline.sh PA2_40kb TCGApaad_wt_mutKRAS 

./run_pipeline.sh PA3_40kb TCGApaad_wt_mutKRAS 

./run_pipeline.sh GSE118588_Panc_beta_40kb TCGApaad_wt_mutKRAS 


#>>> breast
./run_pipeline.sh Barutcu_MCF-10A_40kb TCGAbrca_lum_bas 
./run_pipeline.sh HMEC_40kb TCGAbrca_lum_bas 
./run_pipeline.sh GSE109229_BT474_40kb TCGAbrca_lum_bas 
./run_pipeline.sh GSE109229_SKBR3_40kb TCGAbrca_lum_bas 



#>>> colon
./run_pipeline.sh Rao_HCT-116_2017_40kb TCGAcoad_msi_mss 

./run_pipeline.sh ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss


#>>> kidney
./run_pipeline.sh GSE99051_786_O_40kb TCGAkich_norm_kich

#>>> prostate
./run_pipeline.sh GSE118514_22Rv1_40kb TCGAprad_norm_prad




