#!/usr/bin/bash

######
# usual datasets => run on electron
######
#>>> kidney
./run_pipeline_8fcc.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich  # 
#>>> skin
./run_pipeline_8fcc.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf  #  
./run_pipeline_8fcc.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF  #  
./run_pipeline_8fcc.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1  #  

./run_pipeline_8fcc.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf  # 
./run_pipeline_8fcc.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF   #  
./run_pipeline_8fcc.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1   #  


# prostate
./run_pipeline_8fcc.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad   # 

./run_pipeline_8fcc.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad  # 

#./run_pipeline_8fcc.sh GSE73782_PC3_40kb TCGAprad_norm_prad  # 

#>>> kidney
./run_pipeline_8fcc.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich  # 

#>>> lung
./run_pipeline_8fcc.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR  # 
./run_pipeline_8fcc.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker  #   
./run_pipeline_8fcc.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad  # 
./run_pipeline_8fcc.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS  #  
./run_pipeline_8fcc.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc  #  #  
./run_pipeline_8fcc.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  #   
./run_pipeline_8fcc.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker  #  
./run_pipeline_8fcc.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad  #  
 # => from dataset below, added set.seed in TAD_DE_utils_fasterPermut 
./run_pipeline_8fcc.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS  #  
./run_pipeline_8fcc.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc  #  

#>>> breast
./run_pipeline_8fcc.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas  # 

#./run_pipeline_8fcc.sh GSE75070_MCF-7_shNS_40kb TCGAbrca_lum_bas  # 
./run_pipeline_8fcc.sh Barutcu_MCF-7_40kb TCGAbrca_lum_bas # 

#>>> cerebellum and spinal cord

./run_pipeline_8fcc.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal # 
./run_pipeline_8fcc.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural # 
./run_pipeline_8fcc.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural  # 
./run_pipeline_8fcc.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc  # 
./run_pipeline_8fcc.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal  # 
./run_pipeline_8fcc.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural  # 
./run_pipeline_8fcc.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural  # 
./run_pipeline_8fcc.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc  # 

#>>> colon
./run_pipeline_8fcc.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss  # neu

#>>> liver
./run_pipeline_8fcc.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc  #  neu
./run_pipeline_8fcc.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  #  neu
#./run_pipeline_8fcc.sh GSE58752_liver_40kb TCGAlihc_norm_lihc    # update name
#./run_pipeline_8fcc.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1   # update name  
./run_pipeline_8fcc.sh LI_40kb TCGAlihc_norm_lihc  # neu
./run_pipeline_8fcc.sh LI_40kb TCGAlihc_wt_mutCTNNB1  # neu


#>>> leukemia
./run_pipeline_8fcc.sh K562_40kb TCGAlaml_wt_mutFLT3  #  neu

#>>> pancreas
./run_pipeline_8fcc.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS # neu


