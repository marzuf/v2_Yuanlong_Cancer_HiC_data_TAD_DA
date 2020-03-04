#!/usr/bin/bash

# ./run_pipeline_permut.sh


## LIVER
 ./run_pipeline6.sh GSE105381_HepG2_RANDOMMIDPOS_40kb TCGAlihc_norm_lihc  
 ./run_pipeline6.sh GSE105381_HepG2_RANDOMMIDPOS_40kb TCGAlihc_wt_mutCTNNB1  
 ./run_pipeline6.sh LI_RANDOMMIDPOS_40kb TCGAlihc_norm_lihc 
 ./run_pipeline6.sh LI_RANDOMMIDPOS_40kb TCGAlihc_wt_mutCTNNB1
#BREAST
 ./run_pipeline6.sh ENCSR549MGQ_T47D_RANDOMMIDPOS_40kb TCGAbrca_lum_bas 
 ./run_pipeline6.sh Barutcu_MCF-7_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
 ./run_pipeline6.sh Barutcu_MCF-10A_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
 ./run_pipeline6.sh GSE109229_BT474_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
 ./run_pipeline6.sh GSE109229_SKBR3_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
 ./run_pipeline6.sh HMEC_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
# KIDNEY
 ./run_pipeline6.sh ENCSR079VIJ_G401_RANDOMMIDPOS_40kb TCGAkich_norm_kich 
 ./run_pipeline6.sh ENCSR401TBQ_Caki2_RANDOMMIDPOS_40kb TCGAkich_norm_kich
 ./run_pipeline6.sh GSE99051_786_O_RANDOMMIDPOS_40kb TCGAkich_norm_kich
# SKIN
 ./run_pipeline6.sh ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb TCGAskcm_lowInf_highInf
 ./run_pipeline6.sh ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb TCGAskcm_lowInf_highInf 
 ./run_pipeline6.sh ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb TCGAskcm_wt_mutBRAF
 ./run_pipeline6.sh ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb TCGAskcm_wt_mutBRAF
 ./run_pipeline6.sh ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb TCGAskcm_wt_mutCTNNB1	
 ./run_pipeline6.sh ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb TCGAskcm_wt_mutCTNNB1	
# LUNG
 ./run_pipeline6.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_norm_luad  
 ./run_pipeline6.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_norm_luad    
 ./run_pipeline6.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR 
 ./run_pipeline6.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR 
 ./run_pipeline6.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker 
 ./run_pipeline6.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker 
 ./run_pipeline6.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS 
 ./run_pipeline6.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS 
 ./run_pipeline6.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc  
 ./run_pipeline6.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc  
# PANCREAS
 ./run_pipeline6.sh Panc1_rep12_RANDOMMIDPOS_40kb TCGApaad_wt_mutKRAS 
 ./run_pipeline6.sh GSE118588_Panc_beta_RANDOMMIDPOS_40kb TCGApaad_wt_mutKRAS 
# PROSTATE	
 ./run_pipeline6.sh ENCSR346DCU_LNCaP_RANDOMMIDPOS_40kb TCGAprad_norm_prad    
 ./run_pipeline6.sh GSE118514_22Rv1_RANDOMMIDPOS_40kb TCGAprad_norm_prad       
 ./run_pipeline6.sh GSE118514_RWPE1_RANDOMMIDPOS_40kb TCGAprad_norm_prad       
# GBM
 ./run_pipeline6.sh GSE105194_cerebellum_RANDOMMIDPOS_40kb TCGAgbm_classical_mesenchymal 
 ./run_pipeline6.sh GSE105194_cerebellum_RANDOMMIDPOS_40kb TCGAgbm_classical_neural 
 ./run_pipeline6.sh GSE105194_cerebellum_RANDOMMIDPOS_40kb TCGAgbm_classical_proneural 
 ./run_pipeline6.sh GSE105194_cerebellum_RANDOMMIDPOS_40kb TCGAlgg_IDHwt_IDHmutnc 
 ./run_pipeline6.sh GSE105194_spinal_cord_RANDOMMIDPOS_40kb TCGAgbm_classical_mesenchymal 
 ./run_pipeline6.sh GSE105194_spinal_cord_RANDOMMIDPOS_40kb TCGAgbm_classical_neural  
 ./run_pipeline6.sh GSE105194_spinal_cord_RANDOMMIDPOS_40kb TCGAgbm_classical_proneural 
 ./run_pipeline6.sh GSE105194_spinal_cord_RANDOMMIDPOS_40kb TCGAlgg_IDHwt_IDHmutnc 
# COLORECTAL
 ./run_pipeline6.sh GSE105318_DLD1_RANDOMMIDPOS_40kb TCGAcoad_msi_mss 
 ./run_pipeline6.sh ENCSR504OTV_transverse_colon_RANDOMMIDPOS_40kb TCGAcoad_msi_mss 
 ./run_pipeline6.sh Rao_HCT-116_2017_RANDOMMIDPOS_40kb TCGAcoad_msi_mss 
# LYMPHOBLAST
 ./run_pipeline6.sh K562_RANDOMMIDPOS_40kb TCGAlaml_wt_mutFLT3 
 ./run_pipeline6.sh LG1_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR
 ./run_pipeline6.sh LG1_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker
 ./run_pipeline6.sh LG1_RANDOMMIDPOS_40kb TCGAluad_norm_luad
 ./run_pipeline6.sh LG1_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS
 ./run_pipeline6.sh LG1_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc
 ./run_pipeline6.sh LG2_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR
 ./run_pipeline6.sh LG2_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker
 ./run_pipeline6.sh LG2_RANDOMMIDPOS_40kb TCGAluad_norm_luad
 ./run_pipeline6.sh LG2_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS
 ./run_pipeline6.sh LG2_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc
 ./run_pipeline6.sh PA2_RANDOMMIDPOS_40kb TCGApaad_wt_mutKRAS
 ./run_pipeline6.sh PA3_RANDOMMIDPOS_40kb TCGApaad_wt_mutKRAS








