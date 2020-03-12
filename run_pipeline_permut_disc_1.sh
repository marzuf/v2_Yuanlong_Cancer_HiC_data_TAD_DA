#!/usr/bin/bash

# ./run_pipeline_permut_disc_1.sh



# LIVER
 ./run_pipelineStrict.sh GSE105381_HepG2_RANDOMMIDPOSSTRICT_40kb TCGAlihc_norm_lihc  
 ./run_pipelineStrict.sh GSE105381_HepG2_RANDOMMIDPOSSTRICT_40kb TCGAlihc_wt_mutCTNNB1  
 ./run_pipelineStrict.sh LI_RANDOMMIDPOSSTRICT_40kb TCGAlihc_norm_lihc 
 ./run_pipelineStrict.sh LI_RANDOMMIDPOSSTRICT_40kb TCGAlihc_wt_mutCTNNB1
#BREAST
 ./run_pipelineStrict.sh ENCSR549MGQ_T47D_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas 
 ./run_pipelineStrict.sh Barutcu_MCF-7_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas
 ./run_pipelineStrict.sh Barutcu_MCF-10A_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas
 ./run_pipelineStrict.sh GSE109229_BT474_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas
 ./run_pipelineStrict.sh GSE109229_SKBR3_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas
 ./run_pipelineStrict.sh HMEC_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas
# KIDNEY
 ./run_pipelineStrict.sh ENCSR079VIJ_G401_RANDOMMIDPOSSTRICT_40kb TCGAkich_norm_kich 
 ./run_pipelineStrict.sh ENCSR401TBQ_Caki2_RANDOMMIDPOSSTRICT_40kb TCGAkich_norm_kich
 ./run_pipelineStrict.sh GSE99051_786_O_RANDOMMIDPOSSTRICT_40kb TCGAkich_norm_kich
# SKIN
 ./run_pipelineStrict.sh ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSSTRICT_40kb TCGAskcm_lowInf_highInf
 ./run_pipelineStrict.sh ENCSR862OGI_RPMI-7951_RANDOMMIDPOSSTRICT_40kb TCGAskcm_lowInf_highInf 
 ./run_pipelineStrict.sh ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSSTRICT_40kb TCGAskcm_wt_mutBRAF
 ./run_pipelineStrict.sh ENCSR862OGI_RPMI-7951_RANDOMMIDPOSSTRICT_40kb TCGAskcm_wt_mutBRAF
 ./run_pipelineStrict.sh ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSSTRICT_40kb TCGAskcm_wt_mutCTNNB1	
 ./run_pipelineStrict.sh ENCSR862OGI_RPMI-7951_RANDOMMIDPOSSTRICT_40kb TCGAskcm_wt_mutCTNNB1	
# LUNG
 ./run_pipelineStrict.sh ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAluad_norm_luad  
 ./run_pipelineStrict.sh ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAluad_mutKRAS_mutEGFR 
 ./run_pipelineStrict.sh ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAluad_nonsmoker_smoker 
 ./run_pipelineStrict.sh ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAluad_wt_mutKRAS 
 ./run_pipelineStrict.sh ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAlusc_norm_lusc  
 ./run_pipelineStrict.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAluad_nonsmoker_smoker 
 ./run_pipelineStrict.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAluad_wt_mutKRAS 
 ./run_pipelineStrict.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAlusc_norm_lusc  
 ./run_pipelineStrict.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAluad_norm_luad 
 ./run_pipelineStrict.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAluad_mutKRAS_mutEGFR 
# PANCREAS
 ./run_pipelineStrict.sh Panc1_rep12_RANDOMMIDPOSSTRICT_40kb TCGApaad_wt_mutKRAS 
 ./run_pipelineStrict.sh GSE118588_Panc_beta_RANDOMMIDPOSSTRICT_40kb TCGApaad_wt_mutKRAS 
## PROSTATE	
# ./run_pipelineStrict.sh ENCSR346DCU_LNCaP_RANDOMMIDPOSSTRICT_40kb TCGAprad_norm_prad    
# ./run_pipelineStrict.sh GSE118514_22Rv1_RANDOMMIDPOSSTRICT_40kb TCGAprad_norm_prad       
# ./run_pipelineStrict.sh GSE118514_RWPE1_RANDOMMIDPOSSTRICT_40kb TCGAprad_norm_prad       
## GBM
# ./run_pipelineStrict.sh GSE105194_cerebellum_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_mesenchymal 
# ./run_pipelineStrict.sh GSE105194_cerebellum_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_neural 
# ./run_pipelineStrict.sh GSE105194_cerebellum_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_proneural 
# ./run_pipelineStrict.sh GSE105194_cerebellum_RANDOMMIDPOSSTRICT_40kb TCGAlgg_IDHwt_IDHmutnc 
# ./run_pipelineStrict.sh GSE105194_spinal_cord_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_mesenchymal 
# ./run_pipelineStrict.sh GSE105194_spinal_cord_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_neural  
# ./run_pipelineStrict.sh GSE105194_spinal_cord_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_proneural 
# ./run_pipelineStrict.sh GSE105194_spinal_cord_RANDOMMIDPOSSTRICT_40kb TCGAlgg_IDHwt_IDHmutnc 
## COLORECTAL
# ./run_pipelineStrict.sh GSE105318_DLD1_RANDOMMIDPOSSTRICT_40kb TCGAcoad_msi_mss 
# ./run_pipelineStrict.sh ENCSR504OTV_transverse_colon_RANDOMMIDPOSSTRICT_40kb TCGAcoad_msi_mss 
# ./run_pipelineStrict.sh Rao_HCT-116_2017_RANDOMMIDPOSSTRICT_40kb TCGAcoad_msi_mss 
## LYMPHOBLAST
# ./run_pipelineStrict.sh K562_RANDOMMIDPOSSTRICT_40kb TCGAlaml_wt_mutFLT3 
# ./run_pipelineStrict.sh LG1_RANDOMMIDPOSSTRICT_40kb TCGAluad_mutKRAS_mutEGFR
# ./run_pipelineStrict.sh LG1_RANDOMMIDPOSSTRICT_40kb TCGAluad_nonsmoker_smoker
# ./run_pipelineStrict.sh LG1_RANDOMMIDPOSSTRICT_40kb TCGAluad_norm_luad
# ./run_pipelineStrict.sh LG1_RANDOMMIDPOSSTRICT_40kb TCGAluad_wt_mutKRAS
# ./run_pipelineStrict.sh LG1_RANDOMMIDPOSSTRICT_40kb TCGAlusc_norm_lusc
# ./run_pipelineStrict.sh LG2_RANDOMMIDPOSSTRICT_40kb TCGAluad_mutKRAS_mutEGFR
# ./run_pipelineStrict.sh LG2_RANDOMMIDPOSSTRICT_40kb TCGAluad_nonsmoker_smoker
# ./run_pipelineStrict.sh LG2_RANDOMMIDPOSSTRICT_40kb TCGAluad_norm_luad
# ./run_pipelineStrict.sh LG2_RANDOMMIDPOSSTRICT_40kb TCGAluad_wt_mutKRAS
# ./run_pipelineStrict.sh LG2_RANDOMMIDPOSSTRICT_40kb TCGAlusc_norm_lusc
# ./run_pipelineStrict.sh PA2_RANDOMMIDPOSSTRICT_40kb TCGApaad_wt_mutKRAS
# ./run_pipelineStrict.sh PA3_RANDOMMIDPOSSTRICT_40kb TCGApaad_wt_mutKRAS






