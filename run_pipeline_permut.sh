#!/usr/bin/bash

# ./run_pipeline_permut.sh

# run ele


## LIVER
 ./run_pipeline.sh GSE105381_HepG2_PERMUTG2T_40kb TCGAlihc_norm_lihc  
 ./run_pipeline.sh GSE105381_HepG2_PERMUTG2T_40kb TCGAlihc_wt_mutCTNNB1  
 ./run_pipeline.sh LI_PERMUTG2T_40kb TCGAlihc_norm_lihc 
 ./run_pipeline.sh LI_PERMUTG2T_40kb TCGAlihc_wt_mutCTNNB1
##BREAST
 ./run_pipeline.sh ENCSR549MGQ_T47D_PERMUTG2T_40kb TCGAbrca_lum_bas 
 ./run_pipeline.sh Barutcu_MCF-7_PERMUTG2T_40kb TCGAbrca_lum_bas
 ./run_pipeline.sh Barutcu_MCF-10A_PERMUTG2T_40kb TCGAbrca_lum_bas
 ./run_pipeline.sh GSE109229_BT474_PERMUTG2T_40kb TCGAbrca_lum_bas
 ./run_pipeline.sh GSE109229_SKBR3_PERMUTG2T_40kb TCGAbrca_lum_bas
 ./run_pipeline.sh HMEC_PERMUTG2T_40kb TCGAbrca_lum_bas
## KIDNEY
 ./run_pipeline.sh ENCSR079VIJ_G401_PERMUTG2T_40kb TCGAkich_norm_kich 
 ./run_pipeline.sh ENCSR401TBQ_Caki2_PERMUTG2T_40kb TCGAkich_norm_kich
 ./run_pipeline.sh GSE99051_786_O_PERMUTG2T_40kb TCGAkich_norm_kich
## SKIN
 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_lowInf_highInf
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_lowInf_highInf 
 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_wt_mutBRAF
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_wt_mutBRAF
 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_wt_mutCTNNB1	
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_wt_mutCTNNB1	
## LUNG
 ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_norm_luad  
 ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR 
 ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker 
 ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_wt_mutKRAS 
 ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAlusc_norm_lusc  
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker 
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_wt_mutKRAS 
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAlusc_norm_lusc  
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_norm_luad 
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR 
## PANCREAS
 ./run_pipeline.sh Panc1_rep12_PERMUTG2T_40kb TCGApaad_wt_mutKRAS 
 ./run_pipeline.sh GSE118588_Panc_beta_PERMUTG2T_40kb TCGApaad_wt_mutKRAS 


