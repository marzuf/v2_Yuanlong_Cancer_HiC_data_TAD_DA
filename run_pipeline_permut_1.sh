#!/usr/bin/bash

# ./run_pipeline_permut_1.sh


### LIVER
 ./run_pipeline.sh GSE105381_HepG2_RANDOMMIDPOS_40kb TCGAlihc_norm_lihc  
 ./run_pipeline.sh GSE105381_HepG2_RANDOMMIDPOS_40kb TCGAlihc_wt_mutCTNNB1  
 ./run_pipeline.sh LI_RANDOMMIDPOS_40kb TCGAlihc_norm_lihc 
 ./run_pipeline.sh LI_RANDOMMIDPOS_40kb TCGAlihc_wt_mutCTNNB1
#BREAST
 ./run_pipeline.sh ENCSR549MGQ_T47D_RANDOMMIDPOS_40kb TCGAbrca_lum_bas 
 ./run_pipeline.sh Barutcu_MCF-7_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
 ./run_pipeline.sh Barutcu_MCF-10A_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
 ./run_pipeline.sh GSE109229_BT474_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
 ./run_pipeline.sh GSE109229_SKBR3_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
 ./run_pipeline.sh HMEC_RANDOMMIDPOS_40kb TCGAbrca_lum_bas
# KIDNEY
 ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMMIDPOS_40kb TCGAkich_norm_kich 
 ./run_pipeline.sh ENCSR401TBQ_Caki2_RANDOMMIDPOS_40kb TCGAkich_norm_kich
 ./run_pipeline.sh GSE99051_786_O_RANDOMMIDPOS_40kb TCGAkich_norm_kich
# SKIN
 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb TCGAskcm_lowInf_highInf
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb TCGAskcm_lowInf_highInf 
 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb TCGAskcm_wt_mutBRAF
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb TCGAskcm_wt_mutBRAF
 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb TCGAskcm_wt_mutCTNNB1	
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb TCGAskcm_wt_mutCTNNB1	
# LUNG
 ./run_pipeline.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_norm_luad  
 ./run_pipeline.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR 
 ./run_pipeline.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker 
 ./run_pipeline.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS 
 ./run_pipeline.sh ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc  
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker 
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS 
 ./run_pipeline.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc  
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_norm_luad 
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR 







