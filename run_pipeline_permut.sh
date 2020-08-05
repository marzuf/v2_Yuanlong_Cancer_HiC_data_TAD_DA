#!/usr/bin/bash

# ./run_pipeline_permut.sh

# for graviton: 
#conda create -n Renv r-essentials r-base
# conda activate Renv
# if packages needed -> start R after "conda activate Renv" and "install.packages()" as usual

## LIVER
# ./run_pipeline.sh GSE105381_HepG2_PERMUTG2T_40kb TCGAlihc_norm_lihc   # run single el
# ./run_pipeline.sh GSE105381_HepG2_PERMUTG2T_40kb TCGAlihc_wt_mutCTNNB1  # run pos
# ./run_pipeline.sh LI_PERMUTG2T_40kb TCGAlihc_norm_lihc  # run pos
# ./run_pipeline.sh LI_PERMUTG2T_40kb TCGAlihc_wt_mutCTNNB1 # run pos
##BREAST
# ./run_pipeline.sh ENCSR549MGQ_T47D_PERMUTG2T_40kb TCGAbrca_lum_bas  # run pos
# ./run_pipeline.sh Barutcu_MCF-7_PERMUTG2T_40kb TCGAbrca_lum_bas # run pos
# ./run_pipeline.sh Barutcu_MCF-10A_PERMUTG2T_40kb TCGAbrca_lum_bas # run pos
# ./run_pipeline.sh GSE109229_BT474_PERMUTG2T_40kb TCGAbrca_lum_bas # run pos
# ./run_pipeline.sh GSE109229_SKBR3_PERMUTG2T_40kb TCGAbrca_lum_bas # run pos
# ./run_pipeline.sh HMEC_PERMUTG2T_40kb TCGAbrca_lum_bas # run pos
## KIDNEY
# ./run_pipeline.sh ENCSR079VIJ_G401_PERMUTG2T_40kb TCGAkich_norm_kich  # run pos
# ./run_pipeline.sh ENCSR401TBQ_Caki2_PERMUTG2T_40kb TCGAkich_norm_kich # run pos
# ./run_pipeline.sh GSE99051_786_O_PERMUTG2T_40kb TCGAkich_norm_kich # run pos
## SKIN
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_lowInf_highInf # run pos
# ./run_pipeline.sh ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_lowInf_highInf # run pos 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_wt_mutBRAF # run pos
# ./run_pipeline.sh ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_wt_mutBRAF # run pos
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_wt_mutCTNNB1 # run pos	
# ./run_pipeline.sh ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_wt_mutCTNNB1 # run pos	
## LUNG
# ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_norm_luad   # run pos
## ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_norm_luad     # already run
# ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR # run pos
## ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR # already run
# ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker  # run pos
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker # run pos


# ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_wt_mutKRAS # run grav 
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_wt_mutKRAS # run grav 
# ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAlusc_norm_lusc  # run grav 
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAlusc_norm_lusc  # run grav 
## PANCREAS
# ./run_pipeline.sh Panc1_rep12_PERMUTG2T_40kb TCGApaad_wt_mutKRAS # run grav 
# ./run_pipeline.sh GSE118588_Panc_beta_PERMUTG2T_40kb TCGApaad_wt_mutKRAS # run grav 
## PROSTATE	
# ./run_pipeline.sh ENCSR346DCU_LNCaP_PERMUTG2T_40kb TCGAprad_norm_prad    # run grav 
# ./run_pipeline.sh GSE118514_22Rv1_PERMUTG2T_40kb TCGAprad_norm_prad       # run grav 
# ./run_pipeline.sh GSE118514_RWPE1_PERMUTG2T_40kb TCGAprad_norm_prad       # run grav 


## GBM
 ./run_pipeline.sh GSE105194_cerebellum_PERMUTG2T_40kb TCGAgbm_classical_mesenchymal # run el
 ./run_pipeline.sh GSE105194_cerebellum_PERMUTG2T_40kb TCGAgbm_classical_neural # run el
 ./run_pipeline.sh GSE105194_cerebellum_PERMUTG2T_40kb TCGAgbm_classical_proneural # run el
 ./run_pipeline.sh GSE105194_cerebellum_PERMUTG2T_40kb TCGAlgg_IDHwt_IDHmutnc # run el
 ./run_pipeline.sh GSE105194_spinal_cord_PERMUTG2T_40kb TCGAgbm_classical_mesenchymal # run el
 ./run_pipeline.sh GSE105194_spinal_cord_PERMUTG2T_40kb TCGAgbm_classical_neural  # run el
 ./run_pipeline.sh GSE105194_spinal_cord_PERMUTG2T_40kb TCGAgbm_classical_proneural # run el
 ./run_pipeline.sh GSE105194_spinal_cord_PERMUTG2T_40kb TCGAlgg_IDHwt_IDHmutnc # run el
# COLORECTAL
 ./run_pipeline.sh GSE105318_DLD1_PERMUTG2T_40kb TCGAcoad_msi_mss # run el
 ./run_pipeline.sh ENCSR504OTV_transverse_colon_PERMUTG2T_40kb TCGAcoad_msi_mss # run el
 ./run_pipeline.sh Rao_HCT-116_2017_PERMUTG2T_40kb TCGAcoad_msi_mss # run el
# LYMPHOBLAST
 ./run_pipeline.sh K562_PERMUTG2T_40kb TCGAlaml_wt_mutFLT3 # run el
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR # run el
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker # run el
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_norm_luad # run el
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_wt_mutKRAS # run el
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAlusc_norm_lusc # run el
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR # run el
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker # run el
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_norm_luad # run el
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_wt_mutKRAS # run el
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAlusc_norm_lusc # run el
 ./run_pipeline.sh PA2_PERMUTG2T_40kb TCGApaad_wt_mutKRAS # run el
 ./run_pipeline.sh PA3_PERMUTG2T_40kb TCGApaad_wt_mutKRAS # run el








