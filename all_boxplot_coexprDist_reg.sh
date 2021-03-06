#!/usr/bin/bash

# ./all_boxplot_coexprDist_reg.sh

start_time=$(date -R)   

scriptCoexpr="AUC_coexprDist_reg_boxplot.R"


# should have hicds + expr ds !!!

all_TAD_files_ds=(
#"ENCSR079VIJ_G401_40kb TCGAkich_norm_kich"  
#"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf" 
#"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF"  
#"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1"  
#"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf"
#"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF"
#"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1"  
#"ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad"  
#"GSE118514_RWPE1_40kb TCGAprad_norm_prad"  
#"ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich"  
#"ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR"  
#"ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker"  
#"ENCSR444WCZ_A549_40kb TCGAluad_norm_luad"  
#"ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS"  
#"ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc"  
#"ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR"  
#"ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker"  
#"ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad"  
#"ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS"  
#"ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc"
#"ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas"  
#"Barutcu_MCF-7_40kb TCGAbrca_lum_bas" 
#"GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal"  
#"GSE105194_cerebellum_40kb TCGAgbm_classical_neural" 
#"GSE105194_cerebellum_40kb TCGAgbm_classical_proneural"  
#"GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc"  
#"GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal"  
#"GSE105194_spinal_cord_40kb TCGAgbm_classical_neural"  
#"GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural"  
#"GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc"  
#"GSE105318_DLD1_40kb TCGAcoad_msi_mss"  
#"GSE105381_HepG2_40kb TCGAlihc_norm_lihc"  
#"GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1"  
#"LI_40kb TCGAlihc_norm_lihc"  
#"LI_40kb TCGAlihc_wt_mutCTNNB1"  
#"K562_40kb TCGAlaml_wt_mutFLT3"  
#"Panc1_rep12_40kb TCGApaad_wt_mutKRAS" 
#"LG1_40kb TCGAluad_mutKRAS_mutEGFR"
#"LG1_40kb TCGAluad_nonsmoker_smoker"  
#"LG1_40kb TCGAluad_norm_luad"  
#"LG1_40kb TCGAluad_wt_mutKRAS"  
#"LG1_40kb TCGAlusc_norm_lusc"
#"LG2_40kb TCGAluad_mutKRAS_mutEGFR"
#"LG2_40kb TCGAluad_nonsmoker_smoker"  
#"LG2_40kb TCGAluad_norm_luad"  
#"LG2_40kb TCGAluad_wt_mutKRAS"  
#"LG2_40kb TCGAlusc_norm_lusc"
#"PA2_40kb TCGApaad_wt_mutKRAS" 
#"PA3_40kb TCGApaad_wt_mutKRAS" 
#"GSE118588_Panc_beta_40kb TCGApaad_wt_mutKRAS" 
#"Barutcu_MCF-10A_40kb TCGAbrca_lum_bas" 
#"HMEC_40kb TCGAbrca_lum_bas" 
#"GSE109229_BT474_40kb TCGAbrca_lum_bas" 
#"GSE109229_SKBR3_40kb TCGAbrca_lum_bas" 
#"Rao_HCT-116_2017_40kb TCGAcoad_msi_mss" 
#"ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss"
#"GSE99051_786_O_40kb TCGAkich_norm_kich"
#"GSE118514_22Rv1_40kb TCGAprad_norm_prad"
#"ENCSR079VIJ_G401_RANDOMMIDPOS_40kb TCGAkich_norm_kich"  
#"ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb TCGAskcm_lowInf_highInf" 
#"ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb TCGAskcm_wt_mutBRAF"  
#"ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb TCGAskcm_wt_mutCTNNB1"  
#"ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb TCGAskcm_lowInf_highInf"
#"ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb TCGAskcm_wt_mutBRAF"
#"ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb TCGAskcm_wt_mutCTNNB1"  
#"ENCSR346DCU_LNCaP_RANDOMMIDPOS_40kb TCGAprad_norm_prad"  
#"GSE118514_RWPE1_RANDOMMIDPOS_40kb TCGAprad_norm_prad"  
#"ENCSR401TBQ_Caki2_RANDOMMIDPOS_40kb TCGAkich_norm_kich"  
#"ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR"  
#"ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker"  
#"ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_norm_luad"  
#"ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS"  
#"ENCSR444WCZ_A549_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc"  
#"ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR"  
#"ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker"  
#"ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_norm_luad"  
#"ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS"  
#"ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc"
#"ENCSR549MGQ_T47D_RANDOMMIDPOS_40kb TCGAbrca_lum_bas"  
#"Barutcu_MCF-7_RANDOMMIDPOS_40kb TCGAbrca_lum_bas" 
#"GSE105194_cerebellum_RANDOMMIDPOS_40kb TCGAgbm_classical_mesenchymal"  
#"GSE105194_cerebellum_RANDOMMIDPOS_40kb TCGAgbm_classical_neural" 
#"GSE105194_cerebellum_RANDOMMIDPOS_40kb TCGAgbm_classical_proneural"  
#"GSE105194_cerebellum_RANDOMMIDPOS_40kb TCGAlgg_IDHwt_IDHmutnc"  
#"GSE105194_spinal_cord_RANDOMMIDPOS_40kb TCGAgbm_classical_mesenchymal"  
#"GSE105194_spinal_cord_RANDOMMIDPOS_40kb TCGAgbm_classical_neural"  
#"GSE105194_spinal_cord_RANDOMMIDPOS_40kb TCGAgbm_classical_proneural"  
#"GSE105194_spinal_cord_RANDOMMIDPOS_40kb TCGAlgg_IDHwt_IDHmutnc"  
#"GSE105318_DLD1_RANDOMMIDPOS_40kb TCGAcoad_msi_mss"  
#"GSE105381_HepG2_RANDOMMIDPOS_40kb TCGAlihc_norm_lihc"  
#"GSE105381_HepG2_RANDOMMIDPOS_40kb TCGAlihc_wt_mutCTNNB1"  
#"LI_RANDOMMIDPOS_40kb TCGAlihc_norm_lihc"  
#"LI_RANDOMMIDPOS_40kb TCGAlihc_wt_mutCTNNB1"  
#"K562_RANDOMMIDPOS_40kb TCGAlaml_wt_mutFLT3"  
#"Panc1_rep12_RANDOMMIDPOS_40kb TCGApaad_wt_mutKRAS" 
#"LG1_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR"
#"LG1_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker"  
#"LG1_RANDOMMIDPOS_40kb TCGAluad_norm_luad"  
#"LG1_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS"  
#"LG1_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc"
#"LG2_RANDOMMIDPOS_40kb TCGAluad_mutKRAS_mutEGFR"
#"LG2_RANDOMMIDPOS_40kb TCGAluad_nonsmoker_smoker"  
#"LG2_RANDOMMIDPOS_40kb TCGAluad_norm_luad"  
#"LG2_RANDOMMIDPOS_40kb TCGAluad_wt_mutKRAS"  
#"LG2_RANDOMMIDPOS_40kb TCGAlusc_norm_lusc"
#"PA2_RANDOMMIDPOS_40kb TCGApaad_wt_mutKRAS" 
#"PA3_RANDOMMIDPOS_40kb TCGApaad_wt_mutKRAS" 
#"GSE118588_Panc_beta_RANDOMMIDPOS_40kb TCGApaad_wt_mutKRAS" 
#"Barutcu_MCF-10A_RANDOMMIDPOS_40kb TCGAbrca_lum_bas" 
#"HMEC_RANDOMMIDPOS_40kb TCGAbrca_lum_bas" 
#"GSE109229_BT474_RANDOMMIDPOS_40kb TCGAbrca_lum_bas" 
#"GSE109229_SKBR3_RANDOMMIDPOS_40kb TCGAbrca_lum_bas" 
#"Rao_HCT-116_2017_RANDOMMIDPOS_40kb TCGAcoad_msi_mss" 
#"ENCSR504OTV_transverse_colon_RANDOMMIDPOS_40kb TCGAcoad_msi_mss"
#"GSE99051_786_O_RANDOMMIDPOS_40kb TCGAkich_norm_kich"
#"GSE118514_22Rv1_RANDOMMIDPOS_40kb TCGAprad_norm_prad"
"ENCSR079VIJ_G401_RANDOMMIDPOSSTRICT_40kb TCGAkich_norm_kich"  
"ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSSTRICT_40kb TCGAskcm_lowInf_highInf" 
"ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSSTRICT_40kb TCGAskcm_wt_mutBRAF"  
"ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSSTRICT_40kb TCGAskcm_wt_mutCTNNB1"  
"ENCSR862OGI_RPMI-7951_RANDOMMIDPOSSTRICT_40kb TCGAskcm_lowInf_highInf"
"ENCSR862OGI_RPMI-7951_RANDOMMIDPOSSTRICT_40kb TCGAskcm_wt_mutBRAF"
"ENCSR862OGI_RPMI-7951_RANDOMMIDPOSSTRICT_40kb TCGAskcm_wt_mutCTNNB1"  
"ENCSR346DCU_LNCaP_RANDOMMIDPOSSTRICT_40kb TCGAprad_norm_prad"  
"GSE118514_RWPE1_RANDOMMIDPOSSTRICT_40kb TCGAprad_norm_prad"  
"ENCSR401TBQ_Caki2_RANDOMMIDPOSSTRICT_40kb TCGAkich_norm_kich"  
"ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAluad_mutKRAS_mutEGFR"  
"ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAluad_nonsmoker_smoker"  
"ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAluad_norm_luad"  
"ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAluad_wt_mutKRAS"  
"ENCSR444WCZ_A549_RANDOMMIDPOSSTRICT_40kb TCGAlusc_norm_lusc"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAluad_mutKRAS_mutEGFR"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAluad_nonsmoker_smoker"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAluad_norm_luad"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAluad_wt_mutKRAS"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSSTRICT_40kb TCGAlusc_norm_lusc"
"ENCSR549MGQ_T47D_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas"  
"Barutcu_MCF-7_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas" 
"GSE105194_cerebellum_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_mesenchymal"  
"GSE105194_cerebellum_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_neural" 
"GSE105194_cerebellum_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_proneural"  
"GSE105194_cerebellum_RANDOMMIDPOSSTRICT_40kb TCGAlgg_IDHwt_IDHmutnc"  
"GSE105194_spinal_cord_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_mesenchymal"  
"GSE105194_spinal_cord_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_neural"  
"GSE105194_spinal_cord_RANDOMMIDPOSSTRICT_40kb TCGAgbm_classical_proneural"  
"GSE105194_spinal_cord_RANDOMMIDPOSSTRICT_40kb TCGAlgg_IDHwt_IDHmutnc"  
"GSE105318_DLD1_RANDOMMIDPOSSTRICT_40kb TCGAcoad_msi_mss"  
"GSE105381_HepG2_RANDOMMIDPOSSTRICT_40kb TCGAlihc_norm_lihc"  
"GSE105381_HepG2_RANDOMMIDPOSSTRICT_40kb TCGAlihc_wt_mutCTNNB1"  
"LI_RANDOMMIDPOSSTRICT_40kb TCGAlihc_norm_lihc"  
"LI_RANDOMMIDPOSSTRICT_40kb TCGAlihc_wt_mutCTNNB1"  
"K562_RANDOMMIDPOSSTRICT_40kb TCGAlaml_wt_mutFLT3"  
"Panc1_rep12_RANDOMMIDPOSSTRICT_40kb TCGApaad_wt_mutKRAS" 
"LG1_RANDOMMIDPOSSTRICT_40kb TCGAluad_mutKRAS_mutEGFR"
"LG1_RANDOMMIDPOSSTRICT_40kb TCGAluad_nonsmoker_smoker"  
"LG1_RANDOMMIDPOSSTRICT_40kb TCGAluad_norm_luad"  
"LG1_RANDOMMIDPOSSTRICT_40kb TCGAluad_wt_mutKRAS"  
"LG1_RANDOMMIDPOSSTRICT_40kb TCGAlusc_norm_lusc"
"LG2_RANDOMMIDPOSSTRICT_40kb TCGAluad_mutKRAS_mutEGFR"
"LG2_RANDOMMIDPOSSTRICT_40kb TCGAluad_nonsmoker_smoker"  
"LG2_RANDOMMIDPOSSTRICT_40kb TCGAluad_norm_luad"  
"LG2_RANDOMMIDPOSSTRICT_40kb TCGAluad_wt_mutKRAS"  
"LG2_RANDOMMIDPOSSTRICT_40kb TCGAlusc_norm_lusc"
"PA2_RANDOMMIDPOSSTRICT_40kb TCGApaad_wt_mutKRAS" 
"PA3_RANDOMMIDPOSSTRICT_40kb TCGApaad_wt_mutKRAS" 
"GSE118588_Panc_beta_RANDOMMIDPOSSTRICT_40kb TCGApaad_wt_mutKRAS" 
"Barutcu_MCF-10A_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas" 
"HMEC_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas" 
"GSE109229_BT474_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas" 
"GSE109229_SKBR3_RANDOMMIDPOSSTRICT_40kb TCGAbrca_lum_bas" 
"Rao_HCT-116_2017_RANDOMMIDPOSSTRICT_40kb TCGAcoad_msi_mss" 
"ENCSR504OTV_transverse_colon_RANDOMMIDPOSSTRICT_40kb TCGAcoad_msi_mss"
"GSE99051_786_O_RANDOMMIDPOSSTRICT_40kb TCGAkich_norm_kich"
"GSE118514_22Rv1_RANDOMMIDPOSSTRICT_40kb TCGAprad_norm_prad"
)

# Rscript AUC_coexprDist_tf_boxplot.R c3.tft ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad
# Rscript AUC_coexprDist_tf_boxplot.R chea3_lung ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad
for ds in "${all_TAD_files_ds[@]}"; do
    echo Rscript $scriptCoexpr $ds
    Rscript $scriptCoexpr $ds
done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

