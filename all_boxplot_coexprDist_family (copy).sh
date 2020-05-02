#!/usr/bin/bash

# ./all_boxplot_coexprDist_family.sh

start_time=$(date -R)   

scriptCoexpr="AUC_coexprDist_family_boxplot.R"

# should have hicds + expr ds !!!

all_TAD_files_ds=(
"ENCSR079VIJ_G401_RANDOMMIDPOSDISC_40kb TCGAkich_norm_kich"  
"ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSDISC_40kb TCGAskcm_lowInf_highInf" 
"ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSDISC_40kb TCGAskcm_wt_mutBRAF"  
"ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSDISC_40kb TCGAskcm_wt_mutCTNNB1"  
"ENCSR862OGI_RPMI-7951_RANDOMMIDPOSDISC_40kb TCGAskcm_lowInf_highInf"
"ENCSR862OGI_RPMI-7951_RANDOMMIDPOSDISC_40kb TCGAskcm_wt_mutBRAF"
"ENCSR862OGI_RPMI-7951_RANDOMMIDPOSDISC_40kb TCGAskcm_wt_mutCTNNB1"  
"ENCSR346DCU_LNCaP_RANDOMMIDPOSDISC_40kb TCGAprad_norm_prad"  
"GSE118514_RWPE1_RANDOMMIDPOSDISC_40kb TCGAprad_norm_prad"  
"ENCSR401TBQ_Caki2_RANDOMMIDPOSDISC_40kb TCGAkich_norm_kich"  
"ENCSR444WCZ_A549_RANDOMMIDPOSDISC_40kb TCGAluad_mutKRAS_mutEGFR"  
"ENCSR444WCZ_A549_RANDOMMIDPOSDISC_40kb TCGAluad_nonsmoker_smoker"  
"ENCSR444WCZ_A549_RANDOMMIDPOSDISC_40kb TCGAluad_norm_luad"  
"ENCSR444WCZ_A549_RANDOMMIDPOSDISC_40kb TCGAluad_wt_mutKRAS"  
"ENCSR444WCZ_A549_RANDOMMIDPOSDISC_40kb TCGAlusc_norm_lusc"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSDISC_40kb TCGAluad_mutKRAS_mutEGFR"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSDISC_40kb TCGAluad_nonsmoker_smoker"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSDISC_40kb TCGAluad_norm_luad"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSDISC_40kb TCGAluad_wt_mutKRAS"  
"ENCSR489OCU_NCI-H460_RANDOMMIDPOSDISC_40kb TCGAlusc_norm_lusc"
"ENCSR549MGQ_T47D_RANDOMMIDPOSDISC_40kb TCGAbrca_lum_bas"  
"Barutcu_MCF-7_RANDOMMIDPOSDISC_40kb TCGAbrca_lum_bas" 
"GSE105194_cerebellum_RANDOMMIDPOSDISC_40kb TCGAgbm_classical_mesenchymal"  
"GSE105194_cerebellum_RANDOMMIDPOSDISC_40kb TCGAgbm_classical_neural" 
"GSE105194_cerebellum_RANDOMMIDPOSDISC_40kb TCGAgbm_classical_proneural"  
"GSE105194_cerebellum_RANDOMMIDPOSDISC_40kb TCGAlgg_IDHwt_IDHmutnc"  
"GSE105194_spinal_cord_RANDOMMIDPOSDISC_40kb TCGAgbm_classical_mesenchymal"  
"GSE105194_spinal_cord_RANDOMMIDPOSDISC_40kb TCGAgbm_classical_neural"  
"GSE105194_spinal_cord_RANDOMMIDPOSDISC_40kb TCGAgbm_classical_proneural"  
"GSE105194_spinal_cord_RANDOMMIDPOSDISC_40kb TCGAlgg_IDHwt_IDHmutnc"  
"GSE105318_DLD1_RANDOMMIDPOSDISC_40kb TCGAcoad_msi_mss"  
"GSE105381_HepG2_RANDOMMIDPOSDISC_40kb TCGAlihc_norm_lihc"  
"GSE105381_HepG2_RANDOMMIDPOSDISC_40kb TCGAlihc_wt_mutCTNNB1"  
"LI_RANDOMMIDPOSDISC_40kb TCGAlihc_norm_lihc"  
"LI_RANDOMMIDPOSDISC_40kb TCGAlihc_wt_mutCTNNB1"  
"K562_RANDOMMIDPOSDISC_40kb TCGAlaml_wt_mutFLT3"  
"Panc1_rep12_RANDOMMIDPOSDISC_40kb TCGApaad_wt_mutKRAS" 
"LG1_RANDOMMIDPOSDISC_40kb TCGAluad_mutKRAS_mutEGFR"
"LG1_RANDOMMIDPOSDISC_40kb TCGAluad_nonsmoker_smoker"  
"LG1_RANDOMMIDPOSDISC_40kb TCGAluad_norm_luad"  
"LG1_RANDOMMIDPOSDISC_40kb TCGAluad_wt_mutKRAS"  
"LG1_RANDOMMIDPOSDISC_40kb TCGAlusc_norm_lusc"
"LG2_RANDOMMIDPOSDISC_40kb TCGAluad_mutKRAS_mutEGFR"
"LG2_RANDOMMIDPOSDISC_40kb TCGAluad_nonsmoker_smoker"  
"LG2_RANDOMMIDPOSDISC_40kb TCGAluad_norm_luad"  
"LG2_RANDOMMIDPOSDISC_40kb TCGAluad_wt_mutKRAS"  
"LG2_RANDOMMIDPOSDISC_40kb TCGAlusc_norm_lusc"
"PA2_RANDOMMIDPOSDISC_40kb TCGApaad_wt_mutKRAS" 
"PA3_RANDOMMIDPOSDISC_40kb TCGApaad_wt_mutKRAS" 
"GSE118588_Panc_beta_RANDOMMIDPOSDISC_40kb TCGApaad_wt_mutKRAS" 
"Barutcu_MCF-10A_RANDOMMIDPOSDISC_40kb TCGAbrca_lum_bas" 
"HMEC_RANDOMMIDPOSDISC_40kb TCGAbrca_lum_bas" 
"GSE109229_BT474_RANDOMMIDPOSDISC_40kb TCGAbrca_lum_bas" 
"GSE109229_SKBR3_RANDOMMIDPOSDISC_40kb TCGAbrca_lum_bas" 
"Rao_HCT-116_2017_RANDOMMIDPOSDISC_40kb TCGAcoad_msi_mss" 
"ENCSR504OTV_transverse_colon_RANDOMMIDPOSDISC_40kb TCGAcoad_msi_mss"
"GSE99051_786_O_RANDOMMIDPOSDISC_40kb TCGAkich_norm_kich"
"GSE118514_22Rv1_RANDOMMIDPOSDISC_40kb TCGAprad_norm_prad"
)
#all_TAD_files_ds=(
#"ENCSR489OCU_NCI-H460_RANDOMMIDPOSDISC_40kb TCGAluad_norm_luad"
#"ENCSR489OCU_NCI-H460_RANDOMNBRGENES_40kb TCGAluad_norm_luad"
#"ENCSR489OCU_NCI-H460_RANDOMSHIFT_40kb TCGAluad_norm_luad"
#"ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_norm_luad"
#)

# Rscript AUC_coexprDist_family_boxplot.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad

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
