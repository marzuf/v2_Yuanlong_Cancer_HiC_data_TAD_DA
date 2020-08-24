#!/usr/bin/bash

# ./all_AUC_coexprDist_sortNoDup_otherTADfile_otherFamFile_1.sh # pos

start_time=$(date -R)   

script_name="../Yuanlong_Cancer_HiC_data_TAD_DA/AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R"

# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R ENCSR079VIJ_G401_PERMUTG2T_40kb TCGAkich_norm_kich hgnc
# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 
 
#############################################################
all_data=(
"ENCSR079VIJ_G401_PERMUTG2T_40kb TCGAkich_norm_kich"  
"ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_lowInf_highInf" 
"ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_wt_mutBRAF"  
"ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_wt_mutCTNNB1"  
"ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_lowInf_highInf"
"ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_wt_mutBRAF"
"ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_wt_mutCTNNB1"  
"ENCSR346DCU_LNCaP_PERMUTG2T_40kb TCGAprad_norm_prad"  
"GSE118514_RWPE1_PERMUTG2T_40kb TCGAprad_norm_prad"  
"ENCSR401TBQ_Caki2_PERMUTG2T_40kb TCGAkich_norm_kich"  
"ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR"  
"ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker"  
"ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_norm_luad"  
"ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_wt_mutKRAS"  
"ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAlusc_norm_lusc"  
"ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR"  
"ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker"  
"ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_norm_luad"  
"ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_wt_mutKRAS"  
"ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAlusc_norm_lusc"
"ENCSR549MGQ_T47D_PERMUTG2T_40kb TCGAbrca_lum_bas"  
"Barutcu_MCF-7_PERMUTG2T_40kb TCGAbrca_lum_bas" 
"GSE105194_cerebellum_PERMUTG2T_40kb TCGAgbm_classical_mesenchymal"  
"GSE105194_cerebellum_PERMUTG2T_40kb TCGAgbm_classical_neural" 
"GSE105194_cerebellum_PERMUTG2T_40kb TCGAgbm_classical_proneural"  
"GSE105194_cerebellum_PERMUTG2T_40kb TCGAlgg_IDHwt_IDHmutnc"  
"GSE105194_spinal_cord_PERMUTG2T_40kb TCGAgbm_classical_mesenchymal"  
"GSE105194_spinal_cord_PERMUTG2T_40kb TCGAgbm_classical_neural"  
"GSE105194_spinal_cord_PERMUTG2T_40kb TCGAgbm_classical_proneural"  
"GSE105194_spinal_cord_PERMUTG2T_40kb TCGAlgg_IDHwt_IDHmutnc"  
#"GSE105318_DLD1_PERMUTG2T_40kb TCGAcoad_msi_mss"  
#"GSE105381_HepG2_PERMUTG2T_40kb TCGAlihc_norm_lihc"  
#"GSE105381_HepG2_PERMUTG2T_40kb TCGAlihc_wt_mutCTNNB1"  
#"LI_PERMUTG2T_40kb TCGAlihc_norm_lihc"  
#"LI_PERMUTG2T_40kb TCGAlihc_wt_mutCTNNB1"  
#"K562_PERMUTG2T_40kb TCGAlaml_wt_mutFLT3"  
#"Panc1_rep12_PERMUTG2T_40kb TCGApaad_wt_mutKRAS" 
#"LG1_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR"
#"LG1_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker"  
#"LG1_PERMUTG2T_40kb TCGAluad_norm_luad"  
#"LG1_PERMUTG2T_40kb TCGAluad_wt_mutKRAS"  
#"LG1_PERMUTG2T_40kb TCGAlusc_norm_lusc"
#"LG2_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR"
#"LG2_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker"  
#"LG2_PERMUTG2T_40kb TCGAluad_norm_luad"  
#"LG2_PERMUTG2T_40kb TCGAluad_wt_mutKRAS"  
#"LG2_PERMUTG2T_40kb TCGAlusc_norm_lusc"
#"PA2_PERMUTG2T_40kb TCGApaad_wt_mutKRAS" 
#"PA3_PERMUTG2T_40kb TCGApaad_wt_mutKRAS" 
#"GSE118588_Panc_beta_PERMUTG2T_40kb TCGApaad_wt_mutKRAS" 
#"Barutcu_MCF-10A_PERMUTG2T_40kb TCGAbrca_lum_bas" 
#"HMEC_PERMUTG2T_40kb TCGAbrca_lum_bas" 
#"GSE109229_BT474_PERMUTG2T_40kb TCGAbrca_lum_bas" 
#"GSE109229_SKBR3_PERMUTG2T_40kb TCGAbrca_lum_bas" 
#"Rao_HCT-116_2017_PERMUTG2T_40kb TCGAcoad_msi_mss" 
#"ENCSR504OTV_transverse_colon_PERMUTG2T_40kb TCGAcoad_msi_mss"
#"GSE99051_786_O_PERMUTG2T_40kb TCGAkich_norm_kich"
#"GSE118514_22Rv1_PERMUTG2T_40kb TCGAprad_norm_prad"
)



#############################################################
for data in "${all_data[@]}"; do
    echo "> START for $data"
	echo Rscript $script_name $data
	Rscript $script_name $data
done




###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

