#!/usr/bin/bash

# ./all_AUC_coexprDist_sortNoDup_otherTADfile_otherFamFile_pos.sh

start_time=$(date -R)   

script_name="../Yuanlong_Cancer_HiC_data_TAD_DA/AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R"

# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich hgnc
# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 
 
#############################################################
all_data=(
"GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal"  
"GSE105194_spinal_cord_40kb TCGAgbm_classical_neural"  
"GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural"  
"GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc"  
"GSE105318_DLD1_40kb TCGAcoad_msi_mss"  
"GSE105381_HepG2_40kb TCGAlihc_norm_lihc"  
"GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1"  
"LI_40kb TCGAlihc_norm_lihc"  
"LI_40kb TCGAlihc_wt_mutCTNNB1"  
"K562_40kb TCGAlaml_wt_mutFLT3"  
"Panc1_rep12_40kb TCGApaad_wt_mutKRAS" 
"LG1_40kb TCGAluad_mutKRAS_mutEGFR"
"LG1_40kb TCGAluad_nonsmoker_smoker"  
"LG1_40kb TCGAluad_norm_luad"  
"LG1_40kb TCGAluad_wt_mutKRAS"  
"LG1_40kb TCGAlusc_norm_lusc"
"LG2_40kb TCGAluad_mutKRAS_mutEGFR"
"LG2_40kb TCGAluad_nonsmoker_smoker"  
"LG2_40kb TCGAluad_norm_luad"  
"LG2_40kb TCGAluad_wt_mutKRAS"  
"LG2_40kb TCGAlusc_norm_lusc"
"PA2_40kb TCGApaad_wt_mutKRAS" 
"PA3_40kb TCGApaad_wt_mutKRAS" 
"GSE118588_Panc_beta_40kb TCGApaad_wt_mutKRAS" 
"Barutcu_MCF-10A_40kb TCGAbrca_lum_bas" 
"HMEC_40kb TCGAbrca_lum_bas" 
"GSE109229_BT474_40kb TCGAbrca_lum_bas" 
"GSE109229_SKBR3_40kb TCGAbrca_lum_bas" 
"Rao_HCT-116_2017_40kb TCGAcoad_msi_mss" 
"ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss"
"GSE99051_786_O_40kb TCGAkich_norm_kich"
"GSE118514_22Rv1_40kb TCGAprad_norm_prad"
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

