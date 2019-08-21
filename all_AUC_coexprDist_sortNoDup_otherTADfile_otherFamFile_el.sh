#!/usr/bin/bash

# ./all_AUC_coexprDist_sortNoDup_otherTADfile_otherFamFile_el.sh

start_time=$(date -R)   

script_name="../Yuanlong_Cancer_HiC_data_TAD_DA/AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R"

# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich hgnc
# Rscript AUC_coexprDist_withFam_sortNoDup_otherTADfile_otherFamFile.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 
 
#############################################################
all_data=(
"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf" 
"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF"  
"ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1"  
"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf"
"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF"
"ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1"  
"ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad"  
"GSE118514_RWPE1_40kb TCGAprad_norm_prad"  
"ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich"  
"ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR"  
"ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker"  
"ENCSR444WCZ_A549_40kb TCGAluad_norm_luad"  
"ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS"  
"ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc"  
"ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR"  
"ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker"  
"ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad"  
"ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS"  
"ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc"
"ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas"  
"Barutcu_MCF-7_40kb TCGAbrca_lum_bas" 
"GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal"  
"GSE105194_cerebellum_40kb TCGAgbm_classical_neural" 
"GSE105194_cerebellum_40kb TCGAgbm_classical_proneural"  
"GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc"  
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

