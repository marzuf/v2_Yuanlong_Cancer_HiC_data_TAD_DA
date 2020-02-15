#!/usr/bin/bash

# ./all_prep_data_for_AUC_coexprDist_dist.sh

start_time=$(date -R)   

scriptPrepDist="../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R"

# LAUNCH MANUALLY BECAUSE VERY SLOW
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R ENCSR079VIJ_G401_40kb -> ok
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R ENCSR312KHQ_SK-MEL-5_40kb -> ok
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R ENCSR862OGI_RPMI-7951_40kb -> ok
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R ENCSR346DCU_LNCaP_40kb -> ok
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE118514_RWPE1_40kb -> ok
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R ENCSR401TBQ_Caki2_40kb -> ok
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R ENCSR444WCZ_A549_40kb -> running
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R ENCSR489OCU_NCI-H460_40kb -> runnin
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R ENCSR549MGQ_T47D_40kb -> run
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R Barutcu_MCF-7_40kb -> running
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE105194_cerebellum_40kb -> runnin neu
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE105194_spinal_cord_40kb -> running neu
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE105318_DLD1_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE105381_HepG2_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R LI_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R K562_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R Panc1_rep12_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R LG1_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R LG2_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R PA2_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R PA3_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE118588_Panc_beta_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R Barutcu_MCF-10A_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R HMEC_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE109229_BT474_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE109229_SKBR3_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R Rao_HCT-116_2017_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R ENCSR504OTV_transverse_colon_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE99051_786_O_40kb
#Rscript ../Yuanlong_Cancer_HiC_data_TAD_DA/create_dist_sortNoDup_otherTADfile.R GSE118514_22Rv1_40kb

 


# DS launched 06.07.2019:
#all_TAD_files_ds=(			
#    "ENCSR079VIJ_G401_40kb" 
#    "ENCSR312KHQ_SK-MEL-5_40kb"  
#    "ENCSR862OGI_RPMI-7951_40kb" 
#    "ENCSR346DCU_LNCaP_40kb"   
#    "GSE118514_RWPE1_40kb" 
#    "ENCSR401TBQ_Caki2_40kb" 
#    "ENCSR444WCZ_A549_40kb" 
#    "ENCSR489OCU_NCI-H460_40kb" 
#    "ENCSR549MGQ_T47D_40kb" 
#    "Barutcu_MCF-7_40kb"  
#    "GSE105194_cerebellum_40kb" 
#    "GSE105194_spinal_cord_40kb" 
#    "GSE105318_DLD1_40kb"   
#    "GSE105381_HepG2_40kb" 
#    "LI_40kb" 
#    "K562_40kb" 
#    "Panc1_rep12_40kb"  
#    "LG1_40kb" 
#    "LG2_40kb" 
#    "PA2_40kb" 
#    "PA3_40kb" 
#    "GSE118588_Panc_beta_40kb" 
#    "Barutcu_MCF-10A_40kb" 
#    "HMEC_40kb" 
#    "GSE109229_BT474_40kb" 
#    "GSE109229_SKBR3_40kb" 
#    "Rao_HCT-116_2017_40kb" 
#    "ENCSR504OTV_transverse_colon_40kb" 
#    "GSE99051_786_O_40kb" 
#    "GSE118514_22Rv1_40kb" 
#)
all_TAD_files_ds=(
#"ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb"
#"ENCSR489OCU_NCI-H460_RANDOMNBRGENES_40kb"
"ENCSR489OCU_NCI-H460_RANDOMSHIFT_40kb"
"ENCSR489OCU_NCI-H460_PERMUTG2T_40kb"
)
for ds in "${all_TAD_files_ds[@]}"; do
    echo Rscript $scriptPrepDist $ds
    Rscript $scriptPrepDist $ds
done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

