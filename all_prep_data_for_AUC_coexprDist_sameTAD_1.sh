#!/usr/bin/bash

# ./all_prep_data_for_AUC_coexprDist_sameTAD.sh # done

start_time=$(date -R)   

scriptSameTAD="../Yuanlong_Cancer_HiC_data_TAD_DA/create_sameTAD_sortNoDup_otherTADfile.R" 

all_TAD_files_ds=(
    "ENCSR079VIJ_G401_PERMUTG2T_40kb"
    "ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb" 
    "ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb" 
    "ENCSR346DCU_LNCaP_PERMUTG2T_40kb"   
    "GSE118514_RWPE1_PERMUTG2T_40kb" 
    "ENCSR401TBQ_Caki2_PERMUTG2T_40kb" 
    "ENCSR444WCZ_A549_PERMUTG2T_40kb" 

    "ENCSR549MGQ_T47D_PERMUTG2T_40kb" 
    "Barutcu_MCF-7_PERMUTG2T_40kb"  
    "GSE105194_cerebellum_PERMUTG2T_40kb" 
#    "GSE105194_spinal_cord_PERMUTG2T_40kb" 
#    "GSE105318_DLD1_PERMUTG2T_40kb"   
#    "GSE105381_HepG2_PERMUTG2T_40kb" 
#    "LI_PERMUTG2T_40kb" 
#    "K562_PERMUTG2T_40kb" 
#    "Panc1_rep12_PERMUTG2T_40kb"  
#    "LG1_PERMUTG2T_40kb" 
#    "LG2_PERMUTG2T_40kb" 
#    "PA2_PERMUTG2T_40kb" 
#    "PA3_PERMUTG2T_40kb" 
#    "GSE118588_Panc_beta_PERMUTG2T_40kb" 
#    "Barutcu_MCF-10A_PERMUTG2T_40kb" 
#    "HMEC_PERMUTG2T_40kb" 
#    "GSE109229_BT474_PERMUTG2T_40kb" 
#    "GSE109229_SKBR3_PERMUTG2T_40kb" 
#    "Rao_HCT-116_2017_PERMUTG2T_40kb" 
#    "ENCSR504OTV_transverse_colon_PERMUTG2T_40kb" 
#    "GSE99051_786_O_PERMUTG2T_40kb" 
#    "GSE118514_22Rv1_PERMUTG2T_40kb" 
)


for ds in "${all_TAD_files_ds[@]}"; do

    echo Rscript $scriptSameTAD $ds
    Rscript $scriptSameTAD $ds

done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

