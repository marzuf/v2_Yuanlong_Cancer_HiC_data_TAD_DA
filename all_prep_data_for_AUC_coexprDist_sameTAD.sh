#!/usr/bin/bash

# ./all_prep_data_for_AUC_coexprDist_sameTAD.sh # done

start_time=$(date -R)   

scriptSameTAD="../Yuanlong_Cancer_HiC_data_TAD_DA/create_sameTAD_sortNoDup_otherTADfile.R" 

all_TAD_files_ds=(
#    "ENCSR079VIJ_G401_RANDOMMIDPOS_40kb"
#    "ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOS_40kb" 
#    "ENCSR862OGI_RPMI-7951_RANDOMMIDPOS_40kb" 
#    "ENCSR346DCU_LNCaP_RANDOMMIDPOS_40kb"   
#    "GSE118514_RWPE1_RANDOMMIDPOS_40kb" 
#    "ENCSR401TBQ_Caki2_RANDOMMIDPOS_40kb" 
#    "ENCSR444WCZ_A549_RANDOMMIDPOS_40kb" 
    "ENCSR489OCU_NCI-H460_RANDOMMIDPOSDISC_40kb" 
#    "ENCSR549MGQ_T47D_RANDOMMIDPOS_40kb" 
#    "Barutcu_MCF-7_RANDOMMIDPOS_40kb"  
#    "GSE105194_cerebellum_RANDOMMIDPOS_40kb" 
#    "GSE105194_spinal_cord_RANDOMMIDPOS_40kb" 
#    "GSE105318_DLD1_RANDOMMIDPOS_40kb"   
#    "GSE105381_HepG2_RANDOMMIDPOS_40kb" 
#    "LI_RANDOMMIDPOS_40kb" 
#    "K562_RANDOMMIDPOS_40kb" 
#    "Panc1_rep12_RANDOMMIDPOS_40kb"  
#    "LG1_RANDOMMIDPOS_40kb" 
#    "LG2_RANDOMMIDPOS_40kb" 
#    "PA2_RANDOMMIDPOS_40kb" 
#    "PA3_RANDOMMIDPOS_40kb" 
#    "GSE118588_Panc_beta_RANDOMMIDPOS_40kb" 
#    "Barutcu_MCF-10A_RANDOMMIDPOS_40kb" 
#    "HMEC_RANDOMMIDPOS_40kb" 
#    "GSE109229_BT474_RANDOMMIDPOS_40kb" 
#    "GSE109229_SKBR3_RANDOMMIDPOS_40kb" 
#    "Rao_HCT-116_2017_RANDOMMIDPOS_40kb" 
#    "ENCSR504OTV_transverse_colon_RANDOMMIDPOS_40kb" 
#    "GSE99051_786_O_RANDOMMIDPOS_40kb" 
#    "GSE118514_22Rv1_RANDOMMIDPOS_40kb" 
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

