#!/usr/bin/bash

# ./all_prep_data_for_AUC_coexprDist_sameFam.sh

start_time=$(date -R)   

scriptPrepFam="../Yuanlong_Cancer_HiC_data_TAD_DA/prep_gene_families_TAD_data_otherTADfile.R"

scriptSameFam="../Yuanlong_Cancer_HiC_data_TAD_DA/create_sameFamily_sortNoDup_otherFamFile.R"


# DS launched 06.07.2019:
all_TAD_files_ds=(
    "ENCSR079VIJ_G401_RANDOMMIDPOSDISC_40kb"
    "ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSDISC_40kb" 
    "ENCSR862OGI_RPMI-7951_RANDOMMIDPOSDISC_40kb" 
    "ENCSR346DCU_LNCaP_RANDOMMIDPOSDISC_40kb"   
    "GSE118514_RWPE1_RANDOMMIDPOSDISC_40kb" 
    "ENCSR401TBQ_Caki2_RANDOMMIDPOSDISC_40kb" 
    "ENCSR444WCZ_A549_RANDOMMIDPOSDISC_40kb" 
    "ENCSR489OCU_NCI-H460_RANDOMMIDPOSDISC_40kb" 
    "ENCSR549MGQ_T47D_RANDOMMIDPOSDISC_40kb" 
    "Barutcu_MCF-7_RANDOMMIDPOSDISC_40kb"  
    "GSE105194_cerebellum_RANDOMMIDPOSDISC_40kb" 
    "GSE105194_spinal_cord_RANDOMMIDPOSDISC_40kb" 
    "GSE105318_DLD1_RANDOMMIDPOSDISC_40kb"   
    "GSE105381_HepG2_RANDOMMIDPOSDISC_40kb" 
    "LI_RANDOMMIDPOSDISC_40kb" 
    "K562_RANDOMMIDPOSDISC_40kb" 
    "Panc1_rep12_RANDOMMIDPOSDISC_40kb"  
    "LG1_RANDOMMIDPOSDISC_40kb" 
    "LG2_RANDOMMIDPOSDISC_40kb" 
    "PA2_RANDOMMIDPOSDISC_40kb" 
    "PA3_RANDOMMIDPOSDISC_40kb" 
    "GSE118588_Panc_beta_RANDOMMIDPOSDISC_40kb" 
    "Barutcu_MCF-10A_RANDOMMIDPOSDISC_40kb" 
    "HMEC_RANDOMMIDPOSDISC_40kb" 
    "GSE109229_BT474_RANDOMMIDPOSDISC_40kb" 
    "GSE109229_SKBR3_RANDOMMIDPOSDISC_40kb" 
    "Rao_HCT-116_2017_RANDOMMIDPOSDISC_40kb" 
    "ENCSR504OTV_transverse_colon_RANDOMMIDPOSDISC_40kb" 
    "GSE99051_786_O_RANDOMMIDPOSDISC_40kb" 
    "GSE118514_22Rv1_RANDOMMIDPOSDISC_40kb" 
)


for ds in "${all_TAD_files_ds[@]}"; do

    echo Rscript $scriptPrepFam $ds
    Rscript $scriptPrepFam $ds
    
    echo Rscript $scriptSameFam $ds
    Rscript $scriptSameFam $ds
    
done


###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

