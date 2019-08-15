#!/bin/bash

#./all_cmds_3_assign_genes.sh

# commented names -> names used in Yuanlong_Cancer_HiC_data_TAD_DA

#            all_datasets=(
#            #"ENCSR079VIJ_G401"
#            "ENCSR079VIJ_G401"
#            #"ENCSR105KFX_SK-N-DZ"
#            "ENCSR105KFX_SK-N-DZ"
#            #"ENCSR312KHQ_SK-MEL-5"
#            "ENCSR312KHQ_SK-MEL-5"
#            #"ENCSR346DCU_LNCaP"
#            "ENCSR346DCU_LNCaP"
#            #"ENCSR401TBQ_Caki2"
#            "ENCSR401TBQ_Caki2"
#            #"ENCSR444WCZ_A549"
#            "ENCSR444WCZ_A549"
#            #"ENCSR489OCU_NCI-H460"
#            "ENCSR489OCU_NCI-H460"
#            #"ENCSR549MGQ_T47D"
#            "ENCSR549MGQ_T47D"
#            #"ENCSR834DXR_SK-N-MC"
#            "ENCSR834DXR_SK-N-MC"
#            #"ENCSR862OGI_RPMI-7951"
#            "ENCSR862OGI_RPMI-7951"
#            #"GSE105194_cerebellum"
#            "GSE105194_cerebellum"
#            #"GSE105194_spinal_cord"
#            "GSE105194_spinal_cord"
#            #"GSE105318_DLD1"
#            "GSE105318_DLD1"
#            #"GSE105381_HepG2"
#            "GSE105381_HepG2"
#            #"GSE118514_RWPE1"
#            "GSE118514_RWPE1"
#            #"GSE58752_liver"
#            "LI"
#            #"GSE73782_PC3"  # => bad quality, removed
#            #"GSE75070_MCF-7_shNS"
#            "Barutcu_MCF-7"
#            #"GSM1826481_SK-N-SH"
#            "GSM1826481_SK-N-SH"
#            #"GSM2334832_RPMI-8226_HindIII"
#            "GSM2334832_RPMI-8226_HindIII"
#            #"GSM2334834_U266_MobI"
#            "GSM2334834_U266_MobI"
#            #"Panc1_rep12"
#            "Panc1_rep12 "
#            #"K562"
#            "K562"            
#            )

# additional datasets
all_datasets=(
"LG1"
"LG2"
"PA2"
"PA3"
"GSE118588_Panc_beta"
"Barutcu_MCF-10A"
"HMEC"
"GSE109229_BT474"
"GSE109229_SKBR3"
"Rao_HCT-116_2017"
"ENCSR504OTV_transverse_colon"
"GSE99051_786_O"
"GSE118514_22Rv1"
)



for ds in ${all_datasets[@]}; do

    echo $ds
    ./3_assign_genes.sh $ds

done 


#./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1   						# =>  done 22.01.2019
#parallel -j40 ./3_assign_genes.sh {} ::: ${all_datasets[@]}



