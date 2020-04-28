#!/usr/bin/bash

######## #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### for 6v2 and 9v2

# TCGA data have been prepared according to scripts in folder /mnt/etemp/marie/scripts/TAD_DE_pipeline_v2_TCGAdata

# LIVER
# ./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc  
# ./run_pipeline.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  

# ./run_pipeline.sh LI_40kb TCGAlihc_norm_lihc 
# ./run_pipeline.sh LI_40kb TCGAlihc_wt_mutCTNNB1


# BREAST

# ./run_pipeline.sh ENCSR549MGQ_T47D_40kb TCGAbrca_lum_bas 
# ./run_pipeline.sh Barutcu_MCF-7_40kb TCGAbrca_lum_bas


# KIDNEY
# ./run_pipeline.sh ENCSR079VIJ_G401_40kb TCGAkich_norm_kich 
# ./run_pipeline.sh ENCSR401TBQ_Caki2_40kb TCGAkich_norm_kich


# SKIN
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# ./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf 

# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF
# ./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutBRAF

# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1	
# ./run_pipeline.sh ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1	


# LUNG
# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_norm_luad   
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad   

# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_mutKRAS_mutEGFR 
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR  

# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker 
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker  

# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAluad_wt_mutKRAS 
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_wt_mutKRAS 

# ./run_pipeline.sh ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc  
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc 

	
# PANCREAS
# ./run_pipeline.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS 
	

# PROSTATE	
# ./run_pipeline.sh ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad    
# ./run_pipeline.sh GSE118514_RWPE1_40kb TCGAprad_norm_prad       

# GBM
# ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_mesenchymal 
# ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_neural 
# ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAgbm_classical_proneural 
# ./run_pipeline.sh GSE105194_cerebellum_40kb TCGAlgg_IDHwt_IDHmutnc 

# ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_mesenchymal 
# ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_neural  
# ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAgbm_classical_proneural 
# ./run_pipeline.sh GSE105194_spinal_cord_40kb TCGAlgg_IDHwt_IDHmutnc 


# COLORECTAL
# ./run_pipeline.sh GSE105318_DLD1_40kb TCGAcoad_msi_mss 

# LYMPHOBLAST
# ./run_pipeline.sh K562_40kb TCGAlaml_wt_mutFLT3 


# PERMUT DATA
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_RANDOMSHIFT_40kb TCGAluad_norm_luad
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_norm_luad


# ./run_pipeline.sh ENCSR489OCU_NCI-H460_RANDOMMIDPOS_40kb TCGAluad_norm_luad   

# ./run_pipeline.sh ENCSR489OCU_NCI-H460_RANDOMNBRGENES_40kb TCGAluad_norm_luad   

start_time=$(date -R)    
#set -e

if [[ $# != 2 ]]; then
    echo "invalid # of arguments"
    exit 1
fi

hic_dataset="$1"
expr_dataset="$2"

echo "*** START ***"
echo "... > Hi-C dataset: $hic_dataset"
echo "... > Gene expression dataset: $expr_dataset"

#********************** HARD-CODED SETTINGS FOR THE PIPELINE ********************************************
step1=0     # prepare setting file
step2=1    # run the pipeline

# NB: 1cleanInput is same as 1cleanInputTCGAminCount, except added change rowToKeep based on minCount

#TAD_DE_pipSteps=( "0cleanInputTCGAminCount" "1cleanInputTCGAminCount" "2" "2v2" "3" "4" "5" "6" "7" "8c" "9" "10" "11" "13cleanInput" "14f2" "170revision2EZH2" )
#TAD_DE_pipSteps=( "0cleanInputTCGAminCount" "1cleanInputTCGAminCount" "2" "2v2" "3" "4")
# TAD_DE_pipSteps=( "5" "6" "7" "8c" "9" "10" "11" "13cleanInput" "14f2" "170revision2EZH2" )
#TAD_DE_pipSteps=( "8c" "9" "10" "11" "13cleanInput" "14f2" "170revision2EZH2" )
# done later because slow
#TAD_DE_pipSteps=( "4cond1" "4cond2" )
#TAD_DE_pipSteps=( "9rank" )
#TAD_DE_pipSteps=( "10rank" "11rank" )
# TAD_DE_pipSteps=( "19" )
#TAD_DE_pipSteps=( "5" "6" "9" )
#TAD_DE_pipSteps=( "5" )
#TAD_DE_pipSteps=( "6" "9" "19onlyFC" "8cOnlyRatioDown")
#TAD_DE_pipSteps=( "8cOnlyRatioDown")
#TAD_DE_pipSteps=( "19onlyFC")
#TAD_DE_pipSteps=( "9" "19")
#TAD_DE_pipSteps=( "19")
#TAD_DE_pipSteps=( "10sameNbr" "11sameNbr" )
#TAD_DE_pipSteps=( "8cOnlyRatioDown")
#TAD_DE_pipSteps=( "7" "10")
#TAD_DE_pipSteps=( "19onlyFCandCorr" )
#TAD_DE_pipSteps=( "8cOnlyRatioDown" "19onlyFCandCorr" )
#TAD_DE_pipSteps=( "11sameNbr" )
#TAD_DE_pipSteps=( "11" "11sameNbr" )
#TAD_DE_pipSteps=( "0cleanInputTCGAminCount" )
#TAD_DE_pipSteps=( "1cleanInputTCGAminCount" "3" "4")
#TAD_DE_pipSteps=( "5" "6" "7")
#TAD_DE_pipSteps=( "7")

#TAD_DE_pipSteps=( "8cOnlyRatioDownFastSave" ) ### !!! USE FAST SAVE FOR STEP 8 !!! 
# use fast save and permut version
#TAD_DE_pipSteps=( "5fastSavePermut" "6fastSave")
#TAD_DE_pipSteps=( "11")
#TAD_DE_pipSteps=( "19sameNbrPartial" )
#TAD_DE_pipSteps=( "0cleanInputTCGAminCount" )
#TAD_DE_pipSteps=( "1cleanInputTCGAminCount" "3" "4" )
#TAD_DE_pipSteps=( "5fastSavePermut" "6fastSave")
#TAD_DE_pipSteps=( "5sameNbr" )
#TAD_DE_pipSteps=( "7sameNbr" )
#TAD_DE_pipSteps=( "10sameNbr" )
#TAD_DE_pipSteps=( "8cOnlyFCCfastSave" )
#TAD_DE_pipSteps=( "6fastSave")
#TAD_DE_pipSteps=( "8cOnlyRatioDownFastSave" )
#TAD_DE_pipSteps=( "9" )
#TAD_DE_pipSteps=( "10sameNbr" )
#TAD_DE_pipSteps=( "11sameNbr" )
#TAD_DE_pipSteps=( "19onlyFC" "19sameNbr" )

#TAD_DE_pipSteps=( "9" "10sameNbr" "11sameNbr" "19sameNbr" "19onlyFC" )
#TAD_DE_pipSteps=( "170revision2EZH2" )
# STEP 9; STEP10sameNbr; STEP11sameNbr
#TAD_DE_pipSteps=( "5fastSavePermut" "6fastSave" "5sameNbr" "8cOnlyFCCfastSave" "8cOnlyRatioDownFastSave" "9" "7sameNbr" "10sameNbr" "11sameNbr" "19sameNbr" "19onlyFC" "170revision2EZH2Folders" )
#TAD_DE_pipSteps=( "6fastSave" )
#TAD_DE_pipSteps=( "5sameNbr" "8cOnlyFCCfastSave" "8cOnlyRatioDownFastSave" "9" "7sameNbr" "10sameNbr" "11sameNbr" "19onlyFC" "170revision2EZH2Folders" )
#TAD_DE_pipSteps=( "10sameNbr" "11sameNbr")
#TAD_DE_pipSteps=( "10sameNbr" "11sameNbr" "19sameNbr" "19onlyFC" )

#TAD_DE_pipSteps=( "0cleanInputTCGAminCount")
#TAD_DE_pipSteps=( "0cleanInputTCGAminCount" "1cleanInputTCGAminCount" "3" "4" )
#TAD_DE_pipSteps=(  "5fastSavePermut" "6fastSave" "5sameNbr" "9" "7sameNbr" "10sameNbr" )
#TAD_DE_pipSteps=(  "9sameNbrRescaled" )
#TAD_DE_pipSteps=( "1cleanInputTCGAminCount" "3" "4" )

##TAD_DE_pipSteps=( "5fastSavePermut" "6fastSave" "5sameNbr" "8cOnlyFCCfastSave" "8cOnlyRatioDownFastSave" "9" "7sameNbr" "10sameNbr" "11sameNbr" "170revision2EZH2Folders"  "19sameNbr" "19onlyFC" )
#TAD_DE_pipSteps=( "6fastSave" "5sameNbr" "8cOnlyFCCfastSave" "8cOnlyRatioDownFastSave" "9" "7sameNbr" "10sameNbr" "11sameNbr" "170revision2EZH2Folders"  "19sameNbr" "19onlyFC" )
#TAD_DE_pipSteps=( "6fastSave" "8cOnlyFCCfastSave" "8cOnlyRatioDownFastSave" "9" "7sameNbr" "10sameNbr" "11sameNbr" "170revision2EZH2"  "19sameNbr" "19onlyFC" )
#TAD_DE_pipSteps=( "19onlyFC" "19sameNbr" )
#TAD_DE_pipSteps=( "8cOnlyFCCfastSaveOnlyObs" )

#TAD_DE_pipSteps=( "7sameNbrSpearman" )
#TAD_DE_pipSteps=( "10sameNbrSpearman" "11sameNbrSpearman" )

#TAD_DE_pipSteps=( "5sameNbr" "6sameNbr" "7sameNbr" "8sameNbr" "9sameNbr" "10sameNbr" "11sameNbrSameNbr" )
#TAD_DE_pipSteps=( "9sameNbrDouble" "11sameNbrSameNbrDouble" ) 
#TAD_DE_pipSteps=( "9random" ) 
#TAD_DE_pipSteps=( "9random" "10random" "11random" ) 


# 28.04.2020

#TAD_DE_pipSteps=( "0cleanInputTCGAminCount" "1cleanInputTCGAminCount" "3" "4" "5fastSavePermut" "6fastSave" "5sameNbr" "9" "7sameNbr" "10sameNbr" "11sameNbr" )
#TAD_DE_pipSteps=(  "6fastSave" "5sameNbr" "9" "7sameNbr" "10sameNbr" "11sameNbr" )

# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB20-2_40kb TCGAskcm_lowInf_highInf  # run electron - STOPPED
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB20-3_40kb TCGAskcm_lowInf_highInf  # run electron - STOPPED
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB20-4_40kb TCGAskcm_lowInf_highInf  # run electron - DONE
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB20-5_40kb TCGAskcm_lowInf_highInf  # run electron - run ELECTRON

#TAD_DE_pipSteps=("3" "4" "5fastSavePermut" "6fastSave" "5sameNbr" "9" "7sameNbr" "10sameNbr" "11sameNbr" )

# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB40-2_40kb TCGAskcm_lowInf_highInf  # run positron - DONE
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB40-3_40kb TCGAskcm_lowInf_highInf  # run positron - DONE
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB40-4_40kb TCGAskcm_lowInf_highInf  # run positron - DONE
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB40-5_40kb TCGAskcm_lowInf_highInf  # run positron - DONE

# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB3_40kb TCGAlihc_wt_mutCTNNB1 #  run positron - DONE ?

#TAD_DE_pipSteps=( "0cleanInputTCGAminCount" "1cleanInputTCGAminCount")
TAD_DE_pipSteps=("3" "4" "5fastSavePermut" "6fastSave" "5sameNbr" "9" "7sameNbr" "10sameNbr" "11sameNbr" )

# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB60-2_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB60-3_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB60-4_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB60-5_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB80-2_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB80-3_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB80-4_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB80-5_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB100-2_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB100-3_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB100-4_40kb TCGAskcm_lowInf_highInf  # run 
# ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_RANDOMSUB100-5_40kb TCGAskcm_lowInf_highInf  # run 

# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB1-2_40kb TCGAkich_norm_kich # run POSITRON
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB1-3_40kb TCGAkich_norm_kich # run POS
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB1-4_40kb TCGAkich_norm_kich # run EL
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB1-5_40kb TCGAkich_norm_kich # run EL

# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB2-2_40kb TCGAkich_norm_kich # run 
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB2-3_40kb TCGAkich_norm_kich # run 
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB2-4_40kb TCGAkich_norm_kich # run 
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB2-5_40kb TCGAkich_norm_kich # run 

# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB3-2_40kb TCGAkich_norm_kich # run 
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB3-3_40kb TCGAkich_norm_kich # run 
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB3-4_40kb TCGAkich_norm_kich # run 
# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB3-5_40kb TCGAkich_norm_kich # run 

# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB1-2_40kb TCGAlihc_wt_mutCTNNB1 #  run
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB1-3_40kb TCGAlihc_wt_mutCTNNB1 #  run 
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB1-4_40kb TCGAlihc_wt_mutCTNNB1 #  run 
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB1-5_40kb TCGAlihc_wt_mutCTNNB1 #  run 

# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB2-2_40kb TCGAlihc_wt_mutCTNNB1 #  run 
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB2-3_40kb TCGAlihc_wt_mutCTNNB1 #  run 
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB2-4_40kb TCGAlihc_wt_mutCTNNB1 #  run 
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB2-5_40kb TCGAlihc_wt_mutCTNNB1 #  run 

# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB3-2_40kb TCGAlihc_wt_mutCTNNB1 #  run  
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB3-3_40kb TCGAlihc_wt_mutCTNNB1 #  run 
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB3-4_40kb TCGAlihc_wt_mutCTNNB1 #  run 
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB3-5_40kb TCGAlihc_wt_mutCTNNB1 #  run 

# ./run_pipeline.sh ENCSR079VIJ_G401_RANDOMSUB3-5_40kb TCGAkich_norm_kich

# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB1-5_40kb TCGAlihc_wt_mutCTNNB1
# ./run_pipeline.sh GSE105381_HepG2_RANDOMSUB2-5_40kb TCGAlihc_wt_mutCTNNB1


#TAD_DE_pipSteps=( "0cleanInputTCGAminCountCheckRm" )
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad


runDir="/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA"

TAD_DE_pipDir="/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_TopDom"
TAD_DE_script="./zzz_run_given_step_given_data_v2.sh"
old_inputFolder="/mnt/ed4/marie/scripts/TAD_DE_pipeline/SETTING_FILES_cleanInput"
nPermut="100000"
#nPermut="100"
ncpu="40"

Rexec=`which Rscript`

new_inputFolder="$runDir/PIPELINE/INPUT_FILES/$hic_dataset"
mkdir -p $new_inputFolder
outputFolder="$runDir/PIPELINE/OUTPUT_FOLDER/$hic_dataset/$expr_dataset"
mkdir -p $outputFolder

echo "!!! IMPORTANT HARD-CODED SETTINGS !!!"
echo "... ! step1 = $step1"
echo "... ! step2 = $step2"

echo "... ! old_inputFolder = $old_inputFolder"
echo "... ! nPermut = $nPermut"
echo "... ! ncpu = $ncpu"

if [[ $step2 -eq 1 ]]; then
	echo "... ! TAD_DE_pipDir = $TAD_DE_pipDir"
	echo "... ! TAD_DE_script = $TAD_DE_script"
	echo "... ! TAD_DE_pipStep(s): ${TAD_DE_pipSteps[*]}"
fi

###################################### FUNCTION DEFINITIONS

function mvBack {
  echo "... go back to my folder"
  cd $runDir  
}
trap mvBack EXIT

runCMD() {
  echo "> $1"
  eval $1
}

checkFile() {
if [[ ! -f  $2 ]]; then
    echo "... $1 ($2) does not exist !"
    exit 1
fi
}
############################################################################

old_setting_file="$old_inputFolder/run_settings_${expr_dataset}.R"

checkFile old_setting_file $old_setting_file

new_setting_file="$new_inputFolder/run_settings_${expr_dataset}.R"

if [[ "$step1" -eq 1 ]] ; then
	runCMD "cp $old_setting_file $new_inputFolder"
fi
# !!! TO CHECK FORMAT ZZZZZ
# $hic_dataset/genes2tad/all_assigned_regions.txt
# MCF-7/genes2tad/all_assigned_regions.txt
# TADposDT.txt
#chr1    chr1_TAD1       750001  1300000
new_TADpos_file="$hic_dataset/genes2tad/all_assigned_regions.txt"

# !!! TO CHECK FORMAT ZZZZZ
# $hic_dataset/genes2tad/all_genes_positions.txt
# MCF-7/genes2tad/all_genes_positions.txt
#gene2tadDT.txt
#12893       chr1    761586  762902  chr1_TAD1
new_gene2tad_file="$hic_dataset/genes2tad/all_genes_positions.txt"

checkFile new_TADpos_file $new_TADpos_file

checkFile new_gene2tad_file $new_gene2tad_file

if [[ "$step1" -eq 1 ]] ; then

	cat >> ${new_setting_file} <<- EOM

			# > file edited: `date -R` 

			# path to output folder:
			pipOutFold <- "${outputFolder}"

			# OVERWRITE THE DEFAULT SETTINGS FOR INPUT FILES - use TADs from the current Hi-C dataset 
			TADpos_file <- paste0(setDir, "`realpath $new_TADpos_file`")
							#chr1    chr1_TAD1       750001  1300000
							#chr1    chr1_TAD2       2750001 3650000
							#chr1    chr1_TAD3       3650001 4150000

			gene2tadDT_file <- paste0(setDir, "`realpath $new_gene2tad_file`")
							#LINC00115       chr1    761586  762902  chr1_TAD1
							#FAM41C  chr1    803451  812283  chr1_TAD1
							#SAMD11  chr1    860260  879955  chr1_TAD1
							#NOC2L   chr1    879584  894689  chr1_TAD1

			# overwrite main_settings.R: nCpu <- 25
			nCpu <- ${ncpu}

			# *************************************************************************************************************************
			# ************************************ SETTINGS FOR PERMUTATIONS (5#_, 8c_)
			# *************************************************************************************************************************

			# number of permutations
			nRandomPermut <- $nPermut
			gene2tadAssignMethod <- "maxOverlap"
			nRandomPermutShuffle <- $nPermut
			step8_for_permutGenes <- TRUE
			step8_for_randomTADsFix <- FALSE
			step8_for_randomTADsGaussian <- FALSE
			step8_for_randomTADsShuffle <- FALSE
			step14_for_randomTADsShuffle <- FALSE


            # added here 13.08.2019 to change the filtering of min. read counts
            rm(minCpmRatio)
            min_counts <- 5
            min_sampleRatio <- 0.8


            # to have compatible versions of Rdata
            options(save.defaults = list(version = 2))

	EOM
	echo "WRITTEN: ${new_setting_file}"


echo "> END STEP1:" $(date -R)

checkFile new_setting_file $new_setting_file

fi # end if STEP1

###################################################################################################################################################
### STEP2: MOVE TO THE TAD DE PIPELINE DIRECTORY TO LAUNCH THE TAD DE PIPELINE
###################################################################################################################################################
if [[ "$step2" -eq 1 ]] ; then
	echo "> START STEP2:" $(date -R)
	runCMD "cd $TAD_DE_pipDir"

	echo $TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`
	$TAD_DE_script ${new_setting_file} `echo ${TAD_DE_pipSteps[*]}`

	cd $runDir
	echo "> END STEP2:" $(date -R)
fi

###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################
echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time
exit 0

