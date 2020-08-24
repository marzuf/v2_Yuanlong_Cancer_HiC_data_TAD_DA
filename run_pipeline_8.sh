#!/usr/bin/bash

######## #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### for 6v2 and 9v2

# TCGA data have been prepared according to scripts in folder /mnt/etemp/marie/scripts/TAD_DE_pipeline_v2_TCGAdata

# LIVER

# ./run_pipeline_8.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc # running el


#PIPELINE/OUTPUT_FOLDER/GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/5pF_runPermutationsMedian/permutationsDT.Rdata

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



TAD_DE_pipSteps=( "8cOnlyFCCfastSavePF" )

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

