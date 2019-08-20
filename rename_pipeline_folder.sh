#!/usr/bin/bash

# ./rename_pipeline_folder.sh

pipFolder="PIPELINE/OUTPUT_FOLDER"


oldName="7_runPermutationsMeanTADCorr"
newName="7fastSave_runPermutationsMeanTADCorr"


all_folders=( $(realpath $pipFolder/*/*/$oldName) )

for folder in ${all_folders[@]}; do

    echo "-> mv $folder `dirname $folder`/$newName"
    mv $folder `dirname $folder`/$newName

done 
