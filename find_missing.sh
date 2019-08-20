#!/usr/bin/bash

#./find_missing.sh 6fastSave permDT


#if [[ $# != 2 ]]; then
#    echo "invalid # of arguments"
#    exit 1
#fi

query_prefix="$1"
query_pattern="$2"




mainFolder="PIPELINE/OUTPUT_FOLDER"



ref_prefix="0"

ref_pattern="ipeline_geneList"



ref_datasets=( $( ls $mainFolder/*/*/$ref_prefix_*/*$ref_pattern*.Rdata ) )

#query_datasets=`dirname ls $mainFolder/*/*/$query_prefix_*/*$query_pattern*.Rdata`


for ref in "${ref_datasets[@]}";do
    #echo $ref
    tmp1=$(dirname "$ref")
    tmp2=$(dirname "$tmp1")
    tmp3=${tmp2##$mainFolder}
    newFile=$(ls $mainFolder/$tmp3/${query_prefix}_*/*$query_pattern.Rdata)

done


