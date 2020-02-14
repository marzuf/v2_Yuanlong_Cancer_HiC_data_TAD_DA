#!/usr/bin/bash

# ./3_assign_genes.sh MCF-7
# ./3_assign_genes.sh ENCSR549MGQ_T47D
# ./3_assign_genes.sh MCF-7ENCSR549MGQ_T47D
# (see what has been run in all_cmds_pip.txt)
# ./3_assign_genes.sh ENCSR489OCU_NCI-H460_RANDOMSHIFT

#all_chromo=( "chr"{1..22} "chrX" )  # -> need to do all chromo to have 1 single file at the end
#parallel -j4 echo {} ::: ${all_chromo[@]}
#parallel -j4 echo {} -- ${all_chromo[@]}
#exit 0

start_time=$(date -R)    
#set -e  # comment otherwise do not cat at the end

#if [[ $# != 2 ]]; then
if [[ $# != 1 ]]; then

    echo "invalid # of arguments"
    exit 1
fi

all_chromo=( "chr"{1..22} "chrX" )  # -> need to do all chromo to have 1 single file at the end
#all_chromo=( "chr15" "chr16" "chr17" )
#all_chromo=( "chr15" )

clName="$1"

echo "*** START ***"
echo "... > Cell line: $clName"
echo "... > Chromosome(s): ${all_chromo[*]}"

step1=1		# assign genes to TADs

###################################**** SOME FUNCTIONS
runCMD() {
  echo "> $1"
  eval $1
}
export -f runCMD

checkFile() {
if [[ ! -f  $2 ]]; then
    echo "... $1 ($2) does not exist !"
    exit 1
fi
}
###################################************************


#####**** HARD-CODED SETTINGS

# for all steps
norm="YL"

binSizeKb="40"
Rexec="Rscript"
mainFold="${clName}_${binSizeKb}kb"

maxJobs=200
maxLoad=200

# step1:
g2t_script="../Cancer_HiC_data_TAD_DA/gene2TAD_consensus_version2.R"
infold_genes="/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"
# HARD CODED IN g2t assignment R SCRIPT (used at the end  of STEP 1 for concatenating files)
g2t_prefix="tmp_g2t"
reg_prefix="tmp_assigned"

step1_inputFolder="$mainFold/FINAL_DOMAINS"

step1_outFolder_tmp="$mainFold/genes2tad/tmp"
step1_outFolder_final="$mainFold/genes2tad"

step1_outFile_genes="$step1_outFolder_final/all_genes_positions.txt"
step1_outFile_tads="$step1_outFolder_final/all_assigned_regions.txt"

#####************************

if [[ $step1 -eq 1 ]]; then

runCMD "mkdir -p $step1_outFolder_final"
runCMD "mkdir -p $step1_outFolder_tmp"

# clean tmp folder (cat command at the end !)
if [[ -d $step1_outFolder_tmp ]]; then
echo "delete folder"
	runCMD "rm -r $step1_outFolder_tmp"
fi
# Panc1_rep12/FINAL_DOMAINS/Panc1_rep12_chr12_YL_40kb_final_domains.txt

parallel -j $maxJobs Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000 ::: "${all_chromo[@]}" 



cat $step1_outFolder_tmp/$g2t_prefix* > $step1_outFile_genes
cat $step1_outFolder_tmp/$reg_prefix* > $step1_outFile_tads

#runCMD "cat $step1_outFolder_tmp/$g2t_prefix* > $step1_outFile_genes"
#runCMD "cat $step1_outFolder_tmp/$reg_prefix* > $step1_outFile_tads"

checkFile step1_outFile_genes $step1_outFile_genes
checkFile step1_outFile_tads $step1_outFile_tads

echo "... Has assigned genes to # chromosomes:"
runCMD "cut -f2 $step1_outFile_genes | sort | uniq | wc -l"

echo "... Has assigned genes to the following chromosomes:"
runCMD "cut -f2 $step1_outFile_genes | sort | uniq | paste -sd\" \""

fi

###################################################################################################################################################
########## END ####################################################################################################################################
###################################################################################################################################################

echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time


exit 0
