#!/usr/bin/bash

clName="K562"

runCMD() {
  echo "> $1"
  eval $1
}
export -f runCMD

maxJobs=100
maxLoad=100
all_chromo=( "chr"{1..5} )  # -> need to do all chromo to have 1 single file at the end

chromo=5

runCMD 'echo $chromo'

sh -c "runCMD 'echo $chromo'"

parallel echo {} ::: ${all_chromo[@]}


#parallel sh -c "runCMD 'echo {}'" ::: "${all_chromo[@]}"

# for all steps
norm="YL"
binSizeKb="40"
Rexec="Rscript"
mainFold="${clName}_${binSizeKb}kb"
maxJobs=200
maxLoad=200
g2t_script="../Cancer_HiC_data_TAD_DA/gene2TAD_consensus_version2.R"
infold_genes="/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"
step1_inputFolder="$mainFold/FINAL_DOMAINS"
step1_outFolder_tmp="$mainFold/genes2tad/tmp"

#parallel -i sh -c "runCMD 'echo hello'"  ::: "${all_chromo[@]}"

#parallel sh -c "runCMD 'Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000'"  ::: "${all_chromo[@]}"

parallel Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000 ::: "${all_chromo[@]}" 



#parallel -i -j $maxJobs -l $maxLoad sh -c "runCMD 'Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000'" -- ${all_chromo[@]}
#parallel -i -j $maxJobs -l $maxLoad sh -c "runCMD 'Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000'" ::: ${all_chromo[@]}
#parallel -i -j $maxJobs -l $maxLoad "echo Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000" ::: ${all_chromo[@]}
#parallel echo {} ::: ${all_chromo[@]}
#parallel -i "echo {}" ::: ${all_chromo[@]}
#parallel -i -j $maxJobs -l $maxLoad "runCMD 'Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000'"  ::: ${all_chromo[@]}
#parallel -i -j $maxJobs -l $maxLoad "Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000"  ::: ${all_chromo[@]}
#parallel sh -c "runCMD 'Rscript $g2t_script -f $infold_genes -t $step1_inputFolder/${clName}_{}_${norm}_${binSizeKb}kb_final_domains.txt -c {} -o $step1_outFolder_tmp -b ${binSizeKb}000'"  -- ${all_chromo[@]}
