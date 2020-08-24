
# ./mv_randommidpos_to_pd.sh


# all_files=( $( ls PIPELINE/OUTPUT_FOLDER/*RANDOMMIDPOS_*/*/5_*/*data ) ) # >>> DONE to pd2
#all_files=( $( ls PIPELINE/OUTPUT_FOLDER/*RANDOMMIDPOSDISC_*/*/5_*/*data ) ) # >>> DONE to pd2 and pd4
#all_files=( $( ls PIPELINE/OUTPUT_FOLDER/*RANDOMMIDPOSSTRICT_*/*/5_*/*data ) ) # >>> DONEto pd4
#all_files=( $( ls PIPELINE/OUTPUT_FOLDER/*RANDOMMIDPOS_*/*/6_*/*data ) ) # >>> DONE to pd4
#all_files=( $( ls PIPELINE/OUTPUT_FOLDER/*RANDOMMIDPOSDISC_*/*/6_*/*data ) ) # >>> launched to pd4
all_files=( $( ls PIPELINE/OUTPUT_FOLDER/*RANDOMMIDPOSSTRICT_*/*/6_*/*data ) ) # >>> launched to pd4

mv_folder="/mnt/pd4/marie/"  # first pd2, then full so pd4

for file in "${all_files[@]}"; do
    file_dir=`dirname $file`
	echo $file_dir
	new_folder="$mv_folder/$file_dir"
	echo "mkdir -p $new_folder"
	mkdir -p $new_folder
	echo "mv $file $new_folder"
	mv $file $new_folder
done



#mkdir -p /mnt/pd2/marie/PIPELINE/OUTPUT_FOLDER/Rao_HCT-116_2017_RANDOMMIDPOS_40kb/TCGAcoad_msi_mss/5_runPermutationsMedian/

#mv PIPELINE/OUTPUT_FOLDER/Rao_HCT-116_2017_RANDOMMIDPOS_40kb/TCGAcoad_msi_mss/5_runPermutationsMedian/permutationsDT.Rdata /mnt/pd2/marie/PIPELINE/OUTPUT_FOLDER/Rao_HCT-116_2017_RANDOMMIDPOS_40kb/TCGAcoad_msi_mss/5_runPermutationsMedian/


#mkdir -p /mnt/pd2/marie//PIPELINE/OUTPUT_FOLDER/Panc1_rep12_RANDOMMIDPOS_40kb/TCGApaad_wt_mutKRAS/5_runPermutationsMedian
#mv PIPELINE/OUTPUT_FOLDER/Panc1_rep12_RANDOMMIDPOS_40kb/TCGApaad_wt_mutKRAS/5_runPermutationsMedian/permutationsDT.Rdata /mnt/pd2/marie//PIPELINE/OUTPUT_FOLDER/Panc1_rep12_RANDOMMIDPOS_40kb/TCGApaad_wt_mutKRAS/5_runPermutationsMedian

# rsync -avz --remove-source-files -e ssh /this/dir remoteuser@remotehost:/remote/dir 
#        -u, --update                skip files that are newer on the receiver


mv AUC_COEXPRDIST_WITHREG_BOXPLOT /mnt/pd4/marie/ # 
mv AUC_COEXPRDIST_WITHREG_SORTNODUP /mnt/pd4/marie
mv AUC_COEXPRDIST_WITHTF_BOXPLOT  /mnt/pd4/marie
mv AUC_COEXPRDIST_WITHTF_SORTNODUP /mnt/pd4/marie
mv AUC_COEXPRDIST_PARALOGS_BOXPLOT /mnt/pd4/marie # 
mv AUC_COEXPRDIST_PARALOGS_BOXPLOT_PERMUT /mnt/pd4/marie
