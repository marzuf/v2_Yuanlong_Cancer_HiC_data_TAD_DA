#>>> liver
./run_pipeline_10.sh GSE105381_HepG2_40kb TCGAlihc_norm_lihc  # pos2
./run_pipeline_10.sh GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1  # # pos2
#./run_pipeline_10.sh GSE58752_liver_40kb TCGAlihc_norm_lihc    # update name
#./run_pipeline_10.sh GSE58752_liver_40kb TCGAlihc_wt_mutCTNNB1   # update name  
./run_pipeline_10.sh LI_40kb TCGAlihc_norm_lihc  ## # pos2
./run_pipeline_10.sh LI_40kb TCGAlihc_wt_mutCTNNB1  ## # pos2


#>>> leukemia
./run_pipeline_10.sh K562_40kb TCGAlaml_wt_mutFLT3  # # pos2

#>>> pancreas
./run_pipeline_10.sh Panc1_rep12_40kb TCGApaad_wt_mutKRAS ### pos2 


######
# additional datasets (already run 0-4)
######

#>>> lung
./run_pipeline_10.sh LG1_40kb TCGAluad_mutKRAS_mutEGFR  # pos2
./run_pipeline_10.sh LG1_40kb TCGAluad_nonsmoker_smoker  # pos2
./run_pipeline_10.sh LG1_40kb TCGAluad_norm_luad     # pos2
./run_pipeline_10.sh LG1_40kb TCGAluad_wt_mutKRAS     # pos2
./run_pipeline_10.sh LG1_40kb TCGAlusc_norm_lusc    # pos2
