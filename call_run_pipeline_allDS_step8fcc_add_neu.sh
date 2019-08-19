
######
# additional datasets (already run 0-4) => run on pip1
######

#>>> lung
./run_pipeline_8fcc.sh LG1_40kb TCGAluad_mutKRAS_mutEGFR  # 
./run_pipeline_8fcc.sh LG1_40kb TCGAluad_nonsmoker_smoker  # 
./run_pipeline_8fcc.sh LG1_40kb TCGAluad_norm_luad    # 
./run_pipeline_8fcc.sh LG1_40kb TCGAluad_wt_mutKRAS    # 
./run_pipeline_8fcc.sh LG1_40kb TCGAlusc_norm_lusc    # 

./run_pipeline_8fcc.sh LG2_40kb TCGAluad_mutKRAS_mutEGFR  # 
./run_pipeline_8fcc.sh LG2_40kb TCGAluad_nonsmoker_smoker   # 
./run_pipeline_8fcc.sh LG2_40kb TCGAluad_norm_luad   # 
./run_pipeline_8fcc.sh LG2_40kb TCGAluad_wt_mutKRAS    # 
./run_pipeline_8fcc.sh LG2_40kb TCGAlusc_norm_lusc   # 
 
#>>> pancreas
./run_pipeline_8fcc.sh PA2_40kb TCGApaad_wt_mutKRAS # 

./run_pipeline_8fcc.sh PA3_40kb TCGApaad_wt_mutKRAS # 

./run_pipeline_8fcc.sh GSE118588_Panc_beta_40kb TCGApaad_wt_mutKRAS # 


#>>> breast
./run_pipeline_8fcc.sh Barutcu_MCF-10A_40kb TCGAbrca_lum_bas #  
./run_pipeline_8fcc.sh HMEC_40kb TCGAbrca_lum_bas  #  
./run_pipeline_8fcc.sh GSE109229_BT474_40kb TCGAbrca_lum_bas  #   
./run_pipeline_8fcc.sh GSE109229_SKBR3_40kb TCGAbrca_lum_bas  #  



#>>> colon
./run_pipeline_8fcc.sh Rao_HCT-116_2017_40kb TCGAcoad_msi_mss # 

./run_pipeline_8fcc.sh ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss # 


#>>> kidney
./run_pipeline_8fcc.sh GSE99051_786_O_40kb TCGAkich_norm_kich # 

#>>> prostate 
./run_pipeline_8fcc.sh GSE118514_22Rv1_40kb TCGAprad_norm_prad # 




