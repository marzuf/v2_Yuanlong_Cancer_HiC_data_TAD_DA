#!/usr/bin/bash

# ./run_pipeline_permut.sh

# for graviton: 
#conda create -n Renv r-essentials r-base
# conda activate Renv
# if packages needed -> start R after "conda activate Renv" and "install.packages()" as usual




 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_wt_mutBRAF # run pos
 ./run_pipeline.sh ENCSR312KHQ_SK-MEL-5_PERMUTG2T_40kb TCGAskcm_wt_mutCTNNB1 # run pos	
 ./run_pipeline.sh ENCSR862OGI_RPMI-7951_PERMUTG2T_40kb TCGAskcm_wt_mutCTNNB1 # run pos	
# LUNG
 ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_norm_luad   # run pos


 ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR # run pos


# ./run_pipeline.sh ENCSR444WCZ_A549_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker  # 
# ./run_pipeline.sh ENCSR489OCU_NCI-H460_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker # 





# ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR # 
# ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker # 
# ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_norm_luad # 
# ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_wt_mutKRAS # 
# ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAlusc_norm_lusc # 


# ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR # run el
# ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker # run el
# ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_norm_luad # run el
# ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_wt_mutKRAS # run el
# ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAlusc_norm_lusc # run el
# ./run_pipeline.sh PA2_PERMUTG2T_40kb TCGApaad_wt_mutKRAS # run el
# ./run_pipeline.sh PA3_PERMUTG2T_40kb TCGApaad_wt_mutKRAS # run el









