#!/usr/bin/bash

# => run all pos


######
# additional datasets (already run 0-4)
######

#>>> lung
./run_pipeline_19fc.sh LG1_40kb TCGAluad_mutKRAS_mutEGFR  # pos ok
./run_pipeline_19fc.sh LG1_40kb TCGAluad_nonsmoker_smoker  # run pos batch
./run_pipeline_19fc.sh LG1_40kb TCGAluad_norm_luad    # run pos batch
./run_pipeline_19fc.sh LG1_40kb TCGAluad_wt_mutKRAS    # run pos batch
./run_pipeline_19fc.sh LG1_40kb TCGAlusc_norm_lusc    # run pos batch

./run_pipeline_19fc.sh LG2_40kb TCGAluad_mutKRAS_mutEGFR  # pos ok
./run_pipeline_19fc.sh LG2_40kb TCGAluad_nonsmoker_smoker   # run pos batch
./run_pipeline_19fc.sh LG2_40kb TCGAluad_norm_luad   # run pos batch
./run_pipeline_19fc.sh LG2_40kb TCGAluad_wt_mutKRAS    # neu pip 1 batch
./run_pipeline_19fc.sh LG2_40kb TCGAlusc_norm_lusc   # neu pip 1 batch
 
#>>> pancreas
./run_pipeline_19fc.sh PA2_40kb TCGApaad_wt_mutKRAS # run pip2 neutrino batch

./run_pipeline_19fc.sh PA3_40kb TCGApaad_wt_mutKRAS # run pip2 neutrino batch

./run_pipeline_19fc.sh GSE118588_Panc_beta_40kb TCGApaad_wt_mutKRAS # run pip2 neutrino batch


#>>> breast
./run_pipeline_19fc.sh Barutcu_MCF-10A_40kb TCGAbrca_lum_bas # run on el batch ok
./run_pipeline_19fc.sh HMEC_40kb TCGAbrca_lum_bas  # run on el batch ok
./run_pipeline_19fc.sh GSE109229_BT474_40kb TCGAbrca_lum_bas  # run on el batch ok 
./run_pipeline_19fc.sh GSE109229_SKBR3_40kb TCGAbrca_lum_bas  # run on el batch ok



#>>> colon
./run_pipeline_19fc.sh Rao_HCT-116_2017_40kb TCGAcoad_msi_mss # run el batch ok

./run_pipeline_19fc.sh ENCSR504OTV_transverse_colon_40kb TCGAcoad_msi_mss # run el batch ok


#>>> kidney
./run_pipeline_19fc.sh GSE99051_786_O_40kb TCGAkich_norm_kich # run on el batch

#>>> prostate 
./run_pipeline_19fc.sh GSE118514_22Rv1_40kb TCGAprad_norm_prad # run on el batch




