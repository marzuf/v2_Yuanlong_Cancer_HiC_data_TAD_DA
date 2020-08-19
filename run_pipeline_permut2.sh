#!/usr/bin/bash

# ./run_pipeline_permut2.sh


# run positron ü

## PROSTATE	
 ./run_pipeline.sh ENCSR346DCU_LNCaP_PERMUTG2T_40kb TCGAprad_norm_prad    
 ./run_pipeline.sh GSE118514_22Rv1_PERMUTG2T_40kb TCGAprad_norm_prad       
 ./run_pipeline.sh GSE118514_RWPE1_PERMUTG2T_40kb TCGAprad_norm_prad       
## GBM
 ./run_pipeline.sh GSE105194_cerebellum_PERMUTG2T_40kb TCGAgbm_classical_mesenchymal 
 ./run_pipeline.sh GSE105194_cerebellum_PERMUTG2T_40kb TCGAgbm_classical_neural 
 ./run_pipeline.sh GSE105194_cerebellum_PERMUTG2T_40kb TCGAgbm_classical_proneural 
 ./run_pipeline.sh GSE105194_cerebellum_PERMUTG2T_40kb TCGAlgg_IDHwt_IDHmutnc 
 ./run_pipeline.sh GSE105194_spinal_cord_PERMUTG2T_40kb TCGAgbm_classical_mesenchymal 
 ./run_pipeline.sh GSE105194_spinal_cord_PERMUTG2T_40kb TCGAgbm_classical_neural  
 ./run_pipeline.sh GSE105194_spinal_cord_PERMUTG2T_40kb TCGAgbm_classical_proneural 
 ./run_pipeline.sh GSE105194_spinal_cord_PERMUTG2T_40kb TCGAlgg_IDHwt_IDHmutnc 
# COLORECTAL
 ./run_pipeline.sh GSE105318_DLD1_PERMUTG2T_40kb TCGAcoad_msi_mss 
 ./run_pipeline.sh ENCSR504OTV_transverse_colon_PERMUTG2T_40kb TCGAcoad_msi_mss 
 ./run_pipeline.sh Rao_HCT-116_2017_PERMUTG2T_40kb TCGAcoad_msi_mss 
# LYMPHOBLAST
 ./run_pipeline.sh K562_PERMUTG2T_40kb TCGAlaml_wt_mutFLT3 
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_norm_luad
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAluad_wt_mutKRAS
 ./run_pipeline.sh LG1_PERMUTG2T_40kb TCGAlusc_norm_lusc
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_mutKRAS_mutEGFR
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_nonsmoker_smoker
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_norm_luad
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAluad_wt_mutKRAS
 ./run_pipeline.sh LG2_PERMUTG2T_40kb TCGAlusc_norm_lusc
 ./run_pipeline.sh PA2_PERMUTG2T_40kb TCGApaad_wt_mutKRAS
 ./run_pipeline.sh PA3_PERMUTG2T_40kb TCGApaad_wt_mutKRAS






