
all_analyses = c(
"TCGAkich_norm_kich",
"TCGAskcm_lowInf_highInf", 
"TCGAskcm_wt_mutBRAF",
"TCGAskcm_wt_mutCTNNB1",
"TCGAprad_norm_prad",
"TCGAkich_norm_kich",
"TCGAluad_mutKRAS_mutEGFR", 
"TCGAluad_nonsmoker_smoker",
"TCGAluad_norm_luad",
"TCGAluad_wt_mutKRAS",
"TCGAlusc_norm_lusc",
"TCGAluad_mutKRAS_mutEGFR",
"TCGAluad_nonsmoker_smoker",
"TCGAluad_norm_luad",
"TCGAluad_wt_mutKRAS",
"TCGAlusc_norm_lusc",
"TCGAbrca_lum_bas",
"TCGAskcm_lowInf_highInf",
"TCGAskcm_wt_mutBRAF", 
"TCGAskcm_wt_mutCTNNB1",
"TCGAgbm_classical_mesenchymal",
"TCGAgbm_classical_neural",
"TCGAgbm_classical_proneural",
"TCGAlgg_IDHwt_IDHmutnc",
"TCGAgbm_classical_mesenchymal",
"TCGAgbm_classical_neural",
"TCGAgbm_classical_proneural",
"TCGAlgg_IDHwt_IDHmutnc",
"TCGAcoad_msi_mss",
"TCGAlihc_norm_lihc",
"TCGAlihc_wt_mutCTNNB1",
"TCGAprad_norm_prad",
"TCGAlihc_norm_lihc",
"TCGAlihc_wt_mutCTNNB1",
"TCGAprad_norm_prad",
"TCGAbrca_lum_bas",
"TCGApaad_wt_mutKRAS",
"TCGAlaml_wt_mutFLT3"
)
all_analyses <- unique(all_analyses)

all_cmps = c(
"TCGAkich_norm_kich" = "norm_vs_tumor",
"TCGAskcm_lowInf_highInf" = "subtypes", 
"TCGAskcm_wt_mutBRAF" = "wt_vs_mut",
"TCGAskcm_wt_mutCTNNB1" = "wt_vs_mut",
"TCGAprad_norm_prad" = "norm_vs_tumor",
"TCGAluad_mutKRAS_mutEGFR" = "subtypes", 
"TCGAluad_nonsmoker_smoker" = "subtypes",
"TCGAluad_norm_luad" = "norm_vs_tumor",
"TCGAluad_wt_mutKRAS" = "wt_vs_mut",
"TCGAlusc_norm_lusc" = "norm_vs_tumor",
"TCGAlgg_IDHwt_IDHmutnc" = "wt_vs_mut", 
"TCGAgbm_classical_mesenchymal" = "subtypes",
"TCGAgbm_classical_neural" = "subtypes",
"TCGAgbm_classical_proneural" = "subtypes",
"TCGAcoad_msi_mss" = "subtypes",
"TCGAlihc_norm_lihc" = "norm_vs_tumor",
"TCGAlihc_wt_mutCTNNB1" = "wt_vs_mut",
"TCGAbrca_lum_bas" = "subtypes",
"TCGApaad_wt_mutKRAS" = "wt_vs_mut",
"TCGAlaml_wt_mutFLT3" = "wt_vs_mut"
)

all_cols = c(
"wt_vs_mut" = "red",
"norm_vs_tumor" = "blue",
"subtypes" = "green"
)

addSubtypeLeg <- function(mypos="bottomleft", ...){
  legend(
    mypos,
    col = all_cols,
    legend=names(all_cols),
    lty=1,
    ...
  )
}





