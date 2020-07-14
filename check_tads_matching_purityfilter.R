## => ok
final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
final_dt$regID <- file.path(final_dt$hicds, final_dt$exprds, final_dt$region)
signif_tads <- final_dt$regID[final_dt$adjPvalComb <= 0.01]
disc_tads <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER/EPIC/log10/signifTADs_to_discard.Rdata")

conserved_data <- get(load(file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER/EPIC/log10/conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
conserved_tads <- unlist(conserved_data, use.names = F)

stopifnot(conserved_tads %in% signif_tads)
stopifnot(!conserved_tads %in% disc_tads)

disc_tads <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER//log10/signifTADs_to_discard.Rdata")
conserved_data <- get(load(file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER//log10/conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
conserved_tads <- unlist(conserved_data, use.names = F)

stopifnot(conserved_tads %in% signif_tads)
stopifnot(!conserved_tads %in% disc_tads)