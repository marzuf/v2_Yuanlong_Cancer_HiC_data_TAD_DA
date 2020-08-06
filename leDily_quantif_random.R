#!/usr/bin/Rscript

# Rscript leDily_quantif_random.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


all_patterns <- c("RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT", "RANDOMMIDPOS")



for(rd_patt in all_patterns) {
  
  outFolder <- file.path("LEDILY_QUANTIF_RANDOM", rd_patt)
  dir.create(outFolder, recursive = TRUE)
  
  gene_tad_dt <- get(load(file.path("GENE_RANK_TAD_RANK_RANDOM", rd_patt, "all_gene_tad_signif_dt.Rdata")))
  gene_tad_dt$dataset <- file.path(gene_tad_dt$hicds, gene_tad_dt$exprds)
  
  
  final_dt <- get(load("CREATE_FINAL_TABLE_RANDOM/all_result_dt.Rdata"))
  final_dt <- final_dt[ grepl(paste0(rd_patt, "_40kb"), final_dt$hicds), ]
  final_dt$dataset <- file.path(final_dt$hicds, final_dt$exprds)
  
  stopifnot(length(unique(gene_tad_dt$dataset)) == length(unique(final_dt$dataset)))
  
  
  # A)
  # 100 top and bottom TADs by meanLogFC ("activated_TAD", "repressed_TAD" "nonresp_TAD")
  # % of signif. DE genes in activated and repressed
  # % of of all genes in activated and repressed
  
  # B)
  # % of signif up/down DE genes in activated/repressed TADs
  # % of signif in nonresp TADs
  
  # C)
  # % of TADs with ratioDown >= 75% and <= 25% ("concordant_TAD")
  # % of DE genes in concordant
  
  
  final_dt_withRank <- do.call(rbind, by(final_dt, final_dt$dataset, function(x) {
    x$meanLogFC_upRank <- rank(x$meanLogFC)
    x$meanLogFC_downRank <- rank(-x$meanLogFC)
    x }))
  
  thresh1 <- 0.75
  thresh2 <- 1-thresh1
  
  geneSignifThresh <- 0.01
  tadSignifThresh <- 0.01
  
  gene_tad_dt$geneDElab <- ifelse(gene_tad_dt$adj.P.Val > geneSignifThresh, "notDE_gene",
                                  ifelse(gene_tad_dt$logFC > 0, "upDE_gene", 
                                         ifelse(gene_tad_dt$logFC < 0, "downDE_gene", NA)))
  stopifnot(!is.na(gene_tad_dt$geneDElab))
  
  stopifnot(! (final_dt_withRank$meanLogFC_downRank <= 100 & final_dt_withRank$meanLogFC_upRank <= 100))
  final_dt_withRank$concordanceLab <- ifelse(final_dt_withRank$ratioDown >= thresh1 | final_dt_withRank$ratioDown <= thresh2, "concordant_TAD", "discordant_TAD")
  final_dt_withRank$activationLab <- ifelse(final_dt_withRank$meanLogFC_upRank <= 100, "activated_TAD", 
                                            ifelse(final_dt_withRank$meanLogFC_downRank <= 100, "repressed_TAD", "nonresp_TAD")) 
  
  
  final_dt_withRank$tadDAlab <- ifelse(final_dt_withRank$adjPvalComb > tadSignifThresh, "notDA_TAD",
                                       ifelse(final_dt_withRank$meanLogFC > 0, "upDA_TAD", 
                                              ifelse(final_dt_withRank$meanLogFC < 0, "downDA_TAD", NA)))
  stopifnot(!is.na(final_dt_withRank$tadDAlab))
  
  
  all_merged_dt <- merge(gene_tad_dt[,c("hicds", "exprds", "region", "entrezID", "geneDElab")], 
                         final_dt_withRank[,c("hicds", "exprds", "region", "activationLab", "concordanceLab", "tadDAlab")], 
                         by=c("hicds", "exprds", "region"), all = TRUE)
  stopifnot(!is.na(all_merged_dt))
  
  
  all_merged_dt$dataset <- file.path(all_merged_dt$hicds, all_merged_dt$exprds)
  
  x=all_merged_dt[all_merged_dt$dataset==all_merged_dt$dataset[1],]
  
  all_stats_dt <- do.call(rbind, by(all_merged_dt, all_merged_dt$dataset, function(x) {
    
    ds <- unique(x$dataset)
    stopifnot(length(ds) == 1)
    
    # A)
    # 100 top and bottom TADs by meanLogFC ("activated_TAD", "repressed_TAD" "nonresp_TAD")
    # % of signif. DE genes in activated and repressed
    signifDE_and_resp <- mean(x$geneDElab != "notDE_gene" & x$activationLab != "nonresp_TAD")
    # % of of all genes in activated and repressed
    all_and_resp <- mean(x$activationLab != "nonresp_TAD")
    
    # % of signif. DE genes in DA
    signifDE_and_DA <- mean(x$geneDElab != "notDE_gene" & x$tadDAlab != "notDA_TAD")
    # % of of all genes in activated and repressed
    all_and_DA <- mean(x$tadDAlab != "notDA_TAD")
    
    # B)
    # % of signif up/down DE genes in activated/repressed TADs
    upDE_and_activated <- mean(x$geneDElab == "upDE_gene" & x$activationLab == "activated_TAD")
    upDE_and_upDA <- mean(x$geneDElab == "upDE_gene" & x$tadDAlab == "upDA_TAD")
    downDE_and_repressed <- mean(x$geneDElab == "downDE_gene" & x$activationLab == "repressed_TAD")
    downDE_and_downDA <- mean(x$geneDElab == "downDE_gene" & x$tadDAlab == "downDA_TAD")
    
    # % of signif in nonresp TADs
    signifDE_and_notResp <- mean(x$geneDElab != "notDE_gene" & x$activationLab == "nonresp_TAD")
    signifDE_and_notDA <- mean(x$geneDElab != "notDE_gene" & x$tadDAlab == "notDA_TAD")
    
    # C)
    # % of TADs with ratioDown >= 75% and <= 25% ("concordant_TAD")
    # % of DE genes in concordant
    x_t <- x[,c("dataset", "region", "concordanceLab")]
    x_t <- unique(x_t)
    concordantTAD <- mean(x_t$concordanceLab == "concordant_TAD")
    signifDE_and_concord <- mean(x$geneDElab != "notDE_gene" & x$concordanceLab == "concordant_TAD")
    all_and_concord <-  mean( x$concordanceLab == "concordant_TAD")
    data.frame(
      dataset = ds,
      signifDE_and_resp=signifDE_and_resp,
      all_and_resp=all_and_resp,
      signifDE_and_DA=signifDE_and_DA,
      all_and_DA=all_and_DA,
      upDE_and_activated=upDE_and_activated,
      upDE_and_upDA=upDE_and_upDA,
      downDE_and_repressed=downDE_and_repressed,
      downDE_and_downDA=downDE_and_downDA,
      signifDE_and_notResp=signifDE_and_notResp,
      signifDE_and_notDA=signifDE_and_notDA,
      concordantTAD=concordantTAD,
      signifDE_and_concord=signifDE_and_concord,
      all_and_concord=all_and_concord,
      stringsAsFactors=FALSE
    )
    
    
    
    
  }))
  
  rownames(all_stats_dt) <- NULL
  
  outFile <- file.path(outFolder, "all_stats_dt.Rdata")
  save(all_stats_dt, file=outFile,version=2)
  cat(paste0("... written: ", outFile,"\n"))
  
  outFile <- file.path(outFolder, "signifDE_and_resp_AND_all_and_resp.png")
  png(outFile, height=350, width=550)
  plot_multiDens(
    list(
      signifDE_and_resp = all_stats_dt$signifDE_and_resp,
      all_and_resp = all_stats_dt$all_and_resp
    ), legPos = "topleft", my_xlab="ratio of genes"
  )
  foo <- dev.off()
  
  outFile <- file.path(outFolder, "signifDE_and_DA_AND_all_and_DA.png")
  png(outFile, height=350, width=550)
  plot_multiDens(
    list(
      signifDE_and_DA = all_stats_dt$signifDE_and_DA,
      all_and_DA = all_stats_dt$all_and_DA
    ), legPos = "topright", my_xlab="ratio of genes"
  )
  foo <- dev.off()
  
  outFile <- file.path(outFolder, "upDE_and_activated_AND_upDE_and_upDA.png")
  png(outFile, height=350, width=550)
  plot_multiDens(
    list(
      upDE_and_activated = all_stats_dt$upDE_and_activated,
      upDE_and_upDA = all_stats_dt$upDE_and_upDA
    ), legPos = "topright", my_xlab="ratio of genes"
  )
  foo <- dev.off()
  
  outFile <- file.path(outFolder, "downDE_and_repressed_AND_downDE_and_downDA.png")
  png(outFile, height=350, width=550)
  plot_multiDens(
    list(
      downDE_and_repressed = all_stats_dt$downDE_and_repressed,
      downDE_and_downDA = all_stats_dt$downDE_and_downDA
    ), legPos = "topright", my_xlab="ratio of genes"
  )
  foo <- dev.off()
  
  outFile <- file.path(outFolder, "signifDE_and_notResp_AND_signifDE_and_notDA.png")
  png(outFile, height=350, width=550)
  plot_multiDens(
    list(
      signifDE_and_notResp = all_stats_dt$signifDE_and_notResp,
      signifDE_and_notDA = all_stats_dt$signifDE_and_notDA
    ), legPos = "topright", my_xlab="ratio of genes"
  )
  foo <- dev.off()
  
  outFile <- file.path(outFolder, "signifDE_and_concord_AND_all_and_concord.png")
  png(outFile, height=350, width=550)
  plot_multiDens(
    list(
      signifDE_and_concord = all_stats_dt$signifDE_and_concord,
      all_and_concord = all_stats_dt$all_and_concord
    ), legPos = "topleft", my_xlab="ratio of genes"
  )
  foo <- dev.off()
  
  outFile <- file.path(outFolder, "concordantTAD.png")
  png(outFile, height=350, width=550)
  plot_multiDens(
    list(
      concordantTAD = all_stats_dt$concordantTAD
    ), legPos = "topright", my_xlab="ratio of TADs"
  )
  foo <- dev.off()
  
}



