#!/usr/bin/Rscript

# Rscript leDily_quantif_random.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(reshape2)

buildTable <- T

plotCex <- 1.2

thresh1 <- 0.75
thresh2 <- 1-thresh1

geneSignifThresh <- 0.01
tadSignifThresh <- 0.01



all_patterns <- c("RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT", "RANDOMMIDPOS")



for(rd_patt in all_patterns) {
  
  outFolder <- file.path("LEDILY_QUANTIF_RANDOM", rd_patt)
  dir.create(outFolder, recursive = TRUE)
  

  if(buildTable) {
    
    gene_tad_dt <- get(load(file.path("GENE_RANK_TAD_RANK_RANDOM", rd_patt, "all_gene_tad_signif_dt.Rdata")))
    gene_tad_dt$dataset <- file.path(gene_tad_dt$hicds, gene_tad_dt$exprds)
    
    
    final_dt <- get(load("CREATE_FINAL_TABLE_RANDOM/all_result_dt.Rdata"))
    final_dt <- final_dt[ grepl(paste0(rd_patt, "_40kb"), final_dt$hicds), ]
    final_dt$dataset <- file.path(final_dt$hicds, final_dt$exprds)
    
    stopifnot(length(unique(gene_tad_dt$dataset)) == length(unique(final_dt$dataset)))
    
    
    
    final_dt$dataset <- file.path(final_dt$hicds, final_dt$exprds)
    
    final_dt_withRank <- do.call(rbind, by(final_dt, final_dt$dataset, function(x) {
      x$meanLogFC_upRank <- rank(x$meanLogFC)
      x$meanLogFC_downRank <- rank(-x$meanLogFC)
      x }))
    
    gene_tad_dt$geneDElab <- ifelse(gene_tad_dt$adj.P.Val > geneSignifThresh, "notDE_gene",
                                    ifelse(gene_tad_dt$logFC > 0, "upDE_gene", 
                                           ifelse(gene_tad_dt$logFC < 0, "downDE_gene", NA)))
    stopifnot(!is.na(gene_tad_dt$geneDElab))
    
    stopifnot(! (final_dt_withRank$meanLogFC_downRank <= 100 & final_dt_withRank$meanLogFC_upRank <= 100))
    final_dt_withRank$concordanceLab <- ifelse(final_dt_withRank$ratioDown >= thresh1 | final_dt_withRank$ratioDown <= thresh2, 
                                               "concordant_TAD", "discordant_TAD")
    
    final_dt_withRank$concordanceDirLab <- ifelse(final_dt_withRank$ratioDown >= thresh1,  "concordantDown_TAD",
                                                  ifelse(final_dt_withRank$ratioDown <= thresh2,  "concordantUp_TAD", "discordant_TAD"))
    stopifnot(
      which(final_dt_withRank$concordanceDirLab == "discordant_TAD") == which(final_dt_withRank$concordanceLab == "discordant_TAD")
    )
    
    
    final_dt_withRank$activationLab <- ifelse(final_dt_withRank$meanLogFC_upRank <= 100, "activated_TAD", 
                                              ifelse(final_dt_withRank$meanLogFC_downRank <= 100, "repressed_TAD", "nonresp_TAD")) 
    # final_dt_withRank$dirLab <- ifelse(final_dt_withRank$meanLogFC > 0, "up", "down")
    # stopifnot(!is.na(final_dt_withRank$dirLab))
    
    final_dt_withRank$tadDAlab <- ifelse(final_dt_withRank$adjPvalComb > tadSignifThresh, "notDA_TAD",
                                         ifelse(final_dt_withRank$meanLogFC > 0, "upDA_TAD", 
                                                ifelse(final_dt_withRank$meanLogFC < 0, "downDA_TAD", NA)))
    stopifnot(!is.na(final_dt_withRank$tadDAlab))
    
    all_merged_dt <- merge(gene_tad_dt[,c("hicds", "exprds", "region", "entrezID", "geneDElab")], 
                           final_dt_withRank[,c("hicds", "exprds", "region", "activationLab", "concordanceLab", "tadDAlab", "concordanceDirLab")], 
                           by=c("hicds", "exprds", "region"), all = TRUE)
    stopifnot(!is.na(all_merged_dt))
    
    all_merged_dt$dataset <- file.path(all_merged_dt$hicds, all_merged_dt$exprds)
    
    x=all_merged_dt[all_merged_dt$dataset==all_merged_dt$dataset[1],]
    save(x, file="x.Rdata", version=2)
    
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
      nUpDE_genes <- sum(x$geneDElab == "upDE_gene")
      nDownDE_genes <- sum(x$geneDElab == "downDE_gene")
      # % of signif up/down DE genes in activated/repressed TADs
      # => % of upDE_and_activated | downDE_and_repressed
      # upDE_and_activated <- mean(x$geneDElab == "upDE_gene" & x$activationLab == "activated_TAD") 
      # CORRECTED:
      upDE_and_activated <- sum(x$geneDElab == "upDE_gene" & x$activationLab == "activated_TAD")/nUpDE_genes
      # upDE_and_upDA <- mean(x$geneDElab == "upDE_gene" & x$tadDAlab == "upDA_TAD")
      # CORRECTED:
      upDE_and_upDA <- sum(x$geneDElab == "upDE_gene" & x$tadDAlab == "upDA_TAD")/nUpDE_genes
      # downDE_and_repressed <- mean(x$geneDElab == "downDE_gene" & x$activationLab == "repressed_TAD")
      # CORRECTED:
      downDE_and_repressed <- sum(x$geneDElab == "downDE_gene" & x$activationLab == "repressed_TAD")/nDownDE_genes
      # downDE_and_downDA <- mean(x$geneDElab == "downDE_gene" & x$tadDAlab == "downDA_TAD")
      # CORRECTED:
      downDE_and_downDA <- sum(x$geneDElab == "downDE_gene" & x$tadDAlab == "downDA_TAD")/nDownDE_genes
      
      nRepressed_genes <- sum(x$activationLab == "repressed_TAD")
      nActivated_genes <- sum(x$activationLab == "activated_TAD")
      repressed_and_downDE <- sum(x$geneDElab == "downDE_gene" & x$activationLab == "repressed_TAD")/nRepressed_genes
      activated_and_upDE <- sum(x$geneDElab == "upDE_gene" & x$activationLab == "activated_TAD")/nActivated_genes
      
      nDownDA_genes <- sum(x$tadDAlab  == "downDA_TAD")
      nUpDA_genes <- sum(x$tadDAlab == "upDA_TAD")
      downDA_and_downDE <- sum(x$geneDElab == "downDE_gene" & x$tadDAlab  == "downDA_TAD")/nDownDA_genes
      upDA_and_upDE <- sum(x$geneDElab == "upDE_gene" & x$tadDAlab == "upDA_TAD")/nUpDA_genes
      
      upDownDE_and_activRepr <-
        sum(
          (x$geneDElab == "downDE_gene" & x$activationLab == "repressed_TAD") |
            (x$geneDElab == "upDE_gene" & x$activationLab == "activated_TAD")
        )/
        sum(
          (x$geneDElab == "downDE_gene") | 
            (x$geneDElab == "upDE_gene")
        )
      upDEdownDE_and_activrepr <- mean(
        (x$geneDElab == "upDE_gene" & x$activationLab == "activated_TAD") |
          (x$geneDElab == "downDE_gene" & x$activationLab == "repressed_TAD"))
      upDEdownDE_and_upDAdownDA <- mean(
        (x$geneDElab == "upDE_gene" & x$tadDAlab == "upDA_TAD") |
          (x$geneDElab == "downDE_gene" & x$tadDAlab == "downDA_TAD"))
      # if(upDEdownDE_and_upDAdownDA != upDE_and_upDA+downDE_and_downDA) {
      #   cat("ds=",ds,"\n")
      #   cat("upDEdownDE_and_upDAdownDA=\t", upDEdownDE_and_upDAdownDA, "\n")
      #   cat("upDE_and_upDA=\t", upDE_and_upDA, "\n")
      #   cat("downDE_and_downDA=\t", downDE_and_downDA, "\n")  
      # }
      # cat("upDEdownDE_and_activrepr=\t", upDEdownDE_and_activrepr, "\n")
      # cat("upDE_and_activated=\t", upDEdownDE_and_upDAdownDA, "\n")
      # cat("downDE_and_repressed=\t", downDE_and_repressed, "\n")
      # stopifnot(abs(upDEdownDE_and_upDAdownDA- (upDE_and_upDA+downDE_and_downDA)) <= 10^-4)  # not true for CORRECTED
      # stopifnot(abs(upDEdownDE_and_activrepr - (upDE_and_activated+downDE_and_repressed)) <= 10^-4)   # not true for CORRECTED
      
      # % of signif in nonresp TADs
      signifDE_and_notResp <- mean(x$geneDElab != "notDE_gene" & x$activationLab == "nonresp_TAD")
      signifDE_and_notDA <- mean(x$geneDElab != "notDE_gene" & x$tadDAlab == "notDA_TAD")
      
      # C)
      # % of TADs with ratioDown >= 75% and <= 25% ("concordant_TAD")
      # % of DE genes in concordant
      x_t <- x[,c("dataset", "region", "concordanceLab", "activationLab", "concordanceDirLab")]
      x_t <- unique(x_t)
      concordantTAD <- mean(x_t$concordanceLab == "concordant_TAD")
      signifDE_and_concord <- mean(x$geneDElab != "notDE_gene" & x$concordanceLab == "concordant_TAD")
      all_and_concord <- mean(x$concordanceLab == "concordant_TAD")
      
      # Among TAD genes:
      #   => ratio upDE in activated
      # => ratio downDE in repressed
      # => ratio upDE in nonresp.
      # => ratio downDE in nonresp.
      
      ratioGenes_inTAD_dt <- do.call(rbind, by(x, x$region, function(tad) {
        tadlab <- unique(tad$activationLab)
        dalab <- unique(tad$tadDAlab)
        stopifnot(length(tadlab) == 1)
        stopifnot(length(dalab) == 1)
        ratioUpDE <- mean(tad$geneDElab == "upDE_gene")
        ratioDownDE <- mean(tad$geneDElab == "downDE_gene")
        data.frame(tadlab=tadlab,
                   dalab=dalab,
                   ratioUpDE=ratioUpDE,
                   ratioDownDE=ratioDownDE,
                   stringsAsFactors = FALSE)
      }))
      upDE_in_activatedTAD <- mean(ratioGenes_inTAD_dt$ratioUpDE[ratioGenes_inTAD_dt$tadlab=="activated_TAD"])
      downDE_in_repressedTAD <- mean(ratioGenes_inTAD_dt$ratioDownDE[ratioGenes_inTAD_dt$tadlab=="repressed_TAD"])
      
      upDE_in_upDAtad <- mean(ratioGenes_inTAD_dt$ratioUpDE[ratioGenes_inTAD_dt$dalab=="upDA_TAD"])
      downDE_in_downDAtad <- mean(ratioGenes_inTAD_dt$ratioDownDE[ratioGenes_inTAD_dt$dalab=="downDA_TAD"])
      
      
      upDE_in_nonrespTAD <- mean(ratioGenes_inTAD_dt$ratioUpDE[ratioGenes_inTAD_dt$tadlab=="nonresp_TAD"])
      downDE_in_nonrespTAD <- mean(ratioGenes_inTAD_dt$ratioDownDE[ratioGenes_inTAD_dt$tadlab=="nonresp_TAD"])
      
      upDE_in_notDAtad <- mean(ratioGenes_inTAD_dt$ratioUpDE[ratioGenes_inTAD_dt$dalab=="notDA_TAD"])
      downDE_in_notDAtad <- mean(ratioGenes_inTAD_dt$ratioDownDE[ratioGenes_inTAD_dt$dalab=="notDA_TAD"])
      
      
      # Among all TADs:
      #   => % with ratioDown >= 0.75 in repressed and <= 0.25 in activated
      concordantDown_and_repressed <- mean(x_t$concordanceDirLab == "concordantDown_TAD" & x_t$activationLab == "repressed_TAD")
      concordantUp_and_activated <- mean(x_t$concordanceDirLab == "concordantUp_TAD" & x_t$activationLab == "activated_TAD")
      
      
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
        upDEdownDE_and_activrepr=upDEdownDE_and_activrepr,
        upDEdownDE_and_upDAdownDA=upDEdownDE_and_upDAdownDA,
        signifDE_and_notResp=signifDE_and_notResp,
        signifDE_and_notDA=signifDE_and_notDA,
        concordantTAD=concordantTAD,
        signifDE_and_concord=signifDE_and_concord,
        all_and_concord=all_and_concord,
        concordantDown_and_repressed=concordantDown_and_repressed,
        concordantUp_and_activated=concordantUp_and_activated,
        upDE_in_activatedTAD=upDE_in_activatedTAD,
        downDE_in_repressedTAD=downDE_in_repressedTAD,
        upDE_in_nonrespTAD=upDE_in_nonrespTAD,
        downDE_in_nonrespTAD=downDE_in_nonrespTAD,
        repressed_and_downDE=repressed_and_downDE,
        activated_and_upDE=activated_and_upDE,
        downDA_and_downDE=downDA_and_downDE,
        upDA_and_upDE=upDA_and_upDE,
        upDownDE_and_activRepr=upDownDE_and_activRepr,
        upDE_in_upDAtad=upDE_in_upDAtad,
        downDE_in_downDAtad=downDE_in_downDAtad,
        upDE_in_notDAtad=upDE_in_notDAtad,
        downDE_in_notDAtad=downDE_in_notDAtad,
        
        stringsAsFactors=FALSE
      )
    }))
    
    rownames(all_stats_dt) <- NULL
    
    outFile <- file.path(outFolder, "all_stats_dt.Rdata")
    save(all_stats_dt, file=outFile,version=2)
    cat(paste0("... written: ", outFile,"\n"))
    
  } else {
    outFile <- file.path(outFolder, "all_stats_dt.Rdata")
    all_stats_dt <- get(load(outFile))
  }
  
  
  # boxplot(NUMS ~ GRP, data = ddf, lwd = 2, ylab = 'NUMS')
  # stripchart(NUMS ~ GRP, vertical = TRUE, data = ddf, 
  #            method = "jitter", add = TRUE, pch = 20, col = 'blue')
  
  plot_boxplot <- function(plot_dt, ...) {
    plot_dt_m <- melt(plot_dt,id="dataset")
    
    plot_dt_m$variable <- as.character(plot_dt_m$variable)
    plot_dt_m$variable <- gsub("_and_", "_and\n", plot_dt_m$variable)
    
    par(bty="L")
    boxplot(value ~ variable, data = plot_dt_m, lwd = 2, outline=F,
            cex.main=plotCex, cex.lab=plotCex,cex.axis=plotCex,
            xlab="",...)
    stripchart(value ~ variable, vertical = TRUE, data = plot_dt_m, 
               method = "jitter", add = TRUE, pch = 20, col = 'blue')
    mtext(side=3, text = paste0("# datasets = ", length(unique(plot_dt_m$dataset))))
    
  }
  plot_boxplot(all_stats_dt[,c("dataset", "signifDE_and_resp", "all_and_resp")],
               ylab = 'average ratio of genes')
  
  
  plot_density_2_vars <- function(var1, var2, data, saveF=TRUE, ...) {
    stopifnot(var1 %in% colnames(data))
    stopifnot(var2 %in% colnames(data))
    plot_list <- list(data[,var1],data[,var2])
    names(plot_list) <- c(var1, var2)
    if(saveF) outFile <- file.path(outFolder, paste0(var1,"_AND_", var2, ".png"))
    if(saveF) png(outFile, height=350, width=550)
    plot_multiDens(plot_list, ...)
    if(saveF) foo <- dev.off()
  }
  
  plot_list <- list(
    
    c(xvar = "upDEdownDE_and_activrepr",
      yvar = "upDEdownDE_and_upDAdownDA",
      axlab  = "ratio of genes"
    ),
    c(xvar  = "activated_and_upDE",
      yvar = "repressed_and_downDE",
      axlab  ="ratio of activated/repressed genes"
    ),
    c(xvar = "upDA_and_upDE",
      yvar = "downDA_and_downDE",
      axlab = "ratio of up/down DA genes"
    ),
    c(xvar = "signifDE_and_resp",
      yvar = "all_and_resp",
      axlab = "ratio of genes"
    ),
    c(xvar ="signifDE_and_DA",
      yvar = "all_and_DA",
      axlab = "ratio of genes"
    ),
    c(xvar = "upDE_and_activated",
      yvar = "upDE_and_upDA",
      axlab = "ratio of upDE genes"
    ),
    c(xvar = "downDE_and_repressed",
      yvar = "downDE_and_downDA",
      axlab = "ratio of downDE genes"
    ),
    c(xvar = "signifDE_and_notResp",
      yvar = "signifDE_and_notDA",
      axlab = "ratio of genes"
    ),
    c(xvar = "signifDE_and_concord",
      yvar = "all_and_concord",
      axlab = "ratio of genes"
    ),
    c(xvar = "upDE_and_activated",
      yvar = "downDE_and_repressed",
      axlab = "ratio of up/down DE genes"
    ),
    c(xvar = "concordantDown_and_repressed",
      yvar = "concordantUp_and_activated",
      axlab = "ratio of TADs"
    )
  )
  
  for(i in length(plot_list)){
    
    xvar <- plot_list[[i]][["xvar"]]
    yvar <- plot_list[[i]][["yvar"]]
    axlab <- plot_list[[i]][["axlab"]]
    
    outFile <- file.path(outFolder, paste0(xvar, "_AND_", yvar, ".png"))
    png(outFile, height=350, width=550)
    plot_density_2_vars(xvar, yvar,all_stats_dt, my_xlab= axlab, legPos = "topleft")
    foo <- dev.off()
    
    outFile <- file.path(outFolder, paste0(xvar, "_AND_", yvar, "_box.png"))
    png(outFile, height=400, width=400)
    plot_boxplot(all_stats_dt[,c("dataset",xvar, yvar)],
                 ylab = axlab)
    foo <- dev.off()
    cat(paste0("... ", outFile, "\n"))
    
  }
  
 
  
  axlab <- "ratio of genes in TAD"
  outFile <- file.path(outFolder, "geneDE_within_TADs.png")
  png(outFile, height=350, width=550)
  plot_multiDens(
    list(
      upDE_in_activatedTAD=all_stats_dt$upDE_in_activatedTAD,
      downDE_in_repressedTAD=all_stats_dt$downDE_in_repressedTAD,
      upDE_in_nonrespTAD=all_stats_dt$upDE_in_nonrespTAD,
      downDE_in_nonrespTAD=all_stats_dt$downDE_in_nonrespTAD
    ), legPos = "topright", my_xlab=axlab
  )
  foo <- dev.off()
  
  outFile <- file.path(outFolder, "geneDE_within_TADs_box.png")
  png(outFile, height=400, width=400)
  plot_boxplot(all_stats_dt[,c("dataset","upDE_in_activatedTAD", "downDE_in_repressedTAD", "upDE_in_nonrespTAD", "downDE_in_nonrespTAD")],
               ylab = axlab)
  foo <- dev.off()
  
  cat(paste0("... ", outFile, "\n"))
  
  
  single_vars <- c("concordantTAD", 
                   
                   "activated_and_upDE",
                   "repressed_and_downDE",
                   "upDE_in_nonrespTAD",
                   "downDE_in_nonrespTAD",
                   
                   "upDEdownDE_and_activrepr", 
                   "upDownDE_and_activRepr", 
                   
                   "upDE_in_activatedTAD",
                   "upDE_in_upDAtad",
                   "upDE_and_activated", 
                   
                   "downDE_in_repressedTAD",
                   "downDE_in_downDAtad",
                   "downDE_and_repressed",
                   
                   "upDE_and_upDA", 
                   "downDE_and_downDA")
  single_labs <- c("ratio of TADs", 
                   
                   "ratio of upDE genes",
                   "ratio of downDE genes",
                   "ratio of nonresp. genes",
                   "ratio of nonresp. genes",
                   
                   "ratio of activated|repressed genes",
                   "ratio of activated|repressed genes",
                   "ratio of upDE in activated TADs",
                   "ratio of upDE in upDA TADs",
                   "ratio of activated genes", 
                   "ratio of downDE in repressed TADs",
                   "ratio of downDE in downDA TADs",
                   "ratio of repressed genes",
                   "ratio of upDA genes", 
                   "ratio of downDA genes")
  
  stopifnot(length(single_vars) == length(single_labs))
  
  for(i in 1:length(single_vars)) {
    xvar <- single_vars[i]
    axlab <- single_labs[i]
    outFile <- file.path(outFolder, paste0(xvar, ".png"))
    png(outFile, height=350, width=550)
    plot_multiDens(
      list(
        concordantTAD = all_stats_dt[,paste0(xvar)]
      ), legPos = "topright", my_xlab=axlab
    )
    foo <- dev.off()
    
    outFile <- file.path(outFolder, paste0(xvar, "_box.png"))
    png(outFile, height=400, width=400)
    plot_boxplot(all_stats_dt[,c("dataset",xvar)],
                 ylab = axlab)
    mtext(side=1, text=paste0(xvar))
    foo <- dev.off()
    cat(paste0("... ", outFile, "\n"))
  }
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}



