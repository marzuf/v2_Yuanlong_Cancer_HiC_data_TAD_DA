options(scipen=100)

SSHFS=F

buildData <- TRUE
# FOCUS OSN TAD !
# Rscript rankDiff_activDiff.R <hicds_norm> <hicds_tumor> <exprds>

# Rscript rankDiff_activDiff.R LI_40kb GSE105381_HepG2_40kb TCGAlihc_norm_lihc #=> ok
# Rscript rankDiff_activDiff.R LG1_40kb ENCSR444WCZ_A549_40kb TCGAluad_norm_luad # => ok
# Rscript rankDiff_activDiff.R LG2_40kb ENCSR444WCZ_A549_40kb TCGAluad_norm_luad # => ok
# Rscript rankDiff_activDiff.R LG1_40kb ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad # => ok
# Rscript rankDiff_activDiff.R LG2_40kb ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad # => ok
# Rscript rankDiff_activDiff.R LG1_40kb ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc #=> ok
# Rscript rankDiff_activDiff.R LG2_40kb ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc #=> ok
# Rscript rankDiff_activDiff.R LG1_40kb ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc # => ok
# Rscript rankDiff_activDiff.R LG2_40kb ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc #
# Rscript rankDiff_activDiff.R GSE118514_RWPE1_40kb ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad # => ok
# Rscript rankDiff_activDiff.R GSE118514_RWPE1_40kb GSE118514_22Rv1_40kb TCGAprad_norm_prad  #=> ok


script_name <- "rankDiff_activDiff.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(foreach)
require(doMC)
registerDoMC(40)


hicds_norm <- "LI_40kb"
hicds_tumor <- "GSE105381_HepG2_40kb"
exprds <- "TCGAlihc_norm_lihc"

hicds_norm <- "LG1_40kb"
hicds_tumor <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAlusc_norm_lusc"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

source("plot_lolliTAD_funct.R")
source("my_heatmap.2.R")
source("annotated_TADs.R"); stopifnot(exists("annotated_tads"))

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.2


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
hicds_norm <- args[1]
hicds_tumor <- args[2]
exprds <- args[3]

setDir=""

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)



inFile <- file.path("TAD_MATCHING_ACROSS_HICDS", "all_matching_dt.Rdata")
matching_dt <- get(load(inFile))

outFolder <- file.path("RANKDIFF_ACTIVDIFF")
dir.create(outFolder, recursive = TRUE)

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))


final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

tadSignifThresh <- 0.01

nPlotted <- 10

toplot <- TRUE

matchingCol <- "matchingID_maxOverlapBp"


######################### PREPARE NORM DATA

### prepare TAD rank data for norm hicds

all_norm_files <- list.files(file.path(hicds_norm, "FINAL_DOMAINS_WITH_SCORES"), pattern="final_domains_with_scores.txt$", full.names = TRUE)
# aggregate the rank values
allChr_norm_rankDT <- foreach(normFile = all_norm_files, .combine='rbind') %dopar% {
  dt <- read.delim(normFile, header=F, col.names=c("chromo", "start", "end", "rankValue"))
  dt
}
norm_assigned_region_file <- file.path(hicds_norm, "genes2tad", "all_assigned_regions.txt")
norm_assigned_region_DT <- read.delim(norm_assigned_region_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))

norm_regionFile <- file.path(pipFolder, hicds_norm, exprds, script0_name, "pipeline_regionList.Rdata")
norm_regionList <- get(load(norm_regionFile))

norm_assigned_region_DT <- norm_assigned_region_DT[norm_assigned_region_DT$region %in% norm_regionList,]
stopifnot(setequal(norm_regionList, norm_assigned_region_DT$region))

norm_assigned_region_withRank_DT <- foreach(i=1:nrow(norm_assigned_region_DT),.combine='rbind')%dopar% {
  # norm_assigned_region_DT[i,]
  curr_chromo <- norm_assigned_region_DT$chromo[i]
  curr_region <- norm_assigned_region_DT$region[i]
  curr_TADstart <- norm_assigned_region_DT$start[i]
  curr_TADend <- norm_assigned_region_DT$end[i]
  kpIdx <- which(allChr_norm_rankDT$chromo == curr_chromo &
                   (allChr_norm_rankDT$start == curr_TADstart | allChr_norm_rankDT$end == curr_TADend ))
  stopifnot(length(kpIdx) > 0)
  stopifnot(length(kpIdx) <= 2)
  curr_TADvalue <- mean(allChr_norm_rankDT$rankValue[kpIdx])
  
  curr_meanFC <- final_table_DT$meanLogFC[final_table_DT$hicds == hicds_norm & final_table_DT$exprds == exprds & final_table_DT$region == curr_region]
  stopifnot(length(curr_meanFC) == 1)
  stopifnot(!is.na(curr_meanFC) )
  
  data.frame(
    chromo=curr_chromo,
    region=curr_region,
    start=curr_TADstart,
    end=curr_TADend,
    regionRank = curr_TADvalue,
    normMeanFC = curr_meanFC,
    stringsAsFactors = FALSE
  )
}

### prepare the TAD pvalues for norm hicds
norm_tadPval_file <- file.path(pipFolder, hicds_norm, exprds, script11same_name, "emp_pval_combined.Rdata" )
norm_TAD_Pvals <- get(load(norm_tadPval_file))
norm_TAD_adjPvals <- p.adjust(norm_TAD_Pvals, method="BH")

norm_TAD_adjPvals_dt <- data.frame(refID = names(norm_TAD_adjPvals), adjPval=as.numeric(norm_TAD_adjPvals), stringsAsFactors = FALSE)


### prepare the TAD matching for norm hicds
matching_dt$ref_hicds <- dirname(matching_dt$ref_dataset)
matching_dt$ref_exprds <- basename(matching_dt$ref_dataset)
matching_dt$matching_hicds <- dirname(matching_dt$matching_dataset)
matching_dt$matching_exprds <- basename(matching_dt$matching_dataset)


norm_matching_dt <- matching_dt[matching_dt$ref_exprds == exprds & 
                                  matching_dt$ref_hicds == hicds_norm & 
                                  matching_dt$matching_exprds == exprds &
                                  matching_dt$matching_hicds == hicds_tumor,]

norm_matching_pval_dt <- merge(norm_matching_dt, norm_TAD_adjPvals_dt, by="refID")


norm_matching_pval_dt <- norm_matching_pval_dt[!is.na(norm_matching_pval_dt[,paste0(matchingCol)]),]
norm_matching_pval_dt <- norm_matching_pval_dt[order(norm_matching_pval_dt[, paste0(matchingCol)], decreasing = TRUE),]
# nrow(norm_matching_pval_dt)
stopifnot(!duplicated(norm_matching_pval_dt$refID))

######################### PREPARE TUMOR DATA

### prepare TAD rank data for tumor hicds
all_tumor_files <- list.files(file.path(hicds_tumor, "FINAL_DOMAINS_WITH_SCORES"), pattern="final_domains_with_scores.txt$", full.names = TRUE)
# aggregate the rank values
allChr_tumor_rankDT <- foreach(tumorFile = all_tumor_files, .combine='rbind') %dopar% {
  dt <- read.delim(tumorFile, header=F, col.names=c("chromo", "start", "end", "rankValue"))
  dt
}

tumor_assigned_region_file <- file.path(hicds_tumor, "genes2tad", "all_assigned_regions.txt")
tumor_assigned_region_DT <- read.delim(tumor_assigned_region_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))

tumor_regionFile <- file.path(pipFolder, hicds_tumor, exprds, script0_name, "pipeline_regionList.Rdata")
tumor_regionList <- get(load(tumor_regionFile))

tumor_assigned_region_DT <- tumor_assigned_region_DT[tumor_assigned_region_DT$region %in% tumor_regionList,]
stopifnot(setequal(tumor_regionList, tumor_assigned_region_DT$region))

tumor_assigned_region_withRank_DT <- foreach(i=1:nrow(tumor_assigned_region_DT),.combine='rbind')%dopar% {
  # tumor_assigned_region_DT[i,]
  curr_chromo <- tumor_assigned_region_DT$chromo[i]
  curr_region <- tumor_assigned_region_DT$region[i]
  curr_TADstart <- tumor_assigned_region_DT$start[i]
  curr_TADend <- tumor_assigned_region_DT$end[i]
  kpIdx <- which(allChr_tumor_rankDT$chromo == curr_chromo &
                   (allChr_tumor_rankDT$start == curr_TADstart | allChr_tumor_rankDT$end == curr_TADend ))
  stopifnot(length(kpIdx) > 0)
  stopifnot(length(kpIdx) <= 2)
  curr_TADvalue <- mean(allChr_tumor_rankDT$rankValue[kpIdx])
  
  
  curr_meanFC <- final_table_DT$meanLogFC[final_table_DT$hicds == hicds_tumor & final_table_DT$exprds == exprds & final_table_DT$region == curr_region]
  stopifnot(length(curr_meanFC) == 1)
  stopifnot(!is.na(curr_meanFC) )
  
  
  data.frame(
    chromo=curr_chromo,
    region=curr_region,
    start=curr_TADstart,
    end=curr_TADend,
    regionRank = curr_TADvalue,
    tumorMeanFC = curr_meanFC,
    stringsAsFactors = FALSE
  )
}



### prepare the TAD pvalues for tumor hicds
tumor_tadPval_file <- file.path(pipFolder, hicds_tumor, exprds, script11same_name, "emp_pval_combined.Rdata" )
tumor_TAD_Pvals <- get(load(tumor_tadPval_file))
tumor_TAD_adjPvals <- p.adjust(tumor_TAD_Pvals, method="BH")
tumor_TAD_adjPvals_dt <- data.frame(refID = names(tumor_TAD_adjPvals), adjPval=as.numeric(tumor_TAD_adjPvals), stringsAsFactors = FALSE)

### prepare the TAD matching for tumor hicds
tumor_matching_dt <- matching_dt[matching_dt$ref_exprds == exprds & 
                                  matching_dt$ref_hicds == hicds_tumor & 
                                  matching_dt$matching_exprds == exprds &
                                  matching_dt$matching_hicds == hicds_norm,]


tumor_matching_pval_dt <- merge(tumor_matching_dt, tumor_TAD_adjPvals_dt, by="refID")

tumor_matching_pval_dt <- tumor_matching_pval_dt[!is.na(tumor_matching_pval_dt[,paste0(matchingCol)]),]
tumor_matching_pval_dt <- tumor_matching_pval_dt[order(tumor_matching_pval_dt[, paste0(matchingCol)], decreasing = TRUE),]
# nrow(tumor_matching_pval_dt)
stopifnot(!duplicated(tumor_matching_pval_dt$refID))



###### MERGE RANK AND PVALS - norm
norm_matching_pval_dt <- norm_matching_pval_dt[,c("ref_hicds", "ref_exprds", "matching_hicds", "matching_exprds", "refID", matchingCol, "adjPval")]
colnames(norm_assigned_region_withRank_DT)[colnames(norm_assigned_region_withRank_DT)=="region"] <- "refID"
norm_matching_pval_tadRank_dt <- merge(norm_matching_pval_dt, norm_assigned_region_withRank_DT, by=c("refID"), all.x=TRUE, all.y=FALSE)
colnames(norm_matching_pval_tadRank_dt)[colnames(norm_matching_pval_tadRank_dt) == "regionRank"] <- "refID_rank"

colnames(tumor_assigned_region_withRank_DT)[colnames(tumor_assigned_region_withRank_DT)=="region"] <- paste0(matchingCol)
norm_matching_pval_tadRank_dt <- merge(norm_matching_pval_tadRank_dt, tumor_assigned_region_withRank_DT, by=c(paste0(matchingCol)), all.x=TRUE, all.y=FALSE)
colnames(norm_matching_pval_tadRank_dt)[colnames(norm_matching_pval_tadRank_dt) == "regionRank"] <- "matchingID_rank"
norm_matching_pval_tadRank_dt <- unique(norm_matching_pval_tadRank_dt)

norm_matching_pval_tadRank_dt$rankDiff <- norm_matching_pval_tadRank_dt$refID_rank - norm_matching_pval_tadRank_dt$matchingID_rank


###### MERGE RANK AND PVALS - tumor
tumor_matching_pval_dt <- tumor_matching_pval_dt[,c("ref_hicds", "ref_exprds", "matching_hicds", "matching_exprds", "refID", matchingCol, "adjPval")]
colnames(tumor_assigned_region_withRank_DT)[colnames(tumor_assigned_region_withRank_DT)==paste0(matchingCol)] <- "refID" # was matchingcol because of above
tumor_matching_pval_tadRank_dt <- merge(tumor_matching_pval_dt, tumor_assigned_region_withRank_DT, by=c("refID"), all.x=TRUE, all.y=FALSE)
colnames(tumor_matching_pval_tadRank_dt)[colnames(tumor_matching_pval_tadRank_dt) == "regionRank"] <- "refID_rank"

colnames(norm_assigned_region_withRank_DT)[colnames(norm_assigned_region_withRank_DT)=="refID"] <- paste0(matchingCol) # was refID because of above
tumor_matching_pval_tadRank_dt <- merge(tumor_matching_pval_tadRank_dt, norm_assigned_region_withRank_DT, by=c(paste0(matchingCol)), all.x=TRUE, all.y=FALSE)
colnames(tumor_matching_pval_tadRank_dt)[colnames(tumor_matching_pval_tadRank_dt) == "regionRank"] <- "matchingID_rank"
tumor_matching_pval_tadRank_dt <- unique(tumor_matching_pval_tadRank_dt)

tumor_matching_pval_tadRank_dt$rankDiff <- tumor_matching_pval_tadRank_dt$refID_rank - tumor_matching_pval_tadRank_dt$matchingID_rank

outFile <- file.path(outFolder, "tumor_matching_pval_tadRank_dt.Rdata")
save(tumor_matching_pval_tadRank_dt, file=outFile, version = 2)
cat(paste0("... written: ", outFile, "\n"))
# tumor_matching_pval_tadRank_dt <- get(load(outFile))

outFile <- file.path(outFolder, "norm_matching_pval_tadRank_dt.Rdata")
save(norm_matching_pval_tadRank_dt, file=outFile, version = 2)
cat(paste0("... written: ", outFile, "\n"))
# norm_matching_pval_tadRank_dt <- get(load(outFile))

all_refs <- c("norm", "tumor")
curr_ref <- "norm"

for(curr_ref in all_refs) {
  
  curr_hicds <- eval(parse(text = paste0("hicds_", curr_ref)))
  curr_match <- all_refs[all_refs != curr_ref]
  stopifnot(length(curr_match) == 1)
  match_hicds <- eval(parse(text = paste0("hicds_", curr_match)))
  
  curr_dt <- eval(parse(text=paste0(curr_ref, "_matching_pval_tadRank_dt")))
  
  curr_dt$textCol <- ifelse( file.path(curr_dt$ref_hicds, curr_dt$ref_exprds, curr_dt$refID) %in% annotated_tads, "red", "black" )
  
  myx <- -log10(curr_dt$adjPval)
  myy <- curr_dt$rankDiff
  
  if(toplot) {
    
    # outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_rankDiff_vs_pval_", exprds, "_densplot.", plotType ))
    # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    # densplot(
    #   x = myx,
    #   y = myy,
    #   main=paste0(exprds),
    #   sub=paste0(curr_ref, " as refDS"),
    #   xlab=paste0("-log10 TAD adj. pval"),
    #   ylab=paste0("best matching TAD rank diff. (", curr_ref, "-", curr_match, ")"),
    #   cex.axis=axisCex,
    #   cex.lab=axisCex
    # )
    # addCorr(x = myx, y=myy, bty="n")
    # abline(h=0, lty=2, col="grey")
    # mtext(side=3, paste0(curr_hicds, " matching ", match_hicds))
    # foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
  

  
  
  myx <- curr_dt[, paste0(curr_ref, "MeanFC")]
  
  if(toplot) {
  
    # outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_", curr_ref, "MeanFC_vs_pval_", exprds, "_densplot.", plotType ))
    # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    # densplot(
    #   x = myx,
    #   y = myy,
    #   main=paste0(exprds),
    #   sub=paste0(curr_ref, " as refDS"),
    #   xlab=paste0("mean TAD logFC"),
    #   ylab=paste0("best matching TAD rank diff. (", curr_ref, " - ", curr_match, ")"),
    #   cex.axis=axisCex,
    #   cex.lab=axisCex
    # )
    # addCorr(x = myx, y=myy, bty="n")
    # abline(h=0, lty=2, col="grey")
    # mtext(side=3, paste0(curr_hicds, " matching ", match_hicds))
    # foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  
  }
  
  ###################################################################################### plot only signif TADs  
  
  signif_dt <- curr_dt[curr_dt$adjPval <= tadSignifThresh,]
  
  signif_dt <- signif_dt[order(abs(signif_dt$rankDiff), decreasing = TRUE),]
  
  
  outFile <- file.path(outFolder, paste0( hicds_norm, "_", hicds_tumor, "_", exprds, "_ref_", curr_ref, ".Rdata"))
  save(signif_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  myx <- -log10(signif_dt$adjPval)
  myy <- signif_dt$rankDiff
  
  if(toplot) {
    outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_rankDiff_vs_pval_", exprds, "_densplot_signifTADsOnly.", plotType ))
    # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    # densplot(
    #   x = myx,
    #   y = myy,
    #   main=paste0(exprds),
    #   sub=paste0(curr_ref, " as refDS (tadPval<=", tadSignifThresh, ")" ),
    #   xlab=paste0("-log10 TAD adj. pval"),
    #   ylab=paste0("best matching TAD rank diff. (", curr_ref, "-", curr_match, ")"),
    #   cex.axis=axisCex,
    #   cex.lab=axisCex
    # )
    # # add text label for the top rankDiff
    # text(x=myx[1:nPlotted], y=myy[1:nPlotted], labels = signif_dt$refID[1:nPlotted], cex=0.6, col = signif_dt$textCol[1:nPlotted], pos=3)
    # 
    # 
    # addCorr(x = myx, y=myy, bty="n")
    # abline(h=0, lty=2, col="grey")
    # mtext(side=3, paste0(curr_hicds, " matching ", match_hicds))
    # # add text label for the top pvals
    # signif_dt_sortPval <- signif_dt[order(signif_dt$adjPval),]
    # myx_sortPval <- -log10(signif_dt_sortPval$adjPval)
    # myy_sortPval <- signif_dt_sortPval$rankDiff
    # signif_dt_sortPval$textCol <- ifelse( file.path(signif_dt_sortPval$ref_hicds, signif_dt_sortPval$ref_exprds, signif_dt_sortPval$refID) %in% annotated_tads, "red", "darkgrey" )
    # text(x=myx_sortPval[1:nPlotted], y=myy_sortPval[1:nPlotted], labels = signif_dt_sortPval$refID[1:nPlotted], cex=0.6, col = signif_dt_sortPval$textCol[1:nPlotted], pos=3)
    # 
    # foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
  }
  

  
  
  myx <- signif_dt[, paste0(curr_ref, "MeanFC")]
  
  if(toplot){ 
    outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_", curr_ref, "MeanFC_vs_pval_", exprds, "_densplot_siginfTADsOnly.", plotType ))
    # do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    # densplot(
    #   x = myx,
    #   y = myy,
    #   main=paste0(exprds),
    #   sub=paste0(curr_ref, " as refDS (tadPval<=", tadSignifThresh, ")" ),
    #   xlab=paste0("mean TAD logFC"),
    #   ylab=paste0("best matching TAD rank diff. (", curr_ref, " - ", curr_match, ")"),
    #   cex.axis=axisCex,
    #   cex.lab=axisCex
    # )
    # text(x=myx[1:nPlotted], y=myy[1:nPlotted], labels = signif_dt$refID[1:nPlotted], cex=0.6, col = signif_dt$textCol[1:nPlotted])
    # addCorr(x = myx, y=myy, bty="n")
    # abline(h=0, lty=2, col="grey")
    # mtext(side=3, paste0(curr_hicds, " matching ", match_hicds))
    # 
    # # add text label for the top pvals
    # signif_dt_sortPval <- signif_dt[order(signif_dt$adjPval),]
    # myx_sortPval <- signif_dt_sortPval[, paste0(curr_ref, "MeanFC")]
    # myy_sortPval <- signif_dt_sortPval$rankDiff
    # signif_dt_sortPval$textCol <- ifelse( file.path(signif_dt_sortPval$ref_hicds, signif_dt_sortPval$ref_exprds, signif_dt_sortPval$refID) %in% annotated_tads, "red", "darkgrey" )
    # text(x=myx_sortPval[1:nPlotted], y=myy_sortPval[1:nPlotted], labels = signif_dt_sortPval$refID[1:nPlotted], cex=0.6, col = signif_dt_sortPval$textCol[1:nPlotted], pos=3)
    # 
    
    # foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  

  
  # save the file
  out_signif_dt <- signif_dt[order(abs(signif_dt$rankDiff), decreasing=TRUE),]
  

  outCols <- c("ref_hicds", "matching_hicds", "ref_exprds", "refID", "rankDiff", "adjPval")
  
  out_signif_dt <- out_signif_dt[,outCols]
  
  out_signif_dt$rankDiff <- round(out_signif_dt$rankDiff, 4)
  out_signif_dt$adjPval <- round(out_signif_dt$adjPval, 4)
  
  out_signif_dt$adjPval_rank <- rank(out_signif_dt$adjPval, ties="min")
  
  outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_", exprds, "_signifDT.txt" ))
  write.table(out_signif_dt, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, append=F, file =outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  plotList <- list()
  
  toplot_tads <- out_signif_dt$refID[1:nPlotted]
  
  
  if(toplot) {
    
    pipOutFolder <- file.path("PIPELINE/OUTPUT_FOLDER/")
    # ADDED -> load 1x/dataset
    rnaseqDTFile <- file.path(pipOutFolder, curr_hicds, exprds,script0_name, "rna_rnaseqDT.Rdata")
    initListFile <- file.path(pipOutFolder, curr_hicds, exprds,script0_name, "rna_geneList.Rdata")
    geneListFile <- file.path(pipOutFolder, curr_hicds, exprds,script0_name, "pipeline_geneList.Rdata")
    DE_topTableFile <- file.path(pipOutFolder, curr_hicds, exprds, script1_name, "DE_topTable.Rdata")
    entrezDT_file <- file.path("/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
    gene2tadDT_file <- file.path(curr_hicds, "genes2tad", "all_genes_positions.txt")
    cat(paste0("... load gene2tad_dt\n"))
    curr_gene2tadDT <-  read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
    cat(paste0("... load entrez2symb_dt\n"))
    curr_entrez2symbDT<- entrez2symb_dt
    cat(paste0("... load DE_topTable\n"))
    curr_DE_topTable <- eval(parse(text = load(DE_topTableFile))) 
    cat(paste0("... load rnaseqDT\n"))
    curr_rnaseqDT <- eval(parse(text = load(rnaseqDTFile))) 
    cat(paste0("... load initList\n"))
    curr_initList <- eval(parse(text = load(initListFile)))	
    cat(paste0("... load geneList\n"))
    curr_geneList <- eval(parse(text = load(geneListFile)))      
    
    ############################################# LOLLIPLOT FOR THE TOP RANK DIFF
    # foo <- foreach(i_tad = 1:nPlotted) %dopar% {
    for(i_tad in 1:nPlotted) {
      tad <- toplot_tads[i_tad]
      mytit <- paste0( curr_hicds, " - ", exprds, " - ", tad, "\n", "(# ", i_tad, "; rankDiff=", round(out_signif_dt$rankDiff[i_tad], 2)," (adj. pval rank: ", out_signif_dt$adjPval_rank[i_tad], "/", max(out_signif_dt$adjPval_rank), ")")
      
      
      
      
      
      

      
      
      
      
      
      
      plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                            hicds = curr_hicds,
                                            all_TADs = tad,
                                            
                                            
                                            gene2tadDT=curr_gene2tadDT,entrez2symbDT=curr_entrez2symbDT,
                                            DE_topTable=curr_DE_topTable,
                                            rnaseqDT=curr_rnaseqDT, initList=curr_initList, geneList=curr_geneList,
                                            
                                            
                                            orderByLolli = "startPos", mytitle=mytit)
      
      nPlotted_s <- nPlotted
      nPlotted <- 1
      outFile <- file.path(outFolder, paste0(i_tad, "_", curr_hicds, "_withMatching_", match_hicds, "_", exprds, "_rankDiff_top", nPlotted, "_lolli.", plotType))
      # mytit <- paste0("Top ", nPlotted, " rankDiff signif. TADs (<=",tadSignifThresh , ") - ", curr_hicds, " matching ", match_hicds, " (", exprds, ")")
      mytit <- paste0("Top # ", i_tad, "_", nPlotted, " rankDiff signif. TADs (<=",tadSignifThresh , ") - ", curr_hicds, " matching ", match_hicds, " (", exprds, ")")
      # all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      all_plots <- do.call(grid.arrange, c(plotList[[i_tad]],  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      outHeightGG <- min(c(7 * nPlotted/2, 49))
      outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
      outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
      # ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
      ggsave(filename = outFile, plotList[[i_tad]], width=10, height = 7)
      cat("... written: ", outFile, "\n")
      nPlotted <- nPlotted_s
      
    } # end-for iterating over TADs to plot
    
    # outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_", exprds, "_rankDiff_top", nPlotted, "_lolli.", plotType))
    # 
    # mytit <- paste0("Top ", nPlotted, " rankDiff signif. TADs (<=",tadSignifThresh , ") - ", curr_hicds, " matching ", match_hicds, " (", exprds, ")")
    # all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
    # outHeightGG <- min(c(7 * nPlotted/2, 49))
    # outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
    # outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
    # ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
    # cat("... written: ", outFile, "\n")
    # stop("--ok\n")
  }
  ############################################# LOLLIPLOT FOR THE TOP P VAL
  out_signif_dt$adjPval_rank <- rank(out_signif_dt$adjPval, ties="min")
  out_signif_dt <- out_signif_dt[order(out_signif_dt$adjPval),]
  plotList <- list()
  toplot_tads <- out_signif_dt$refID[1:nPlotted]
  out_signif_dt$rankDiff_rank <- rank(abs(out_signif_dt$rankDiff), ties="min")
  
  
   
  
  
  if(toplot) {
    
    
    pipOutFolder <- file.path("PIPELINE/OUTPUT_FOLDER/")
    # ADDED -> load 1x/dataset
    rnaseqDTFile <- file.path(pipOutFolder, curr_hicds, exprds,script0_name, "rna_rnaseqDT.Rdata")
    initListFile <- file.path(pipOutFolder, curr_hicds, exprds,script0_name, "rna_geneList.Rdata")
    geneListFile <- file.path(pipOutFolder, curr_hicds, exprds,script0_name, "pipeline_geneList.Rdata")
    DE_topTableFile <- file.path(pipOutFolder, curr_hicds, exprds, script1_name, "DE_topTable.Rdata")
    entrezDT_file <- file.path("/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
    gene2tadDT_file <- file.path(curr_hicds, "genes2tad", "all_genes_positions.txt")
    cat(paste0("... load gene2tad_dt\n"))
    curr_gene2tadDT <-  read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
    cat(paste0("... load entrez2symb_dt\n"))
    curr_entrez2symbDT<- entrez2symb_dt
    cat(paste0("... load DE_topTable\n"))
    curr_DE_topTable <- eval(parse(text = load(DE_topTableFile))) 
    cat(paste0("... load rnaseqDT\n"))
    curr_rnaseqDT <- eval(parse(text = load(rnaseqDTFile))) 
    cat(paste0("... load initList\n"))
    curr_initList <- eval(parse(text = load(initListFile)))	
    cat(paste0("... load geneList\n"))
    curr_geneList <- eval(parse(text = load(geneListFile)))   
    
    
    for(i_tad in 1:nPlotted) {
      tad <- toplot_tads[i_tad]
      mytit <- paste0( curr_hicds, " - ", exprds, " - ", tad, "\n", 
                       "adj. pval=", round(out_signif_dt$adjPval[i_tad], 2)," (# ", i_tad, "; abs. rank diff. rank: ", out_signif_dt$rankDiff_rank[i_tad], "/", max(out_signif_dt$rankDiff_rank), ")")
      
      

      
      
      plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                            hicds = curr_hicds,
                                            all_TADs = tad,
                                            
                                            gene2tadDT=curr_gene2tadDT,entrez2symbDT=curr_entrez2symbDT,
                                            DE_topTable=curr_DE_topTable,
                                            rnaseqDT=curr_rnaseqDT, initList=curr_initList, geneList=curr_geneList,
                                            
                                            
                                            orderByLolli = "startPos", mytitle=mytit)
      
      nPlotted_s <- nPlotted
      nPlotted <- 1
      outFile <- file.path(outFolder, paste0(i_tad, "_", curr_hicds, "_withMatching_", match_hicds, "_", exprds, "_adjPval_top", nPlotted, "_lolli.", plotType))
      # mytit <- paste0("Top ", nPlotted, " rankDiff signif. TADs (<=",tadSignifThresh , ") - ", curr_hicds, " matching ", match_hicds, " (", exprds, ")")
      mytit <- paste0("Top # ", i_tad, " rankDiff signif. TADs (<=",tadSignifThresh , ") - ", curr_hicds, " matching ", match_hicds, " (", exprds, ")")
      # all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      all_plots <- do.call(grid.arrange, c(plotList[[i_tad]],  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      outHeightGG <- min(c(7 * nPlotted/2, 49))
      outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
      outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
      # ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
      ggsave(filename = outFile, plotList[[i_tad]] , width=10, height = 7)
      cat("... written: ", outFile, "\n")
      nPlotted <- nPlotted_s
      
    } # end-for iterating over TADs to plot
    
    # outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_", exprds, "_adjPval_top", nPlotted, "_lolli.", plotType))
    # mytit <- paste0("Top ", nPlotted, " rankDiff signif. TADs (<=",tadSignifThresh , ") - ", curr_hicds, " matching ", match_hicds, " (", exprds, ")")
    # all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
    # outHeightGG <- min(c(7 * nPlotted/2, 49))
    # outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
    # outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
    # ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
    # cat("... written: ", outFile, "\n")
  }

  # stop("--ok\n")
  
}




######################################################################################
######################################################################################
######################################################################################
# cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



