options(scipen=100)

SSHFS=F

buildData <- TRUE
# FOCUS OSN TAD !

# Rscript report_figure6_report.R GSE118514_RWPE1_40kb GSE118514_22Rv1_40kb TCGAprad_norm_prad 

script_name <- "report_figure6_report.R"

plotLolli <- FALSE

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(ggpubr)
require(ggrepel)

hicds_norm <- "LI_40kb"
hicds_tumor <- "GSE105381_HepG2_40kb"
exprds <- "TCGAlihc_norm_lihc"

hicds_norm <- "LG1_40kb"
hicds_tumor <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAlusc_norm_lusc"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

source("plot_lolliTAD_funct.R")
source("my_heatmap.2.R")
# source("annotated_TADs.R"); stopifnot(exists("annotated_tads"))
annotated_tads <- "GSE118514_RWPE1_40kb/GSE118514_22Rv1_40kb/TCGAprad_norm_prad/chr19_TAD172"


script0_name <- "0_prepGeneData"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.2

myHeightGG <- 7
myWidthGG <- 7

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
hicds_norm <- args[1]
hicds_tumor <- args[2]
exprds <- args[3]


inFile <- file.path("TAD_MATCHING_ACROSS_HICDS", "all_matching_dt.Rdata")
matching_dt <- get(load(inFile))

outFolder <- file.path("REPORT_FIGURE6_REPORT")
dir.create(outFolder, recursive = TRUE)

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))


final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

tadSignifThresh <- 0.01

nPlotted <- 10

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
norm_TAD_adjPvals <- get(load(norm_tadPval_file))
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
tumor_TAD_adjPvals <- get(load(tumor_tadPval_file))
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
  
  curr_dt$textCol <- ifelse( file.path(curr_dt$ref_hicds, curr_dt$matching_hicds,curr_dt$ref_exprds, curr_dt$refID) %in% annotated_tads, "red", "black" )
  
  myx <- -log10(curr_dt$adjPval)
  myy <- curr_dt$rankDiff
  
  outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_rankDiff_vs_pval_", exprds, "_densplot.", plotType ))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = myx,
    y = myy,
    main=paste0(exprds),
    sub=paste0(curr_ref, " as refDS"),
    xlab=paste0("-log10 TAD adj. pval"),
    ylab=paste0("best matching TAD rank diff. (", curr_ref, "-", curr_match, ")"),
    cex.axis=axisCex,
    cex.lab=axisCex
  )
  addCorr(x = myx, y=myy, bty="n")
  abline(h=0, lty=2, col="grey")
  mtext(side=3, paste0(curr_hicds, " matching ", match_hicds))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  myx <- curr_dt[, paste0(curr_ref, "MeanFC")]
  
  outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_", curr_ref, "MeanFC_vs_pval_", exprds, "_densplot.", plotType ))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = myx,
    y = myy,
    main=paste0(exprds),
    sub=paste0(curr_ref, " as refDS"),
    xlab=paste0("mean TAD logFC"),
    ylab=paste0("best matching TAD rank diff. (", curr_ref, " - ", curr_match, ")"),
    cex.axis=axisCex,
    cex.lab=axisCex
  )
  addCorr(x = myx, y=myy, bty="n")
  abline(h=0, lty=2, col="grey")
  mtext(side=3, paste0(curr_hicds, " matching ", match_hicds))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  ###################################################################################### plot only signif TADs  

  
  signif_dt <- curr_dt[curr_dt$adjPval <= tadSignifThresh,]
  
  signif_dt <- signif_dt[order(abs(signif_dt$rankDiff), decreasing = TRUE),]
  
  signif_dt[,paste0("adjPval_log10")] <- -log10(signif_dt$adjPval)
  
  myx <- -log10(signif_dt$adjPval)
  myy <- signif_dt$rankDiff
  
  mylab <- paste0("best matching TAD rank diff. (", curr_ref, "-", curr_match, ")")
  
  
  outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_rankDiff_vs_pval_", exprds, "_densplot_signifTADsOnly.", plotType ))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = myx,
    y = myy,
    main=paste0(exprds),
    sub=paste0(curr_ref, " as refDS (tadPval<=", tadSignifThresh, ")" ),
    xlab=paste0("-log10 TAD adj. pval"),
    ylab=mylab,
    cex.axis=axisCex,
    cex.lab=axisCex
  )
  # add text label for the top rankDiff
  text(x=myx[1:nPlotted], y=myy[1:nPlotted], labels = signif_dt$refID[1:nPlotted], cex=0.6, col = signif_dt$textCol[1:nPlotted], pos=3)
  
  
  addCorr(x = myx, y=myy, bty="n")
  abline(h=0, lty=2, col="grey")
  mtext(side=3, paste0(curr_hicds, " matching ", match_hicds))
  # add text label for the top pvals
  signif_dt_sortPval <- signif_dt[order(signif_dt$adjPval),]
  myx_sortPval <- -log10(signif_dt_sortPval$adjPval)
  myy_sortPval <- signif_dt_sortPval$rankDiff
  signif_dt_sortPval$textCol <- ifelse( file.path(signif_dt_sortPval$ref_hicds, signif_dt_sortPval$ref_exprds, signif_dt_sortPval$refID) %in% annotated_tads, "red", "darkgrey" )
  text(x=myx_sortPval[1:nPlotted], y=myy_sortPval[1:nPlotted], labels = signif_dt_sortPval$refID[1:nPlotted], cex=0.6, col = signif_dt_sortPval$textCol[1:nPlotted], pos=3)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  myx <- signif_dt[, paste0(curr_ref, "MeanFC")]
  
  outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_", curr_ref, "MeanFC_vs_pval_", exprds, "_densplot_siginfTADsOnly.", plotType ))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = myx,
    y = myy,
    main=paste0(exprds),
    sub=paste0(curr_ref, " as refDS (tadPval<=", tadSignifThresh, ")" ),
    xlab=paste0("mean TAD logFC"),
    ylab=paste0("best matching TAD rank diff. (", curr_ref, " - ", curr_match, ")"),
    cex.axis=axisCex,
    cex.lab=axisCex
  )
  text(x=myx[1:nPlotted], y=myy[1:nPlotted], labels = signif_dt$refID[1:nPlotted], cex=0.6, col = signif_dt$textCol[1:nPlotted])
  addCorr(x = myx, y=myy, bty="n")
  abline(h=0, lty=2, col="grey")
  mtext(side=3, paste0(curr_hicds, " matching ", match_hicds))
  
  # add text label for the top pvals
  signif_dt_sortPval <- signif_dt[order(signif_dt$adjPval),]
  myx_sortPval <- signif_dt_sortPval[, paste0(curr_ref, "MeanFC")]
  myy_sortPval <- signif_dt_sortPval$rankDiff
  signif_dt_sortPval$textCol <- ifelse( file.path(signif_dt_sortPval$ref_hicds, signif_dt_sortPval$ref_exprds, signif_dt_sortPval$refID) %in% annotated_tads, "red", "darkgrey" )
  text(x=myx_sortPval[1:nPlotted], y=myy_sortPval[1:nPlotted], labels = signif_dt_sortPval$refID[1:nPlotted], cex=0.6, col = signif_dt_sortPval$textCol[1:nPlotted], pos=3)
  
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  save(signif_dt, file="signif_dt.Rdata", version=2)
  
  
  mytit <- paste0("TCGA PRAD")
  mysubtit <- paste0("normal vs. tumor")
  mylab <- paste0("\u0394 TAD rank score")
  
  p <- ggscatter(signif_dt,
                 title=paste0(mytit),
                 # subtitle=paste0(paste0(curr_hicds, " matching ", match_hicds), "\n", curr_ref, " as refDS (tadPval<=", tadSignifThresh, ")" ),
                 subtitle=paste0(mysubtit),
                 col=signif_dt$textCol,
                 xlab = mylab,
                 ylab = paste0("TAD adj. p-val. [-log10]"),
                 y="adjPval_log10", 
                 x="rankDiff")+ 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8))                    
  
  
  p <- p + geom_vline(xintercept=0, linetype=2)
  
  p <- p + theme(plot.title = element_text(hjust=0.5, face="bold", size=16),
                 plot.subtitle = element_text(hjust=0.5, face="italic", size=14),
                 axis.title = element_text(size=14),
                 axis.text = element_text(size=12)
                 )
  add_dt_1 <- signif_dt[order(signif_dt$adjPval), ][1:nPlotted,]
  add_dt_2 <- signif_dt[order(abs(signif_dt$rankDiff), decreasing=TRUE),][1:nPlotted,]
  add_dt <- rbind(add_dt_1, add_dt_2)
  add_dt <- unique(add_dt)
  
  p <- p+
    geom_text_repel(data = add_dt,
                    aes(y=adjPval_log10,x=rankDiff,label=refID), col=add_dt$textCol)
  # p
  # 
  # 
  # signif_dt <- signif_dt[order(signif_dt$adjPval),]
  # p <- p+
  #   geom_text_repel(data = signif_dt[1:nPlotted,],
  #                   aes(y=adjPval_log10,x=rankDiff,label=refID), col=signif_dt$textCol[1:nPlotted])
  # 
  # signif_dt <- signif_dt[order(abs(signif_dt$rankDiff), decreasing=TRUE),]
  # p <- p+
  #   geom_text_repel(data = signif_dt[1:nPlotted,],
  #                   aes(y=adjPval_log10,x=rankDiff,label=refID), col=signif_dt$textCol[1:nPlotted])
  # 
  # p
  
    
  outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_rankDiff_vs_pval_", exprds, "_ggrepel_signifTADsOnly.", plotType ))
  ggsave(plot = p, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  # stop("-ok")
  
  
  
  
  
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
  
  if(plotLolli) {
    
    plotList <- list()
    
    toplot_tads <- out_signif_dt$refID[1:nPlotted]
    
    ############################################# LOLLIPLOT FOR THE TOP RANK DIFF
    # foo <- foreach(i_tad = 1:nPlotted) %dopar% {
    for(i_tad in 1:nPlotted) {
      tad <- toplot_tads[i_tad]
      mytit <- paste0( curr_hicds, " - ", exprds, " - ", tad, "\n", "rankDiff=", round(out_signif_dt$rankDiff[i_tad], 2)," (adj. pval rank: ", out_signif_dt$adjPval_rank[i_tad], "/", max(out_signif_dt$adjPval_rank), ")")
      plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                            hicds = curr_hicds,
                                            all_TADs = tad,
                                            orderByLolli = "startPos", mytitle=mytit)
    } # end-for iterating over TADs to plot
    
    outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_", exprds, "_rankDiff_top", nPlotted, "_lolli.", plotType))
    
    mytit <- paste0("Top ", nPlotted, " rankDiff signif. TADs (<=",tadSignifThresh , ") - ", curr_hicds, " matching ", match_hicds, " (", exprds, ")")
    all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
    outHeightGG <- min(c(7 * nPlotted/2, 49))
    outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
    outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
    ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
    cat("... written: ", outFile, "\n")
    # stop("--ok\n")
    
    ############################################# LOLLIPLOT FOR THE TOP P VAL
    out_signif_dt$adjPval_rank <- rank(out_signif_dt$adjPval, ties="min")
    out_signif_dt <- out_signif_dt[order(out_signif_dt$adjPval),]
    plotList <- list()
    toplot_tads <- out_signif_dt$refID[1:nPlotted]
    out_signif_dt$rankDiff_rank <- rank(abs(out_signif_dt$rankDiff), ties="min")
    
    for(i_tad in 1:nPlotted) {
      tad <- toplot_tads[i_tad]
      mytit <- paste0( curr_hicds, " - ", exprds, " - ", tad, "\n", 
                       "adj. pval=", round(out_signif_dt$adjPval[i_tad], 2)," (abs. rank diff. rank: ", out_signif_dt$rankDiff_rank[i_tad], "/", max(out_signif_dt$rankDiff_rank), ")")
      plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                            hicds = curr_hicds,
                                            all_TADs = tad,
                                            orderByLolli = "startPos", mytitle=mytit)
    } # end-for iterating over TADs to plot
    
    outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_", exprds, "_adjPval_top", nPlotted, "_lolli.", plotType))
    
    mytit <- paste0("Top ", nPlotted, " rankDiff signif. TADs (<=",tadSignifThresh , ") - ", curr_hicds, " matching ", match_hicds, " (", exprds, ")")
    all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
    outHeightGG <- min(c(7 * nPlotted/2, 49))
    outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
    outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
    ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
    cat("... written: ", outFile, "\n")
    # stop("--ok\n")
    
    
  }
  
  
}




######################################################################################
######################################################################################
######################################################################################
# cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt
... load //mnt/ed4/marie/other_datasets/TCGAprad_norm_prad/norm_ID.Rdata
... load //mnt/ed4/marie/other_datasets/TCGAprad_norm_prad/prad_ID.Rdata
... load //mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE118514_22Rv1_40kb/TCGAprad_norm_prad/0_prepGeneData/rna_rnaseqDT.Rdata
... load //mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE118514_22Rv1_40kb/TCGAprad_norm_prad/0_prepGeneData/rna_geneList.Rdata
... load //mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/GSE118514_22Rv1_40kb/TCGAprad_norm_prad/0_prepGeneData/pipeline_geneList.Rdata
... start plotting
... written:  REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_TCGAprad_norm_prad_adjPval_top10_lolli.png 
There were 42 warnings (use warnings() to see them)
*** DONE
2019-10-21 17:18:47
2019-10-21 17:20:56
Error: object 'ata' not found
Execution halted
kmarie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA\[0;31mmarie@electron[00m:[01;32m/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA[00m$ Rscript report_figure6_report.R GSE118514_RWPE1_40kb GSE118514_22Rv1_40kb TCGAprad_norm_prad 
> START  report_figure6_report.R 
Loading required package: foreach
Loading required package: doMC
Loading required package: iterators
Loading required package: parallel
Loading required package: ggplot2
Loading required package: ggpubr
Loading required package: magrittr
Loading required package: ggrepel
Loading required package: gtools
Warning message:
In dir.create(outFolder, recursive = TRUE) :
  'REPORT_FIGURE6_REPORT' already exists
... written: REPORT_FIGURE6_REPORT/tumor_matching_pval_tadRank_dt.Rdata
... written: REPORT_FIGURE6_REPORT/norm_matching_pval_tadRank_dt.Rdata
Error: unexpected symbol in:
"                 subtitle=paste0(mysubtit)
                 col"
Execution halted
kmarie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA\[0;31mmarie@electron[00m:[01;32m/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA[00m$ Rscript report_figure6_report.R GSE118514_RWPE1_40kb GSE118514_22Rv1_40kb TCGAprad_norm_prad 
> START  report_figure6_report.R 
Loading required package: foreach
Loading required package: doMC
Loading required package: iterators
Loading required package: parallel
Loading required package: ggplot2
Loading required package: ggpubr
Loading required package: magrittr
Loading required package: ggrepel
Loading required package: gtools
Warning message:
In dir.create(outFolder, recursive = TRUE) :
  'REPORT_FIGURE6_REPORT' already exists
... written: REPORT_FIGURE6_REPORT/tumor_matching_pval_tadRank_dt.Rdata
... written: REPORT_FIGURE6_REPORT/norm_matching_pval_tadRank_dt.Rdata
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_normMeanFC_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_normMeanFC_vs_pval_TCGAprad_norm_prad_densplot_siginfTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_ggrepel_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_TCGAprad_norm_prad_signifDT.txt
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_tumorMeanFC_vs_pval_TCGAprad_norm_prad_densplot.RWPE1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_tumorMeanFC_vs_pval_TCGAprad_norm_prad_densplot_siginfTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_ggrepel_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_TCGAprad_norm_prad_signifDT.txt
Warning messages:
1: In if (color %in% names(data) & is.null(add.params$color)) add.params$color <- color :
  the condition has length > 1 and only the first element will be used
2: In if (color %in% names(data) & is.null(add.params$color)) add.params$color <- color :
  the condition has length > 1 and only the first element will be used
*** DONE
2019-10-21 17:43:54
2019-10-21 17:44:30
kmarie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA\[0;31mmarie@electron[00m:[01;32m/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA[00m$ Rscript report_figure6_report.R GSE118514_RWPE1_40kb GSE118514_22Rv1_40kb TCGAprad_norm_prad 
> START  report_figure6_report.R 
Loading required package: foreach
Loading required package: doMC
Loading required package: iterators
Loading required package: parallel
Loading required package: ggplot2
Loading required package: ggpubr
Loading required package: magrittr
Loading required package: ggrepel
Loading required package: gtools
Warning message:
In dir.create(outFolder, recursive = TRUE) :
  'REPORT_FIGURE6_REPORT' already exists
... written: REPORT_FIGURE6_REPORT/tumor_matching_pval_tadRank_dt.Rdata
... written: REPORT_FIGURE6_REPORT/norm_matching_pval_tadRank_dt.Rdata
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_normMeanFC_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_normMeanFC_vs_pval_TCGAprad_norm_prad_densplot_siginfTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_ggrepel_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_TCGAprad_norm_prad_signifDT.txt
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_tumorMeanFC_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_tumorMeanFC_vs_pval_TCGAprad_norm_prad_densplot_siginfTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_ggrepel_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_TCGAprad_norm_prad_signifDT.txt
Warning messages:
1: In if (color %in% names(data) & is.null(add.params$color)) add.params$color <- color :
  the condition has length > 1 and only the first element will be used
2: In if (color %in% names(data) & is.null(add.params$color)) add.params$color <- color :
  the condition has length > 1 and only the first element will be used
*** DONE
2019-10-21 17:46:33
2019-10-21 17:47:08
kmarie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA\[0;31mmarie@electron[00m:[01;32m/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA[00m$ Rscript report_figure6_report.R GSE118514_RWPE1_40kb GSE118514_22Rv1_40kb TCGAprad_norm_prad 
> START  report_figure6_report.R 
Loading required package: foreach
Loading required package: doMC
Loading required package: iterators
Loading required package: parallel
Loading required package: ggplot2
Loading required package: ggpubr
Loading required package: magrittr
Loading required package: ggrepel
Loading required package: gtools
Warning message:
In dir.create(outFolder, recursive = TRUE) :
  'REPORT_FIGURE6_REPORT' already exists
... written: REPORT_FIGURE6_REPORT/tumor_matching_pval_tadRank_dt.Rdata
... written: REPORT_FIGURE6_REPORT/norm_matching_pval_tadRank_dt.Rdata
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_normMeanFC_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_normMeanFC_vs_pval_TCGAprad_norm_prad_densplot_siginfTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_ggrepel_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_RWPE1_40kb_withMatching_GSE118514_22Rv1_40kb_TCGAprad_norm_prad_signifDT.txt
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_tumorMeanFC_vs_pval_TCGAprad_norm_prad_densplot.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_densplot_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_tumorMeanFC_vs_pval_TCGAprad_norm_prad_densplot_siginfTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_rankDiff_vs_pval_TCGAprad_norm_prad_ggrepel_signifTADsOnly.png
... written: REPORT_FIGURE6_REPORT/GSE118514_22Rv1_40kb_withMatching_GSE118514_RWPE1_40kb_TCGAprad_norm_prad_signifDT.txt
Warning messages:
1: In if (color %in% names(data) & is.null(add.params$color)) add.params$color <- color :
  the condition has length > 1 and only the first element will be used
2: In if (color %in% names(data) & is.null(add.params$color)) add.params$color <- color :
  the condition has length > 1 and only the first element will be used
*** DONE
2019-10-21 17:55:28
2019-10-21 17:56:05
kmarie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA\[0;31mmarie@electron[00m:[01;32m/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA[00m$ git add *.R
kmarie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA\[0;31mmarie@electron[00m:[01;32m/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA[00m$ (reverse-i-search)`': [K[94@c': Rscript report_figure6_report.R GSE118514_RWPE1_40kb GSE118514_22Rv1_40kb TCGAprad_norm_prad[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[Co': java -Xmx2g -jar  /mnt/ed2/shared/TADcompare/Software/juicer/juicer_tools.jar pre -n -d -r 25000 -c chr6 -q 0 GM12878_chr6_25kb_matrix.pre GM12878_chr6_25kb_matrix.hic chr6.size[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[1@m[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[Cm': git commit -m "19.10.19"[K[0;31mmarie@electron[00m:[01;32m/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA[00m$ git commit -m "19.10.19"[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[1P.10.19"2.10.19"[1P.10.19"[1P.10.19"2.10.19"1.10.19"
[master ad3ef84] 21.10.19
 3 files changed, 1123 insertions(+), 3 deletions(-)
 create mode 100644 go_specificity_geneLevel_tadLevel_intersectDiff_report.R
 create mode 100644 report_figure6_report.R
kmarie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA\[0;31mmarie@electron[00m:[01;32m/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA[00m$ git push
Username for 'https://github.com': marzuf
Password for 'https://marzuf@github.com': 
Counting objects: 7, done.
Delta compression using up to 44 threads.
Compressing objects:  20% (1/5)   Compressing objects:  40% (2/5)   Compressing objects:  60% (3/5)   Compressing objects:  80% (4/5)   Compressing objects: 100% (5/5)   Compressing objects: 100% (5/5), done.
Writing objects:  20% (1/5)   Writing objects:  40% (2/5)   Writing objects:  60% (3/5)   Writing objects:  80% (4/5)   Writing objects: 100% (5/5)   Writing objects: 100% (5/5), 9.47 KiB | 0 bytes/s, done.
Total 5 (delta 2), reused 0 (delta 0)
remote: Resolving deltas:   0% (0/2)[Kremote: Resolving deltas:  50% (1/2)[Kremote: Resolving deltas: 100% (2/2)[Kremote: Resolving deltas: 100% (2/2), completed with 2 local objects.[K
To https://github.com/marzuf/v2_Yuanlong_Cancer_HiC_data_TAD_DA.git
   6afb45d..ad3ef84  master -> master
kmarie@electron:/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA\[0;31mmarie@electron[00m:[01;32m/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA[00m$ exit
exit

Script done on Thu 24 Oct 2019 11:55:19 AM CEST
