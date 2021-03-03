options(scipen=100)

SSHFS=F

set.seed(010321) # for sampling

# Rscript revision_rankDiff_activDiff_v2.R

buildPlotList <- F#TRUE
buildTable <- F#T
runPermut <- T

nRandom_stats <- 100000

# copied from OLDER_v2_Yuanlong..._DA/report_figure6_all.R
script_name <- "revision_rankDiff_activDiff_v2.R"

plotLolli <- FALSE

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(ggpubr)
require(ggrepel)
require(ggsci)
require(gridExtra)

require(e1071)


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

source("revision_settings.R")

script0_name <- "0_prepGeneData"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.2

myHeightGG <- 7
myWidthGG <- 7

# args <- commandArgs(trailingOnly = TRUE)
# stopifnot(length(args) == 3)
# hicds_norm <- args[1]
# hicds_tumor <- args[2]
# exprds <- args[3]

mypair=all_pairs[1]

outFolder <- file.path("REVISION_RANKDIFF_ACTIVDIFF_V2")
dir.create(outFolder, recursive = TRUE)

tadSignifThresh <- 0.01

nPlotted <- 10

quantileThresh <- 0.9
qtcol <- "red"

matchingCol <- "matchingID_maxOverlapBp"

outHeightGG <- 3*3
outWidthGG <- 3*4

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)

if(buildTable){
  # inFile <- file.path("TAD_MATCHING_ACROSS_HICDS", "all_matching_dt.Rdata")
  inFile <- file.path("TAD_MATCHING_ACROSS_HICDS_NOPERMUT", "all_matching_dt.Rdata")
  matching_dt <- get(load(inFile))
  
  all_plots <- list()

 matching_data <- foreach(mypair = all_pairs) %dopar% {
    hicds_norm <- dirname(dirname(mypair))
    hicds_tumor <- basename(dirname(mypair))
    exprds <- basename(mypair)
    
    pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
    stopifnot(dir.exists(pipFolder))
    
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
    norm_TAD_pvals <- get(load(norm_tadPval_file))
    norm_TAD_adjPvals <- p.adjust(norm_TAD_pvals, method="BH")
    norm_TAD_adjPvals_dt <- data.frame(refID = names(norm_TAD_adjPvals), adjPval=as.numeric(norm_TAD_adjPvals), stringsAsFactors = FALSE)
    
    
    stopifnot(sum(norm_TAD_adjPvals_dt$adjPval <= tadSignifThresh) == 
      sum(final_table_DT$adjPvalComb[final_table_DT$hicds==hicds_norm &
                                       final_table_DT$exprds == exprds ] <= tadSignifThresh))
    
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
    tumor_TAD_pvals <- get(load(tumor_tadPval_file))
    tumor_TAD_adjPvals <- p.adjust(tumor_TAD_pvals, method="BH")
    tumor_TAD_adjPvals_dt <- data.frame(refID = names(tumor_TAD_adjPvals), 
                                        adjPval=as.numeric(tumor_TAD_adjPvals), stringsAsFactors = FALSE)
    
    stopifnot(sum(tumor_TAD_adjPvals_dt$adjPval <= tadSignifThresh) == 
                sum(final_table_DT$adjPvalComb[final_table_DT$hicds==hicds_tumor &
                                                 final_table_DT$exprds == exprds ] <= tadSignifThresh))
    
    
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
    curr_ref = "norm"

    if(buildPlotList) {

      for(curr_ref in all_refs) {

        curr_hicds <- eval(parse(text = paste0("hicds_", curr_ref)))
        curr_match <- all_refs[all_refs != curr_ref]
        stopifnot(length(curr_match) == 1)
        match_hicds <- eval(parse(text = paste0("hicds_", curr_match)))

        curr_dt <- eval(parse(text=paste0(curr_ref, "_matching_pval_tadRank_dt")))

        curr_dt$textCol <- ifelse( file.path(curr_dt$ref_hicds, curr_dt$matching_hicds, curr_dt$ref_exprds, curr_dt$refID) %in% annotated_tads, "red", "black" )

        myx <- -log10(curr_dt$adjPval)
        myy <- curr_dt$rankDiff

        myx <- curr_dt[, paste0(curr_ref, "MeanFC")]

        ###################################################################################### plot only signif TADs

        signif_dt <- curr_dt[curr_dt$adjPval <= tadSignifThresh,]

        signif_dt <- signif_dt[order(abs(signif_dt$rankDiff), decreasing = TRUE),]

        signif_dt[,paste0("adjPval_log10")] <- -log10(signif_dt$adjPval)

        myx <- -log10(signif_dt$adjPval)
        myy <- signif_dt$rankDiff

        mylab <- paste0("best matching TAD rank diff. (", curr_ref, "-", curr_match, ")")
        save(signif_dt, file="signif_dt.Rdata", version=2)


        p <- ggscatter(signif_dt,
                       title=paste0(exprds),
                       ytics = pretty(signif_dt[,"adjPval_log10"], n=20),
                       xtics = pretty(signif_dt[,"rankDiff"], n=20),
                       subtitle=paste0(paste0(curr_hicds, " matching ", match_hicds), "\n", curr_ref, " as refDS (tadPval<=", tadSignifThresh, ")" ),
                       col=signif_dt$textCol,
                       xlab = mylab,
                       ylab = paste0("-log10 TAD adj. pval"),
                       y="adjPval_log10",
                       x="rankDiff")

        p <- p + geom_vline(xintercept=0, linetype=2)

        p <- p + theme(plot.title = element_text(hjust=0.5, face="bold", size=16),
                       plot.subtitle = element_text(hjust=0.5, face="italic", size=14)
        )
        add_dt_1 <- signif_dt[order(signif_dt$adjPval), ][1:nPlotted,]
        add_dt_2 <- signif_dt[order(abs(signif_dt$rankDiff), decreasing=TRUE),][1:nPlotted,]
        add_dt <- rbind(add_dt_1, add_dt_2)
        add_dt <- unique(add_dt)

        p <- p+
          geom_text_repel(data = add_dt,
                          aes(y=adjPval_log10,x=rankDiff,label=refID), col=add_dt$textCol)

        outFile <- file.path(outFolder, paste0(curr_hicds, "_withMatching_", match_hicds, "_rankDiff_vs_pval_", exprds, "_ggrepel_signifTADsOnly.", plotType ))
        ggsave(plot = p, file=outFile, height=myHeightGG, width=myWidthGG)
        cat(paste0("... written: ", outFile, "\n"))

        to_save_name <- paste0(mypair, "_", curr_hicds, "_asRef")

        all_plots[[paste0(to_save_name)]] <- p
      }
    outFile <- file.path(outFolder, "all_plots.Rdata")
    save(all_plots, file =outFile, version=2)
    cat("... written: ", outFile, "\n")
  } else { # if/else buildplotlist
    outFile <- file.path(outFolder, "all_plots.Rdata")
    all_plots <- get(load(outFile))
  }
  all_plots_sub <- all_plots[seq(1, length(all_plots), by=2)]

  out_p <- do.call(grid.arrange, c(all_plots_sub,  list(ncol=4)))

  outFile <- file.path(outFolder, paste0("all_normAsRef_pval_vs_rank.", plotType))
  ggsave(filename = outFile, out_p, width=outWidthGG*2.5, height = outHeightGG*2.5)
  cat("... written: ", outFile, "\n")
  
  list(
    norm_matching_pval_tadRank_dt=norm_matching_pval_tadRank_dt,
    tumor_matching_pval_tadRank_dt=tumor_matching_pval_tadRank_dt#, too heavy to store the plots !!!
    # all_plots=all_plots
  )
  
 } # end foreach  
 names(matching_data) <- all_pairs
 outFile <- file.path(outFolder, "matching_data.Rdata")
 save(matching_data, file=outFile, version=2)
 cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "matching_data.Rdata")
  # matching_data <- get(load("REVISION_RANKDIFF_ACTIVDIFF/matching_data.Rdata"))
  matching_data <- get(load(file=outFile))
}

ds1_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["norm_matching_pval_tadRank_dt"]]))
unique(ds1_matching_dt$ref_hicds)
# [1] "LI_40kb"              "LG1_40kb"             "LG2_40kb"             "GSE118514_RWPE1_40kb"
ds2_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["tumor_matching_pval_tadRank_dt"]]))
unique(ds2_matching_dt$ref_hicds)
# [1] "GSE105381_HepG2_40kb"      "ENCSR444WCZ_A549_40kb"     "ENCSR489OCU_NCI-H460_40kb"
# [4] "ENCSR346DCU_LNCaP_40kb"    "GSE118514_22Rv1_40kb"     


matching_withRank_dt <- rbind(ds1_matching_dt, ds2_matching_dt)
rownames(matching_withRank_dt) <- NULL

stopifnot(matching_withRank_dt$matching_exprds == matching_withRank_dt$ref_exprds )

matching_withRank_dt$ref_region_ID <- file.path(matching_withRank_dt$ref_hicds,
                                                matching_withRank_dt$ref_exprds,
                                                matching_withRank_dt$refID
                                                )
matching_withRank_dt$matching_region_ID <- file.path(matching_withRank_dt$matching_hicds,
                                                matching_withRank_dt$matching_exprds,
                                                matching_withRank_dt$matchingID_maxOverlapBp)

stopifnot(matching_withRank_dt$ref_region_ID %in% names(regionID_pvals))
matching_withRank_dt$ref_region_pval <- regionID_pvals[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_pval))
stopifnot(round(matching_withRank_dt$ref_region_pval, 6) == round(matching_withRank_dt$adjPval, 6))

stopifnot(matching_withRank_dt$matching_region_ID %in% names(regionID_pvals))
matching_withRank_dt$matching_region_pval <- regionID_pvals[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_pval))

matching_withRank_dt$ref_tadSignif <- ifelse(matching_withRank_dt$adjPval <= tadSignifThresh, "signif.", "not signif.")
my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(matching_withRank_dt$ref_tadSignif))

legTitle <- ""

matching_withRank_dt$abs_rankDiff <- abs(matching_withRank_dt$rankDiff)
qt_rankdiff <- quantile(matching_withRank_dt$abs_rankDiff, probs=quantileThresh)

##################
plot_dt <- matching_withRank_dt
outFile <- file.path(outFolder, "plot_dt.Rdata")
save(plot_dt, file=outFile)
cat(paste0("... written: ", outFile, "\n"))
# plot_dt=get(load("REVISION_RANKDIFF_ACTIVDIFF/plot_dt.Rdata"))

# ######################################################################################## PERMUTATION TEST - V1

# 

nsignif <- sum(plot_dt$ref_tadSignif == "signif.")
nnonsignif <- sum(plot_dt$ref_tadSignif == "not signif.")

col_var <- "abs_rankDiff"

cat(paste0("... start permut\n"))
# if(runPermut) {
  permut_data <- foreach(i= 1:nRandom_stats) %dopar% {
    
    i_rd_signif <- sample(seq_len(nrow(plot_dt)), size=nsignif)
    i_rd_nonsignif <- setdiff(seq_len(nrow(plot_dt)), i_rd_signif)
    stopifnot(c(1:nrow(plot_dt)) ==  sort(c(i_rd_signif, i_rd_nonsignif)))
    
    rd_signif <- plot_dt[i_rd_signif, paste0(col_var)]
    stopifnot(length(rd_signif) == nsignif)
    
    rd_nonsignif <- plot_dt[i_rd_nonsignif, paste0(col_var)]
    stopifnot(length(rd_nonsignif) == nnonsignif)
    
    list(
        # all_signif = rd_signif,  # too heavy
        #  all_nonsignif = rd_nonsignif,
      all_signif = sample(rd_signif, size=1),  # too heavy
       all_nonsignif = sample(rd_nonsignif, size=1),
        mean_signif = mean(rd_signif),
        mean_nonsignif = mean(rd_nonsignif),
         geQt_signif=sum(rd_signif >= qt_rankdiff), 
         geQt_nonsignif= sum(rd_nonsignif >= qt_rankdiff))
  }
  # save(permut_data, file=file.path(outFolder, "permut_data.Rdata"), version=2)
# } else {
  # permut_data <- get(load(file.path(outFolder, "permut_data.Rdata")))
# }
  cat(paste0("... permut done\n"))

tmp_ctg_dt <- matching_withRank_dt
tmp_ctg_dt$aboveQt <- tmp_ctg_dt$abs_rankDiff >= qt_rankdiff
ctg_mat <- table(tmp_ctg_dt$aboveQt, tmp_ctg_dt$ref_tadSignif)
names(dimnames(ctg_mat)) <- c("GeQt", "TAD signif." )
print(ctg_mat)

all_permut_geQt_signif <- unlist(lapply(permut_data, function(x)x[["geQt_signif"]]))
all_permut_geQt_nonsignif <- unlist(lapply(permut_data, function(x)x[["geQt_nonsignif"]]))

obs_geQt_signif <- ctg_mat["TRUE", "signif."]
obs_geQt_nonsignif <- ctg_mat["TRUE", "not signif."]

empPval_signif <- (sum(all_permut_geQt_signif >= obs_geQt_signif) + 1)/(nRandom_stats+1)
empPval_nonsignif <- (sum(all_permut_geQt_nonsignif >= obs_geQt_nonsignif) + 1)/(nRandom_stats+1)


all_cmps <- unique(file.path(matching_withRank_dt$matching_hicds, matching_withRank_dt$matching_exprds,
                             matching_withRank_dt$ref_hicds, matching_withRank_dt$ref_exprds))

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_withRank_dt), 
                " (signif.: ", sum(matching_withRank_dt$adjPval <= tadSignifThresh), ")")

mylineTxt <- paste(
  paste0(quantileThresh, "-qt for all TADs"), 
  paste0(round(100*obs_geQt_signif/sum(ctg_mat[, "signif."]) , 2), "% of signif. TADs"),
  paste0(round(100*obs_geQt_nonsignif/sum(ctg_mat[, "not signif."]) , 2), "% of not signif. TADs"),
  paste0("Fisher test p-val = ", formatC(fisher.test(ctg_mat)$p.value, format = "e", digits = 4)),
  paste0("emp. p-val signif. (n=", nRandom_stats,") = ", formatC(empPval_signif, format = "e", digits = 4)),
  paste0("emp. p-val nonsignif. (n=", nRandom_stats,") = ", formatC(empPval_nonsignif, format = "e", digits = 4)),
  sep="\n")

plotTit <- paste0("TAD abs. rank diff. and signif.")

p3 <- ggdensity(matching_withRank_dt,
                x = paste0("abs_rankDiff"),
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = paste0("Abs. rank diff. with matched TAD"),
                # add = "median",                  # Add median line.
                rug = FALSE,                      # Add marginal rug
                color = "ref_tadSignif",
                fill = "ref_tadSignif",
                palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  scale_color_manual(values=my_cols)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  mytheme +
  geom_vline(xintercept=qt_rankdiff, linetype=2, col=qtcol)


  txtY <- ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range[2]
  
  p3 <- p3 +   annotate("text", 
                  label=mylineTxt, 
                  x=qt_rankdiff, 
                  y =txtY, color = qtcol, hjust = 0, vjust=1)
    
  
outFile <- file.path(outFolder, paste0("tad_absRankDiff_signif_notsignif_density.", plotType))
ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


tmp_dt <- matching_withRank_dt
tmp_dt$ref_tadSignif <- "all"

matching_withRank_dt_all <- rbind(matching_withRank_dt, tmp_dt)

matching_withRank_dt_all$ref_tadSignif <- factor(matching_withRank_dt_all$ref_tadSignif, levels =c("signif.", "not signif.", "all"))
stopifnot(!is.na(matching_withRank_dt_all$ref_tadSignif))


my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(matching_withRank_dt_all$ref_tadSignif))

plotTit <- paste0("TAD abs. rank diff. and signif. with all")

p3 <- ggdensity(matching_withRank_dt_all,
                x = paste0("abs_rankDiff"),
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = paste0("Abs. rank diff. with matched TAD"),
                # add = "median",                  # Add median line.
                rug = FALSE,                      # Add marginal rug
                color = "ref_tadSignif",
                fill = "ref_tadSignif",
                palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  scale_color_manual(values=my_cols)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  mytheme+
  geom_vline(xintercept=qt_rankdiff, linetype=2, col=qtcol)


txtY <- ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range[2]

p3 <- p3 +   annotate("text", 
                      label=mylineTxt, 
                      x=qt_rankdiff, 
                      y =txtY, color = qtcol,
                      hjust = 0, vjust=1)


outFile <- file.path(outFolder, paste0("tad_absRankDiff_signif_notsignif_withAll_density.", plotType))
ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

# save(permut_data, file=file.path(outFolder, "permut_data.Rdata"), version=2)


all_permut_geQt_signif <- unlist(lapply(permut_data, function(x)x[["geQt_signif"]]))
all_permut_geQt_nonsignif <- unlist(lapply(permut_data, function(x)x[["geQt_nonsignif"]]))

# rd_dt <- data.frame(
#   permut_signif =  unlist(lapply(permut_data, function(x)x[["all_signif"]])),
#   permut_nonsignif =  unlist(lapply(permut_data, function(x)x[["all_nonsignif"]])))
# require(reshape2)
# rd_dt_m <- melt(rd_dt)
rd_dt_m1 <- data.frame(variable="signif.",
                       value =  unlist(lapply(permut_data, function(x)x[["all_signif"]])))
rd_dt_m2 <- data.frame(variable="not signif.",
                      value =  unlist(lapply(permut_data, function(x)x[["all_nonsignif"]])))
# rd_dt_m1 <- data.frame(variable="signif.", 
#                        value =  unlist(lapply(permut_data, function(x)x[["mean_signif"]])))
# rd_dt_m2 <- data.frame(variable="not signif.", 
#                        value =  unlist(lapply(permut_data, function(x)x[["mean_nonsignif"]])))


rd_dt_m <- rbind(rd_dt_m1, rd_dt_m2)

my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(rd_dt_m$variable))

mylineTxt <- paste(
paste0(quantileThresh, "-qt for all obs. TADs"), 
# paste0(round(sum(rd_dt_m1$value >= qt_rankdiff)/nrow(rd_dt_m1), 2), "% of rd signif. TADs"),
# paste0(round(sum(rd_dt_m2$value >= qt_rankdiff)/nrow(rd_dt_m2), 2), "% of rd not signif. TADs"),
paste0(round(sum(rd_dt_m1$value >= qt_rankdiff)/nrow(rd_dt_m1), 2), "% of rd mean signif. TADs"),
paste0(round(sum(rd_dt_m2$value >= qt_rankdiff)/nrow(rd_dt_m2), 2), "% of rd mean not signif. TADs"),
sep="\n")

# plotTit <- paste0("TAD abs. rank diff. and signif. - random data")
# plotTit <- paste0("TAD abs. rank diff. and signif. - random data mean")
plotTit <- paste0("TAD abs. rank diff. and signif. - random data (1 value per permut)")

p3 <- ggdensity(rd_dt_m,
                x = paste0("value"),
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = paste0("Abs. rank diff. with matched TAD"),
                # add = "median",                  # Add median line.
                rug = FALSE,                      # Add marginal rug
                color = "variable",
                fill = "variable",
                palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  scale_color_manual(values=my_cols)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  mytheme+
  geom_vline(xintercept=qt_rankdiff, linetype=2, col=qtcol)


txtY <- ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range[2]

p3 <- p3 +   annotate("text", 
                      label=mylineTxt, 
                      x=qt_rankdiff, 
                      y =txtY, color = qtcol,
                      hjust = 0, vjust=1)


outFile <- file.path(outFolder, paste0("tad_absRankDiff_signif_notsignif_randomData_density.", plotType))
ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))



########################## rd above thresh - signif.

my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(rd_dt_m$variable))

mylineTxt <- paste("observed #")

for(plotvar in c("signif", "nonsignif")) {
  
  
  varqt <- eval(parse(text=paste0("obs_geQt_", plotvar)))
  
  rd_dt <- data.frame(
    permut_above =  unlist(lapply(permut_data, function(x)x[[paste0("geQt_", plotvar)]]))
  )
  
  plotTit <- paste0("Random ", plotvar, " TADs above threshold")
  
  mySub <- paste0("# values = ", nrow(rd_dt), " (# permut = ",nRandom_stats,")")
  
  p3 <- ggdensity(rd_dt,
                  x = paste0("permut_above"),
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0("# of TADs above qt. thresh (", plotvar, ".)"),
                  # add = "median",                  # Add median line.
                  rug = FALSE                    # Add marginal rug
  ) +
    ggtitle(plotTit, subtitle = mySub)+
    labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    mytheme+
    geom_vline(xintercept=varqt, linetype=2, col=qtcol)
  
  
  txtY <- ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range[2]
  
  p3 <- p3 +   annotate("text", 
                        label=mylineTxt, 
                        x=varqt, 
                        y =txtY, color = qtcol,
                        hjust = 0, vjust=1)
  
  outFile <- file.path(outFolder, paste0("tad_absRankDiff_", plotvar, "_aboveThresh_randomData_density.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG*1.2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
}


######################################################################################
######################################################################################
######################################################################################
# cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


