options(scipen=100)

SSHFS=F

buildData <- TRUE

# Rscript rankDiff_activDiff.R <hicds_norm> <hicds_tumor> <exprds>
# Rscript rankDiff_activDiff.R LI_40kb GSE105381_HepG2_40kb TCGAlihc_norm_lihc


script_name <- "rankDiff_activDiff.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(foreach)
require(doMC)
registerDoMC(40)


hicds_norm <- "LI_40kb"
hicds_tumor <- "GSE105381_HepG2_40kb"
exprds <- "TCGAlihc_norm_lihc"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

script0_name <- "0_prepGeneData"
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


inFile <- "TAD_MATCHING_ACROSS_HICDS/all_matching_dt.Rdata"
matching_dt <- get(load(inFile))

outFolder <- file.path("RANKDIFF_ACTIVDIFF")
dir.create(outFolder, recursive = TRUE)

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
  data.frame(
    chromo=curr_chromo,
    region=curr_region,
    start=curr_TADstart,
    end=curr_TADend,
    regionRank = curr_TADvalue,
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
matching_dt$matching_exprds <- basename(matching_dt$ref_dataset)


norm_matching_dt <- matching_dt[matching_dt$ref_exprds == exprds & 
                                  matching_dt$ref_hicds == hicds_norm & 
                                  matching_dt$matching_exprds == exprds &
                                  matching_dt$matching_hicds == hicds_tumor,]

norm_matching_pval_dt <- merge(norm_matching_dt, norm_TAD_adjPvals_dt, by="refID")




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
  data.frame(
    chromo=curr_chromo,
    region=curr_region,
    start=curr_TADstart,
    end=curr_TADend,
    regionRank = curr_TADvalue,
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

matchingCol <- "matchingID_maxOverlapBp"

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

outFile <- file.path(outFolder, paste0(hicds_norm, "_withMatching_", hicds_tumor, "_rankDiff_vs_pval_", exprds, "_densplot.", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = -log10(norm_matching_pval_tadRank_dt$adjPval),
  y = norm_matching_pval_tadRank_dt$rankDiff,
  main=paste0(exprds),
  sub=paste0("norm as refDS"),
  xlab=paste0("-log10 TAD adj. pval"),
  ylab=paste0("best matching TAD rank diff."),
  cex.axis=axisCex,
  cex.lab=axisCex
)
abline(h=0, lty=2, col="grey")
mtext(side=3, paste0(hicds_norm, " matching ", hicds_tumor))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


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

outFile <- file.path(outFolder, paste0(hicds_tumor, "_withMatching_", hicds_norm, "_rankDiff_vs_pval_", exprds, "_densplot.", plotType ))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = -log10(tumor_matching_pval_tadRank_dt$adjPval),
  y = tumor_matching_pval_tadRank_dt$rankDiff,
  main=paste0(exprds),
  sub=paste0("tumor as refDS"),
  xlab=paste0("-log10 TAD adj. pval"),
  ylab=paste0("best matching TAD rank diff."),
  cex.axis=axisCex,
  cex.lab=axisCex
)
abline(h=0, lty=2, col="grey")
mtext(side=3, paste0(hicds_tumor, " matching ", hicds_norm))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


