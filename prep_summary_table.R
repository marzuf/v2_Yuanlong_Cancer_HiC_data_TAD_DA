source("../MANUSCRIPT_FIGURES/full_dataset_names.R")

dtA <- get(load("../MANUSCRIPT_FIGURES/FIG_2/SUMMARY_PLOT/wide_dt_fig2D_supp_fig2C_signifTADs_signifGenes_topGenes_dt.Rdata"))

geneSignif <- 0.01
dtB_all <- get(load("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))
dtB1 <- aggregate(adj.P.Val ~ hicds + exprds, data=dtB_all, FUN=function(x) sum(x <= geneSignif))
colnames(dtB1)[3] <- "nSignifGenes"
dtB_t <- unique(dtB_all[,c("hicds", "exprds", "region")])
dtB2 <- aggregate(region ~ hicds + exprds, data=dtB_t, FUN=length)
colnames(dtB2)[3] <- "nTADs"
dtB <- merge(dtB1, dtB2, by=c("hicds", "exprds"))
dtB$dataset_name <- file.path(hicds_names[dtB$hicds], exprds_names[dtB$exprds])
stopifnot(setequal(dtA$dataset_name, dtB$dataset_name))


dtC <- get(load("CMP_SIGNIF_PURITYFLAGGED_FINAL_RANDOM/aran/CPE/log10/nSignif_purityTagged_obsRandom_withSamp.Rdata"))
dtC <- dtC[,c("dataset", "nSignif_obs", "nTot_obs", "nPurityFlagged_obs", "nSignifAndFlagged_obs", "nSamp1", "nSamp2", "cond1", "cond2")]
dtC <- unique(dtC)
colnames(dtC) <- gsub("_obs", "", colnames(dtC))
stopifnot(!duplicated(dtC$dataset))
dtC$dataset_name <- file.path(hicds_names[dirname(dtC$dataset)], exprds_names[basename(dtC$dataset)])
stopifnot(dtC$dataset_name %in% dtB$dataset_name)

dtBC <-  merge(dtC, dtB, by=c("dataset_name"), all=TRUE)
dtABC <- merge(dtA, dtBC, by=c("dataset_name"), all=TRUE)

stopifnot(na.omit(dtABC)$nSignifTADs == na.omit(dtABC)$nSignif)
stopifnot(na.omit(dtABC)$nTADs == na.omit(dtABC)$nTot)

dtABC$nSignif_obs <- NULL
dtABC$nTot_obs <- NULL
dtABC$dataset <- NULL

all_summary_dt <- dtABC[,c("dataset_name", "hicds", "exprds",
                           "nSamp1", "nSamp2", "cond1", "cond2",
                           "nTADs", "nSignifTADs",
                           "nPurityFlagged", "nSignifAndFlagged",
                           "nSignifGenes", "nSignifTADs_withSignifGenes", "nSignifTADs_withTopGenes")]


x_dt <- na.omit(all_summary_dt)

stopifnot(x_dt$nPurityFlagged <= x_dt$nTADs)
stopifnot(x_dt$nSignifAndFlagged <= x_dt$nSignifTADs)
stopifnot(x_dt$nSignifTADs_withSignifGenes <= x_dt$nSignifTADs)
stopifnot(x_dt$nSignifTADs_withTopGenes <= x_dt$nSignifTADs)


outFile <- "PREP_SUMMARY_TABLE/all_summary_dt.Rdata"
dir.create(dirname(outFile))
save(all_summary_dt, file=outFile, version=2)


