########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript allTADs_and_purity.R\n"))

script_name <- "allTADs_and_purity.R"


####### !!!!!!!!!!!  FOR THE MOMENT I DO NOT RESTRICT 
# the datasets in which I look at expression (I do not ensure that is a DS where the conserved region is signif. DA)
# rationale: I want to look if this set of genes is related to purity, irrespective of DA


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))




# Rscript purity_available_samples_pkg.R 

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- T


mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")



# to quickly retrieve tad-level stat.
all_result_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
} else{
  purity_ds <- ""
}

script0_name <- "0_prepGeneData"


if(purity_ds == "") {
  file_suffix <- ""
  purity_file <- file.path("tcga_purity_aran2015.csv")
  purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
   pm <- purity_metrics[1]
   purity_plot_name <- "aran"
  # all the ranks are between 1 and 0
} else if(purity_ds == "EPIC") {
  file_suffix <- "_EPIC"
  purity_file <- file.path("../OLDER_v2_Yuanlong_Cancer_HiC_data_TAD_DA/EPIC_PURITY/all_epic_purity_data.Rdata")
  epic_purity_data <- get(load(purity_file))
  purity_dt <- as.data.frame(do.call(rbind, c(lapply(epic_purity_data, function(x) x[["infiltration_fraction"]]))))
  all_pm_metrics <- colnames(purity_dt) #"Bcells"      "CAFs"        "CD4_Tcells"  "CD8_Tcells"  "Endothelial" "Macrophages" "NKcells"     "otherCells" 
  pm <- "otherCells"
  purity_dt$Sample.ID <- rownames(purity_dt)
  purity_dt$Sample.ID <- gsub("\\.", "-", purity_dt$Sample.ID)
  purity_plot_name <- "EPIC"
} else {
  stop("--invalid purity_ds\n")
}
# 
# pm <- "CPE"
# purity_plot_name <- "aran - CPE"

library(TCGAbiolinks)
purity_dt <- Tumor.purity
purity_dt <- purity_dt[,c("Sample.ID", paste0(pm))]
purity_dt <- na.omit(purity_dt)
purity_dt$Sample.ID <- sapply(purity_dt$Sample.ID, function(x) substr(x, 1, 15))
agg_purity_dt <- aggregate(.~Sample.ID, data=purity_dt, FUN=mean)
purity_dt <- agg_purity_dt

stopifnot(!duplicated(purity_dt$Sample.ID))
purity_values <- setNames(purity_dt[,paste0(pm)], purity_dt$Sample.ID)

outFolder <- file.path("PURITY_AVAILABLE_SAMPLES_PKG", purity_ds, pm)
dir.create(outFolder, recursive = TRUE)

cat(paste0("!!! > purity metric: ", pm, "\n"))

all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

cat(paste0(">>> purity metric\t=\t", pm, "\n"))

mainFolder <- "."

# all_ds=all_ds[1]

ex_hicds <- "ENCSR489OCU_NCI-H460_40kb"
ex_exprds <- "TCGAluad_norm_luad"
ex_TAD <- "chr11_TAD390"

# all_ds = file.path(ex_hicds, ex_exprds)

if(buildTable){
  
  all_ds_avPurity_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
    
    cat(paste0("... start: ", ds, "\n"))
    
    hicds <- dirname(ds)
    exprds <- basename(ds)
    
    settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingFile))
    source(settingFile)
    
    samp1 <- get(load(file.path(setDir, sample1_file)))
    samp2 <- get(load(file.path(setDir, sample2_file)))
    
    pur_samp1 <- samp1[samp1 %in% purity_dt$Sample.ID | paste0(samp1, "A") %in% purity_dt$Sample.ID]
    cat(paste0("For ", cond1, " - available samples:\t", length(pur_samp1), "/", length(samp1), "\n"))
    
    pur_samp2 <- samp2[samp2 %in% purity_dt$Sample.ID | paste0(samp2, "A") %in% purity_dt$Sample.ID]
    cat(paste0("For ", cond2, " - available samples:\t", length(pur_samp2), "/", length(samp2), "\n"))
    
    if(length(pur_samp1) == 0 & length(pur_samp2) == 0) return(NULL)
    
    pur2_samp1 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp1  | purity_dt$Sample.ID %in% paste0(samp1, "A") ]
    stopifnot(length(pur2_samp1) == length(pur_samp1))
    
    pur2_samp2 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp2  | purity_dt$Sample.ID %in% paste0(samp2, "A") ]
    stopifnot(length(pur2_samp2) == length(pur_samp2))
    
    av_samp2 <- samp2[paste0(samp2, "A") %in% purity_dt$Sample.ID | samp2 %in% purity_dt$Sample.ID]
    match_samp2 <- sapply(av_samp2, function(x) {
      sampmatch <- purity_dt$Sample.ID[purity_dt$Sample.ID == x]
      if(length(sampmatch) == 0) {
        sampmatch <- purity_dt$Sample.ID[purity_dt$Sample.ID == paste0(x,"A")]
      }
      stopifnot(length(sampmatch) == 1)
      return(sampmatch)
    })
    missing_samp2 <- samp2[! (paste0(samp2, "A") %in% purity_dt$Sample.ID | samp2 %in% purity_dt$Sample.ID)]
    
    
    av_samp1 <- samp1[paste0(samp1, "A") %in% purity_dt$Sample.ID | samp1 %in% purity_dt$Sample.ID]
    match_samp1 <- sapply(av_samp1, function(x) {
      sampmatch <- purity_dt$Sample.ID[purity_dt$Sample.ID == x]
      if(length(sampmatch) == 0) {
        sampmatch <- purity_dt$Sample.ID[purity_dt$Sample.ID == paste0(x,"A")]
      }
      stopifnot(length(sampmatch) == 1)
      return(sampmatch)
      })
    if(length(match_samp1) == 0) match_samp1 <- c()
    missing_samp1 <- samp1[! (paste0(samp1, "A") %in% purity_dt$Sample.ID | samp1 %in% purity_dt$Sample.ID)]
    
    stopifnot(length(missing_samp1) + length(av_samp1) == length(samp1))
    stopifnot(length(missing_samp2) + length(av_samp2) == length(samp2))
    
    stopifnot(length(match_samp1) == length(av_samp1))
    stopifnot(length(match_samp2) == length(av_samp2))
    
    id_samp <-c(av_samp1, missing_samp1, av_samp2, missing_samp2)
    aran_samp <-  c(match_samp1, rep(NA, length(missing_samp1)), match_samp2, rep(NA, length(missing_samp2)))
    all_conds <- c(rep(cond1, length(samp1)), rep(cond2, length(samp2)))
    
    stopifnot(length(all_conds) == length(aran_samp))
    stopifnot(length(id_samp) == length(aran_samp))
    
    out_dt <- data.frame(
      hicds=hicds,
      exprds=exprds,
      id_samp = id_samp,
      aran_samp =aran_samp,
      cond = all_conds,
      stringsAsFactors = FALSE
    )

	out_dt$purity <- purity_values[paste0(out_dt$aran_samp)]

    out_dt

  }
  outFile <- file.path(outFolder,"all_ds_avPurity_dt.Rdata" )
  save(all_ds_avPurity_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder,"all_ds_avPurity_dt.Rdata" )
  all_ds_avPurity_dt <- get(load(outFile))
}

######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

