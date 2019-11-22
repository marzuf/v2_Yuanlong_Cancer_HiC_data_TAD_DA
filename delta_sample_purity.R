########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript delta_sample_purity.R\n"))

script_name <- "delta_sample_purity.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript delta_sample_purity.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

outFolder <- file.path("DELTA_SAMPLE_PURITY")
dir.create(outFolder, recursive = TRUE)

purity_file <- file.path("tcga_purity_aran2015.csv")
purity_dt <- read.delim(purity_file, header=TRUE, sep="\t")
purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
pm <- purity_metrics[1]
# all the ranks are between 1 and 0

all_hicds <- list.files(pipFolder)
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL

purity_delta_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  
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
  
  all_pm_deltaValues <- sapply(purity_metrics, function(pm) {
    if(cond1 == "norm") {
      stopifnot(length(pur_samp1) == 0)
      purValues_cond1 <- 1
    } else {
      purValues_cond1 <- purity_dt[purity_dt$Sample.ID %in% pur_samp1 | purity_dt$Sample.ID %in% paste0(pur_samp1, "A"),paste0(pm)]  
    }
    stopifnot(cond2 != "norm")
    # stopifnot(length(pur_samp2) > 0)
    purValues_cond2 <- purity_dt[purity_dt$Sample.ID %in% pur_samp2 | purity_dt$Sample.ID %in% paste0(pur_samp2, "A"),paste0(pm)]
    abs(mean(purValues_cond1, na.rm=TRUE)-mean(purValues_cond2, na.rm=TRUE))
  })

  out_dt <- t(data.frame(all_pm_deltaValues))
  rownames(out_dt) <- ds
  out_dt
  
}
outFile <- file.path(outFolder, "purity_delta_dt.Rdata")
save(purity_delta_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

