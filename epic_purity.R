########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript epic_purity.R\n"))

script_name <- "epic_purity.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

library(EPIC)

# Rscript epic_purity.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- TRUE

myHeightGG <- 7
myWidthGG <- 9
plotType <- "png"

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

outFolder <- file.path("EPIC_PURITY")
dir.create(outFolder, recursive = TRUE)


all_hicds <- list.files(pipFolder)
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

setDir <- ifelse(SSHFS, "/media/electron", "")
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


if(buildTable) {
  
  all_epic_purity_data <- foreach(ds = all_ds) %dopar% {
    
    
    hicds <- dirname(ds)
    exprds <- basename(ds)
    
    cat(paste0("... start\t", hicds, "-", exprds, "\n"))
    
    
    settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingFile))
    source(settingFile)
    
    samp1 <- get(load(file.path(setDir, sample1_file)))
    samp2 <- get(load(file.path(setDir, sample2_file)))
    
    fpkm_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")
    stopifnot(file.exists(fpkm_file))
    fpkm_dt <- get(load(fpkm_file))  
    
    # no ID conversion supported
    cat(paste0("... retain genes with available symbols:\t", nrow(fpkm_dt), " -> "))
    fpkm_dt <- fpkm_dt[rownames(fpkm_dt) %in% names(entrez2symb),]
    cat(paste0(nrow(fpkm_dt), "\n"))
    
    newnames <- entrez2symb[rownames(fpkm_dt)]
    newnames_dup <- newnames[duplicated(newnames)]
    newnames_nodup <- newnames[!newnames %in% newnames_dup]
    
    cat(paste0("... retain only no duplicated genes:\t", nrow(fpkm_dt), " -> "))
    fpkm_dt <- fpkm_dt[rownames(fpkm_dt) %in% names(newnames_nodup),]
    cat(paste0(nrow(fpkm_dt), "\n"))
    
    rownames(fpkm_dt) <- entrez2symb[rownames(fpkm_dt)]
    stopifnot(setequal(rownames(fpkm_dt), newnames_nodup))
    
    
    # tumor purity = fraction of tumor cells
    # in EPIC: = fraction of "otherCells"
    cat(paste0("... start\t", hicds, "-", exprds, "\t EPIC \n"))
    out_epic <- EPIC(bulk = fpkm_dt[,c(samp1,samp2)])
    list(
      infiltration_fraction = out_epic[["mRNAProportions"]]
    )
  } # end-foreach iterating over DS  
  names(all_epic_purity_data) <- all_ds
  outFile <- file.path(outFolder, "all_epic_purity_data.Rdata")
  save(all_epic_purity_data, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  } else {
  outFile <- file.path(outFolder, "all_epic_purity_data.Rdata")
  all_epic_purity_data <- get(load(outFile))
}
  
    
    





##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))
