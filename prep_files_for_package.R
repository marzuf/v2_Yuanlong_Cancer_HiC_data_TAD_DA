# Rscript prep_files_for_package.R

# cp /mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/6_runPermutationsMeanLogFC/meanLogFC_permDT.Rdata .                           
# mv meanLogFC_permDT.Rdata ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_meanLogFC_permDT.Rdata
# -> but this is too big file for the package (so slow to build and load...)

require(foreach)
require(doMC)
registerDoMC(40)

outFolder <- "/mnt/etemp/marie/MANUSCRIPT_FIGURES/COCODATA/data"

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"
nPermCut <- 1000

##############################################################################################
cat(paste0("... prep corr values\n"))
corr_type <- "meanCorr"
script7sameNbr_name <- "7sameNbr_runPermutationsMeanTADCorr"

### RETRIEVE ALL THE FILES IN THE FOLDER !!!
mainPipFold <- file.path("PIPELINE", "OUTPUT_FOLDER")
all_sampleCorr_files <- list.files(mainPipFold, pattern="meanCorr_sample_around_TADs_sameNbr.Rdata", full.names = TRUE, recursive = TRUE)

all_sampleCorr_files <- all_sampleCorr_files[grepl(script7sameNbr_name, all_sampleCorr_files)]  ### added 26.11.2019 otherwise match also the files for partial corr.
all_sampleCorr_files <- all_sampleCorr_files[!grepl("RANDOM", all_sampleCorr_files) & !grepl("PERMUT", all_sampleCorr_files)]

all_hicds <- list.files(mainPipFold)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(mainPipFold, x)))
stopifnot(length(all_sampleCorr_files) == length(unlist(all_exprds)))

# in the script, use combine='c' but for the package keep as list
all_sample_corrValues <- foreach(corr_file = all_sampleCorr_files) %dopar% {
  stopifnot(file.exists(corr_file))
  corr_data <- eval(parse(text = load(corr_file)))
  all_samp_corrs <- as.numeric(sapply(corr_data, function(x) x[[paste0(corr_type)]]))
  stopifnot(!is.null(all_samp_corrs))
  all_samp_corrs <- na.omit(all_samp_corrs)  
  all_samp_corrs
}

outFile <- file.path(outFolder, "all_sample_corrValues.Rdata")
save(all_sample_corrValues, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

##############################################################################################
cat(paste0("... load permut DT\n"))
permutationsDT <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "5_runPermutationsMedian/permutationsDT.Rdata")))
permutationsDT_cut <- permutationsDT[,1:nPermCut]
permutationsDT <- permutationsDT_cut
outFile <- file.path(outFolder, paste0("cut", nPermCut, "_", hicds, "_", exprds, "_permutationsDT.Rdata"))
save(permutationsDT, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


##############################################################################################
cat(paste0("... load logFC permut DT\n"))
meanLogFC_dt <- get(load(file.path("PIPELINE/OUTPUT_FOLDER", hicds, exprds, "6_runPermutationsMeanLogFC/meanLogFC_permDT.Rdata")))
meanLogFC_dt_cut <- meanLogFC_dt[,1:nPermCut]
meanLogFC_permDT <- meanLogFC_dt_cut
outFile <- file.path(outFolder, paste0("cut", nPermCut, "_", hicds, "_", exprds, "_meanLogFC_permDT.Rdata"))
save(meanLogFC_permDT, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

