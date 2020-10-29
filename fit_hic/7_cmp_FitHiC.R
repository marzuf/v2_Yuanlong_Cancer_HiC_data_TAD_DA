

# Rscript 7_cmp_FitHiC.R RUN_FITHIC/KARPAS_DMSO/high_resol RUN_FITHIC/KARPAS_DMSO/low_resol CMP_FITHIC/KARPAS_DMSO/low_vs_high

require(foreach)
require(doMC)
registerDoMC(40)

options(scipen=100)

script_name <- "7_cmp_FitHiC.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")


signifThresh <- 0.05
signifCol <- "q_value"


outfolder1 <- "RUN_FITHIC/KARPAS_DMSO/high_resol"
outfolder2 <- "RUN_FITHIC/KARPAS_DMSO/low_resol"
#/chr21/.spline_pass2.significances.txt.gz"
output_folder <- file.path("CMP_FITHIC", "KARPAS_DMSO", "low_vs_high")

fitHiCsuffix <- ".spline_pass2.significances.txt.gz"

args <- commandArgs(trailingOnly = TRUE)
outfolder1 <- args[1]
outfolder2 <- args[2]
output_folder <- args[3]
stopifnot(dir.exists(outfolder1))
stopifnot(dir.exists(outfolder2))

dir.create(output_folder, recursive = TRUE)

logFile <- file.path(output_folder, "cmp_FitHiC_logfile.txt")
file.remove(logFile)

all_files1 <- list.files(outfolder1)
all_files2 <- list.files(outfolder2)

commonChromo <- intersect(all_files1, all_files2)
stopifnot(length(commonChromo) > 0)

txt <- paste0(" !!! HARD-CODED !!! \n")
cat(txt)
cat(txt, file=logFile,append=TRUE)
txt <- paste0("... fitHiCsuffix\t:\t", fitHiCsuffix, "\n")
cat(txt)
cat(txt, file=logFile,append=TRUE)
txt <- paste0("... signifThresh\t:\t", signifThresh, "\n")
cat(txt)
cat(txt, file=logFile,append=TRUE)
txt <- paste0("... signifCol\t:\t",signifCol, "\n")
cat(txt)
cat(txt, file=logFile,append=TRUE)


txt <- paste0("... # chromo found in outfolder1\t:\t", length(all_files1), "\n")
cat(txt)
cat(txt, file=logFile,append=TRUE)
txt <- paste0("... # chromo found in outfolder2\t:\t", length(all_files2), "\n")
cat(txt)
cat(txt, file=logFile,append=TRUE)
txt <- paste0("... # common chromo found\t:\t", length(commonChromo), "\n")
cat(txt)
cat(txt, file=logFile,append=TRUE)

########################################################################################################
########################################################################################################
########################################################################################################
chromo = commonChromo[1]
all_chr_dt <- foreach(chromo = commonChromo, .combine='rbind') %dopar% {
  
  file1 <- file.path(outfolder1,  chromo, fitHiCsuffix)
  cat(file1, "\n")
  stopifnot(file.exists(file1))
  
  file2 <- file.path(outfolder2,  chromo, fitHiCsuffix)
  cat(file1, "\n")
  stopifnot(file.exists(file2))
  
  
  fithic1 <- read.table(gzfile(file1), header=T)
  fithic2 <- read.table(gzfile(file2), header=T)
  
  stopifnot(signifCol %in% colnames(fithic1))
  stopifnot(signifCol %in% colnames(fithic2))
  
  stopifnot(fithic1$chr1 == fithic1$chr2)
  stopifnot(fithic2$chr1 == fithic2$chr1)
  stopifnot(fithic1$chr1 == chromo)
  stopifnot(fithic2$chr1 == chromo)  
  
  fithic1$interaction <- paste0(fithic1$fragmentMid1, "_", fithic1$fragmentMid2)
  fithic2$interaction <- paste0(fithic2$fragmentMid1, "_", fithic2$fragmentMid2)
  
  stopifnot(is.numeric(fithic1[,paste0(signifCol)]))
  signif1_interactions <- fithic1$interaction[fithic1[,paste0(signifCol)] <= signifThresh]
  signif2_interactions <- fithic2$interaction[fithic2[,paste0(signifCol)] <= signifThresh]
  
  commonInteractions <- intersect(fithic1$interaction, fithic2$interaction)
  commonSignifInter <- intersect(signif1_interactions, signif2_interactions)
  
  data.frame(
    chromo = chromo,
    nInteractions_1 = nrow(fithic1),
    nInteractions_2 = nrow(fithic2),
    nInteractions_12 = length(commonInteractions),
    nSignifInter_1 = length(signif1_interactions),
    nSignifInter_2 = length(signif2_interactions),
    nSignifInter_12 = length(signif1_interactions),
    stringsAsFactors = FALSE
  )

} # end-for iterating over commmon chromo

outFile <- file.path(output_folder, "all_chr_dt.txt")
write.table(all_chr_dt, file = outFile, sep="\t", quote=FALSE, append=FALSE)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(output_folder, "all_chr_dt.Rdata")
save(all_chr_dt, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

cat(paste0("... written: ", logFile, "\n"))

########################################################################################################
########################################################################################################
########################################################################################################
cat("*** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
