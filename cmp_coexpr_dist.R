# Rscript cmp_coexpr_dist.R

startTime <- Sys.time()

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

script_name <- "cmp_coexpr_dist.R"

outFolder <- "CMP_COEXPR_DIST"
dir.create(outFolder)

all_sampTypes <- c("sameKb", "fixKb", "sameNbr")


script4_name <- "4_runMeanTADCorr"
script7sameNbr_name <- "7sameNbr_runPermutationsMeanTADCorr"


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight * 1.2


samp_type <- "sameNbr"
fixKbSize <- ""


pipFolder <- file.path(".")
pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")

all_obsCorr_files <- list.files(pipOutFolder, pattern="all_meanCorr_TAD.Rdata", recursive=TRUE, full.names=TRUE)
all_obsCorr_files <- all_obsCorr_files[grepl(script4_name, all_obsCorr_files)]

all_permutCorr_files <- list.files(pipOutFolder, pattern="meanCorr_sample_around_TADs_sameNbr.Rdata", recursive=TRUE, full.names=TRUE)
all_permutCorr_files <- all_permutCorr_files[grepl(script7sameNbr_name, all_permutCorr_files)]

stopifnot(file.exists(all_permutCorr_files))
stopifnot(file.exists(all_obsCorr_files))

nDS <- length(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), list.files)))

stopifnot(length(all_permutCorr_files) == nDS)
stopifnot(length(all_obsCorr_files) == nDS)


### PREPARE MEAN CORR DATA
sampCorr_DT <- foreach(c_file = all_permutCorr_files, .combine='rbind') %dopar% {

  corr_data <- get(load(c_file))
  
  sampCorr_values <- unlist(lapply(corr_data, function(x)x[["meanCorr"]]))

  data.frame(
    meanCorr=as.numeric(sampCorr_values),
    meanCorr_type="sample_sameNbr",
    stringsAsFactors = FALSE
  )
}

### PREPARE OBSERVED DATA
meanCorr_DT <- foreach(c_file = all_obsCorr_files, .combine='rbind') %dopar% {
  stopifnot(file.exists(c_file))
  tad_corr <- eval(parse(text = load(c_file)))
  data.frame(
#    region = names(tad_corr),
    meanCorr = tad_corr,
    meanCorr_type = "observed",
    stringsAsFactors = FALSE
  )
}


all_corr_DT <- rbind(meanCorr_DT, sampCorr_DT)

# outFile <- file.path(outFolder, paste0("all_corr_DT.Rdata"))
# save(all_corr_DT, file = outFile)
# cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("meanCorr_cmp_multidens.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  split(all_corr_DT$meanCorr, all_corr_DT$meanCorr_type),
  plotTit = paste0("meanCorr comparison"),
  my_xlab = paste0("meanCorr")
  )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




