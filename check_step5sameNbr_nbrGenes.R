# Rscript check_step5sameNbr_nbrGenes.R

startTime <- Sys.time()

options(scipen=100)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")



script_name <- "check_step5sameNbr_nbrGenes.R"

outFolder <- "CHECK_STEP5SAMENBR_NBRGENES"
dir.create(outFolder, recursive = TRUE)

all_sampTypes <- c("sameKb", "fixKb", "sameNbr")

script0_name <- "0_prepGeneData"
script7sameNbr_name <- "7sameNbr_runPermutationsMeanTADCorr"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight 


samp_type <- "sameNbr"
fixKbSize <- ""

minGenes <- 3


pipFolder <- file.path(".")
pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")

all_gene_files <- list.files(pipOutFolder, pattern="pipeline_geneList.Rdata", recursive=TRUE, full.names=TRUE)
all_gene_files <- all_gene_files[grepl(script0_name, all_gene_files)]

all_permutCorr_files <- list.files(pipOutFolder, pattern="meanCorr_sample_around_TADs_sameNbr.Rdata", recursive=TRUE, full.names=TRUE)
all_permutCorr_files <- all_permutCorr_files[grepl(script7sameNbr_name, all_permutCorr_files)]

stopifnot(file.exists(all_permutCorr_files))
stopifnot(file.exists(all_gene_files))

nDS <- length(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), list.files)))

stopifnot(length(all_permutCorr_files) == nDS)
stopifnot(length(all_gene_files) == nDS)


### PREPARE MEAN CORR DATA
nGenes_sampCorr_DT <- foreach(c_file = all_permutCorr_files, .combine='rbind') %dopar% {
  corr_data <- get(load(c_file))
  nGenes_values <- unlist(lapply(corr_data, function(x)x[["nGenes"]]))
  data.frame(
    hicds = basename(dirname(dirname(dirname(c_file)))),
    exprds = basename(dirname(dirname(c_file))),
    region = names(nGenes_values),
    nGenes_permut = as.numeric(nGenes_values),
    stringsAsFactors = FALSE
  )
}

### PREPARE OBSERVED DATA
nGenes_DT <- foreach(c_file = all_gene_files, .combine='rbind') %dopar% {
  stopifnot(file.exists(c_file))
  geneList <- get(load(c_file))
  hicds <- basename(dirname(dirname(dirname(c_file))))
  g2tfile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
  stopifnot(file.exists(g2tfile))
  g2tdt <- read.delim(g2tfile, header=F, stringsAsFactors = F, col.names=c("entrezID", "chromo", "start", "end", "region"))
  g2tdt$entrezID <- as.character(g2tdt$entrezID)
  stopifnot(geneList %in% g2tdt$entrezID)
  g2tdt <- g2tdt[g2tdt$entrezID %in% geneList,]
  nGenes <- table(g2tdt$region)
  stopifnot(nGenes >= minGenes)
  data.frame(
    hicds = hicds,
    exprds = basename(dirname(dirname(c_file))),
    region = names(nGenes),    
    nGenes_obs = as.numeric(nGenes),
    stringsAsFactors = FALSE
  )
}

stopifnot(nrow(nGenes_sampCorr_DT) == nrow(nGenes_DT))

merged_dt <- merge(nGenes_sampCorr_DT, nGenes_DT, by=c("hicds", "exprds", "region"), all=TRUE)
stopifnot(!is.na(merged_dt))

totDS <- length(unique(paste0(merged_dt$hicds, merged_dt$exprds)))

x_var <- "nGenes_obs"
y_var <- "nGenes_permut"

myx <- merged_dt[,paste0(x_var)]
myy <- merged_dt[,paste0(y_var)]

outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_all_ds_densplot", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x=myx,
  y=myy,
  xlab = paste0(x_var),
  ylab = paste0(y_var),
  main = paste0( y_var, " vs. ",x_var)
)
mtext(side=3, text = paste0("all DS - n=", totDS))
addCorr(x = myx, y = myy, bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




