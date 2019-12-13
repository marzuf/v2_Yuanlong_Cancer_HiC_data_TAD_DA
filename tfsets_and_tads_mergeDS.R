startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript tfsets_and_tads_mergeDS.R crisp 
# Rscript tfsets_and_tads_mergeDS.R c3.mir
# Rscript tfsets_and_tads_mergeDS.R c3.tft
# Rscript tfsets_and_tads_mergeDS.R c3.all
# Rscript tfsets_and_tads_mergeDS.R trrust
# Rscript tfsets_and_tads_mergeDS.R tftg
# Rscript tfsets_and_tads_mergeDS.R motifmap 


plotType <- "png"
myHeight <- 400
myWidth <- 600
plotCex <- 1.4


dsIn <- "crisp"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1 | length(args) == 2)
dsIn <- args[1]
if(length(args) == 2) {
  all_hicds <- args[2]
} else {
  all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
}

stopifnot(dsIn %in% c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap"))

outFolder <- file.path(paste0("TFSETS_AND_TADS_MERGEDS_", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)

inFolder <- file.path(paste0("TFSETS_AND_TADS_", toupper(dsIn)))

all_obs_files <- list.files(inFolder, pattern="obs_dt.Rdata", recursive = TRUE, full.names = TRUE)
all_permut_files <- list.files(inFolder, pattern="permut_dt.Rdata", recursive = TRUE, full.names = TRUE)

cat(paste0("... length(all_obs_files)\t=\t", length(all_obs_files), "\n"))
cat(paste0("... length(all_permut_files)\t=\t", length(all_permut_files), "\n"))

nDS <- length(all_permut_files)
stopifnot(length(all_obs_files) == nDS)

# build the observed data

all_obs_data <- foreach(obs_file = all_obs_files, .combine='rbind') %dopar%{
  obs_dt <- get(load(obs_file))
  obs_dt
}

all_permut_data <- foreach(permut_file = all_permut_files, .combine='rbind') %dopar%{
  permut_dt <- get(load(permut_file))
  permut_dt$TAD_ratio <- permut_dt$reg_nTADs_permut/permut_dt$reg_nGenes
  # med_permut_dt <- aggregate( .~reg_nGenes, data = permut_dt, FUN=median)
  med_permut_dt <- aggregate( .~reg_nGenes, data = permut_dt, FUN=median)
  colnames(med_permut_dt)[colnames(med_permut_dt) != "reg_nGenes"] <- paste0(colnames(med_permut_dt)[colnames(med_permut_dt) != "reg_nGenes"], "_med")
  med_permut_dt
}

outFile <- file.path(outFolder, "all_obs_data.Rdata")
save(all_obs_data, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "all_permut_data.Rdata")
save(all_permut_data, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


all_permut_data <- all_permut_data[order(all_permut_data$reg_nGenes),]



outFile <- file.path(outFolder, paste0("all_ds_med_obs_permut_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(reg_nTADs_permut_med ~ reg_nGenes, data = all_permut_data, 
        xlab = "# regulated genes ",
        ylab = "# of TADs",
        main = paste0("all DS"), 
        at = unique(all_permut_data$reg_nGenes))
points(
  x = all_obs_data$reg_nGenes,
  y = all_obs_data$reg_nTADs,
  cex.lab = plotCex,
  cex.axis = plotCex,
  pch = 16,
  col ="red"
)
legend("topleft", 
       legend=c("obs.", "permut."),
       col = c("red", "black"),
       pch=16,
       bty="n")

mtext(side=3, text = paste0("# DS = ", nDS))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("all_ds_med_obs_permut_ratio_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(TAD_ratio_med ~ reg_nGenes, data = all_permut_data, 
        xlab = "# regulated genes ",
        ylab = "# TADs/# genes",
        main = paste0("all DS"), 
        at = unique(all_permut_data$reg_nGenes))
points(
  x = all_obs_data$reg_nGenes,
  y = all_obs_data$reg_nTADs/all_obs_data$reg_nGenes,
  cex.lab = plotCex,
  cex.axis = plotCex,
  pch = 16,
  col ="red"
)
legend("topleft", 
       legend=c("obs.", "permut."),
       col = c("red", "black"),
       pch=16,
       bty="n")

mtext(side=3, text = paste0("# DS = ", nDS))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))





outFile <- file.path(outFolder, paste0("all_ds_med_obs_permut_multidens.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot_multiDens(
  list(obs = all_obs_data$reg_nTADs_obs/all_obs_data$reg_nGenes,
       permut = all_permut_data$TAD_ratio_med),
  plotTit = paste0("all DS")
)
mtext(side=3, text = paste0("# DS = ", nDS))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
  
  
# the same, but only for TF that have <= 100 targets  

outFile <- file.path(outFolder, paste0("all_ds_med_obs_permut_boxplot_to100.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(reg_nTADs_permut_med ~ reg_nGenes, data = all_permut_data[all_permut_data$reg_nGenes <= 100,], 
        xlab = "# regulated genes ",
        ylab = "# of TADs",
        main = paste0("all DS"), 
        at = unique(all_permut_data$reg_nGenes[all_permut_data$reg_nGenes <= 100]))
points(
  x = all_obs_data$reg_nGenes,
  y = all_obs_data$reg_nTADs,
  cex.lab = plotCex,
  cex.axis = plotCex,
  pch = 16,
  col ="red"
)
mtext(paste0("reg_nGenes <= 100; # DS = ", nDS), side=3)
legend("topleft", 
       legend=c("obs.", "permut."),
       col = c("red", "black"),
       pch=16,
       bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



outFile <- file.path(outFolder, paste0("all_ds_med_obs_permut_ratio_boxplot_to100.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(TAD_ratio_med ~ reg_nGenes, data = all_permut_data[all_permut_data$reg_nGenes <= 100,], 
        xlab = "# regulated genes ",
        ylab = "# TADs/# genes",
        main = paste0("all DS"), 
        at = unique(all_permut_data$reg_nGenes[all_permut_data$reg_nGenes <= 100]))
points(
  x = all_obs_data$reg_nGenes,
  y = all_obs_data$reg_nTADs_obs/all_obs_data$reg_nGenes,
  cex.lab = plotCex,
  cex.axis = plotCex,
  pch = 16,
  col ="red"
)
mtext(paste0("reg_nGenes <= 100; # DS = ", nDS), side=3)
legend("topleft", 
       legend=c("obs.", "permut."),
       col = c("red", "black"),
       pch=16,
       bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("obs_permut_multidens_to100.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot_multiDens(
  list(obs = all_obs_data$reg_nTADs_obs[all_obs_data$reg_nGenes <= 100]/all_obs_data$reg_nGenes[all_obs_data$reg_nGenes <= 100],
       permut = all_permut_data$TAD_ratio_med[all_permut_data$reg_nGenes <= 100]),
  plotTit = paste0("all DS")
)
mtext(paste0("reg_nGenes <= 100; # DS = ", nDS), side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))
