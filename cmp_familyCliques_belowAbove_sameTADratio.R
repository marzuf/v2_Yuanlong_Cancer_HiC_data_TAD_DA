# Rscript cmp_familyCliques_belowAbove_sameTADratio.R

inFolderSame <- "FAMILYCLIQUES_RUNTADMEANCORRRATIODOWN_SAMETAD"
inFolderDiff <- "FAMILYCLIQUES_RUNTADMEANCORRRATIODOWN"
nMaxTADsize <- 2
maxSameTAD <- 0.5

plotType <- "svg"
myHeight <- 5
myWidth <- 7

# "sameTAD": 
# discard if  if(max(table(corr_cpt_gene2tad_dt$region)/nrow(corr_cpt_gene2tad_dt)) <= maxSameTAD) # keep if > 0.5 from sameTAD


diffTAD_data <- get(load(file.path(inFolderDiff, nMaxTADsize, "all_meanCorr_ratioDown.Rdata")))
sameTAD_data <- get(load(file.path(inFolderSame, nMaxTADsize, "all_meanCorr_ratioDown.Rdata")))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- file.path("CMP_FAMILYCLIQUES_BELOWABOVE_SAMETADRATIO", nMaxTADsize)
dir.create(outFolder, recursive = TRUE)

sameTAD_rdCmpnt <- lapply(sameTAD_data, function(sublist) lapply(sublist, function(x) x[["famClique_ratioDown"]]))
diffTAD_rdCmpnt <- lapply(diffTAD_data, function(sublist) lapply(sublist, function(x) x[["famClique_ratioDown"]]))

plot_list <- list(
  sameTAD=unlist(sameTAD_rdCmpnt),
  diffTAD=unlist(diffTAD_rdCmpnt)
)

outFile <- file.path(outFolder, paste0("sameTAD_diffTAD_ratioDown_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  plot_list,
  plotTit = paste0("TAD ratioDown (family components; max. size=", nMaxTADsize, "*mean TAD size)"),
  legPos = "topleft"
)
mtext(side=3, text=paste0("sameTAD: > ", maxSameTAD, " same TAD; diffTAD: <=", maxSameTAD, " sameTAD"), font=3)
foo <- dev.off()



sameTAD_corrCmpnt <- lapply(sameTAD_data, function(sublist) lapply(sublist, function(x) x[["famClique_meanCorr"]]))
diffTAD_corrCmpnt <- lapply(diffTAD_data, function(sublist) lapply(sublist, function(x) x[["famClique_meanCorr"]]))

plot_list <- list(
  sameTAD=unlist(sameTAD_corrCmpnt),
  diffTAD=unlist(diffTAD_corrCmpnt)
)

outFile <- file.path(outFolder, paste0("sameTAD_diffTAD_meanCorr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  plot_list,
  plotTit = paste0("TAD meanCorr (family cliques; max. size=", nMaxTADsize, "*mean TAD size)")
)
mtext(side=3, text=paste0("sameTAD: > ", maxSameTAD, " same TAD; diffTAD: <=", maxSameTAD, " sameTAD"), font=3)
foo <- dev.off()


