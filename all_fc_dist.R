
options(scipen=100)

# Rscript all_fc_dist.R

script_name <- "all_fc_dist.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

buildData <- TRUE

require(flux)
require(foreach)
require(doMC)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ggsci)
registerDoMC(40)

outFolder <- file.path("ALL_FC_DIST")
dir.create(outFolder, recursive = TRUE)

runFolder <- file.path("..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA")

pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"


plotCex <- 1.4
plotType <- "png"
myHeight <- 400
myWidth <- 600



all_hicds <- list.files(file.path(pipFolder))
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

hicds="LG1_40kb"
exprds="TCGAlusc_norm_lusc"


if(buildData) {
  all_fc_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    ds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat(paste0("... start ", hicds, " - ", exprds, "\n"))
      

            
      de_dt <- get(load(file.path(pipFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")))
      
      geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
      
      
      
      stopifnot(names(geneList) %in% de_dt$genes)
      all_fc <- setNames(de_dt$logFC[de_dt$genes %in% names(geneList)], de_dt$genes[de_dt$genes %in% names(geneList)])
      
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        entrez_DE = names(all_fc),
        FC = as.numeric(all_fc),
        stringsAsFactors = FALSE
      )

    } #end-foreach exprds  
    ds_dt
    
  } #end-foreach hicds
  outFile <- file.path(outFolder,  "all_fc_dt.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(all_fc_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_fc_dt.Rdata")
  all_fc_dt <- get(load(outFile))
}


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
# source("cancer_cols.R")
all_fc_dt$cmp_type <- as.character(all_cmps[paste0(all_fc_dt$exprds)])
stopifnot(!is.na(all_fc_dt$cmp_type))

ratioDown <- c(by(all_fc_dt, all_fc_dt$cmp_type, function(dt) sum(dt$FC <0)/nrow(dt)))
ratioFC <- c(by(all_fc_dt, all_fc_dt$cmp_type, function(dt) sum(abs(dt$FC[dt$FC <0]))/sum(abs(dt$FC))))

outFile <- file.path(outFolder, paste0("geneFC_dist_allDS_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

par(bty="L")
plot_multiDens(
  split(x = all_fc_dt$FC, all_fc_dt$cmp_type)
)
legend("topleft",
       legend = paste0(names(ratioDown), "=", round(ratioDown,4)),
       bty="n",
       title.adj = 0.5,
       title = "ratioDown:")

legend("bottomleft",
       legend = paste0(names(ratioFC), "=", round(ratioFC,4)),
       bty="n",
       title.adj = 0.5,
       title = "ratioFC:")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


