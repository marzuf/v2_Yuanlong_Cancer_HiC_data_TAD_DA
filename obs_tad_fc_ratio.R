
options(scipen=100)

# Rscript obs_tad_fc_ratio.R

script_name <- "obs_tad_fc_ratio.R"

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

outFolder <- file.path("OBS_TAD_FC_RATIO")
dir.create(outFolder, recursive = TRUE)

runFolder <- file.path("..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA")

pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"


all_hicds <- list.files(file.path(pipFolder))
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

hicds="LG1_40kb"
exprds="TCGAlusc_norm_lusc"


if(buildData) {
  foo1 <- foreach(hicds = all_hicds) %do%{
    foo2 <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      cat(paste0("... start ", hicds, " - ", exprds, "\n"))
      
      g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2t_file))
      
      g2t_dt <- read.delim(g2t_file, sep="\t", col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE, header=FALSE)
      g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
      
      de_dt <- get(load(file.path(pipFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")))
      
      geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
      all(names(geneList) %in% de_dt$genes)
      
      
      stopifnot(geneList %in% g2t_dt$entrezID)
      g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      
      all_tads <- unique(as.character(g2t_dt$region))
      
      curr_tad = all_tads[1]
      all_obs_ratioFC <- foreach(curr_tad = all_tads, .combine='c') %dopar% {
        
        # get FC
        tad_entrez <- g2t_dt$entrezID[g2t_dt$region == curr_tad]
        stopifnot(tad_entrez %in% geneList)
        tad_entrez <- unique(tad_entrez)
        
        tad_entrez_de <- names(geneList)[geneList %in% tad_entrez]
        tad_entrez_de <- unique(tad_entrez_de)
        stopifnot(tad_entrez_de %in% de_dt$genes)
        all_tad_fc <- de_dt$logFC[de_dt$genes %in% tad_entrez_de]

        tad_negFC <- all_tad_fc[all_tad_fc < 0]
        
        sum(abs(tad_negFC))/sum(abs(all_tad_fc))
        
      } # end-foreach TAD
      names(all_obs_ratioFC) <- all_tads
      outFile <- file.path(outFolder, hicds, exprds, "all_obs_ratioFC.Rdata")
      dir.create(dirname(outFile), recursive = TRUE)
      save(all_obs_ratioFC, file=outFile, version=2)
      cat(paste0("... written: ", outFile, "\n"))
    } #end-foreach exprds  
  } #end-foreach hicds
  
} else {
  outFile <- file.path(outFolder, hicds, exprds, "all_obs_ratioFC.Rdata")
  all_obs_ratioFC <- get(load(outFile))
}


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



