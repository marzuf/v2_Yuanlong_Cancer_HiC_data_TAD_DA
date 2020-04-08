
options(scipen=100)

# Rscript permG2t_tad_fc_ratio.R

script_name <- "permG2t_tad_fc_ratio.R"

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

outFolder <- file.path("PERMG2T_TAD_FC_RATIO")
dir.create(outFolder, recursive = TRUE)

runFolder <- file.path("..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA")

pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script5_name <- "5_runPermutationsMedian"

nPermut <- 1000


all_hicds <- list.files(file.path(pipFolder))
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

# hicds_toplot <- c("ENCSR504OTV_transverse_colon_RANDOMMIDPOSSTRICT_40kb","ENCSR312KHQ_SK-MEL-5_RANDOMMIDPOSSTRICT_40kb")
# exprds_toplot <- c("TCGAcoad_msi_mss", "TCGAskcm_wt_mutCTNNB1")
hicds= "ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAluad_mutKRAS_mutEGFR"
hicds="ENCSR504OTV_transverse_colon_40kb"
exprds="TCGAcoad_msi_mss"
hicds="ENCSR312KHQ_SK-MEL-5_40kb"
exprds="TCGAskcm_wt_mutCTNNB1"
myhicds=hicds
myexprds=exprds

if(buildData) {
  foo1 <- foreach(hicds = all_hicds) %do%{
    foo2 <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      if(hicds != myhicds | exprds != myexprds) return(NULL)
      
      cat(paste0("... start ", hicds, " - ", exprds, "\n"))
      
      cat(paste0("... load permutation data\n"))
      permut_dt <-  get(load(file.path(pipFolder, hicds, exprds, script5_name, "permutationsDT.Rdata")))
      # permut_dt <-  get(load(paste0(hicds, "_", exprds, "_1000Permut_permDT.Rdata")))
      cat(paste0("... loaded\n"))
      
      geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
      stopifnot(setequal(geneList, rownames(permut_dt)))

      
      de_dt <- get(load(file.path(pipFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")))

      all(names(geneList) %in% de_dt$genes)      
      
      permut_dt <- permut_dt[,1:nPermut]
      
      all_tads <- sort(unique(permut_dt[,1]))
      
      permut_ratio_fc_dt <- foreach(i = 1:ncol(permut_dt), .combine = 'cbind') %dopar% {
        
        g2t_dt <- data.frame(
          entrezID = as.character(rownames(permut_dt)),
          region = as.character(permut_dt[,i]),
          stringsAsFactors = FALSE
          )
        
        stopifnot(setequal(unique(g2t_dt$region), all_tads))
        
        
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
        
        all_obs_ratioFC
        
      }
      checkRN <- rownames(permut_ratio_fc_dt)
      rownames(permut_ratio_fc_dt) <- all_tads
      # colnames(permut_ratio_fc_dt) <- paste0("result.", 1:ncol(permut_ratio_fc_dt))
      
      ratio_permDT <- permut_ratio_fc_dt
      
 
    
      outFile <- file.path(outFolder, hicds, exprds, paste0("ratioFC_", nPermut, "Permut_permDT.Rdata"))
      dir.create(dirname(outFile), recursive = TRUE)
      save(ratio_permDT, file=outFile, version=2)
      cat(paste0("... written: ", outFile, "\n"))
      
      stopifnot(checkRN==rownames(permut_ratio_fc_dt))
      
    } #end-foreach exprds  
  } #end-foreach hicds
  
} else {
  outFile <- file.path(outFolder, hicds, exprds, paste0("ratioFC_", nPermut, "Permut_permDT.Rdata"))
  ratio_permDT <- get(load(outFile))
}


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



