#!/usr/bin/Rscript

startTime <- Sys.time()

# Rscript familyCliques_runTADmeanCorrRatioDown.R

script_name <- "familyCliques_runTADmeanCorrRatioDown.R"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "svg"
myHeight <- 5
myWidth <- 7

# don't add at TADs the end and beginning -> I just loose half TADs, and poor quality data at extremity

require(doMC)
require(foreach)
registerDoMC(40)
require(reshape2)
require(igraph)

runFolder <- "."
pipFolder <- file.path(runFolder, "PIPELINE", "OUTPUT_FOLDER")
familyVar <- "hgnc_family_short"


withDiag <- FALSE
minCmpntSize <- 3
minGenes <- 3
maxSameTAD <- 0.5
corrMethod <- "pearson"

nMaxSize <- 1


outFolder <- file.path("FAMILYCLIQUES_RUNTADMEANCORRRATIODOWN", nMaxSize)

inFolder <- file.path("PREP_FAMILYCLIQUES", nMaxSize)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

hicds = "Barutcu_MCF-10A_40kb"
exprds="TCGAbrca_lum_bas"

buildData <- TRUE

# all_hicds=all_hicds[1:5]
# all_hicds=all_hicds[1]

if(buildData) {
  
  all_meanCorr_ratioDown <- foreach(hicds = all_hicds) %do%{
    cat(paste0("... start: ", hicds, "\n"))
    exprds_corr <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
      
      # retrieve file 
      
      famMod_file <- file.path(inFolder, hicds, exprds, "all_fams_dt.Rdata")
      
      
      stopifnot(file.exists(famMod_file))
      fam_data <- get(load(famMod_file))
      
      fam_dt <- do.call(rbind, lapply(fam_data, function(x) x[["fam_cl_dt"]]))
      
      fam_dt$entrezID <- as.character(fam_dt$entrezID)
      fam_dt$clique <- as.character(fam_dt$clique)
      
      # INPUT DATA
      gene2tadDT_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(gene2tadDT_file))
      gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
      gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
      all_gene2tadDT <- gene2tadDT
      gene2tadDT <- gene2tadDT[grepl("_TAD", gene2tadDT$region),]
      
      
      pipeline_geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
      rna_geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_geneList.Rdata")))
      norm_rnaseqDT <- eval(parse(text = load(file.path(pipFolder, hicds, exprds,
                                                        "0_prepGeneData", "rna_qqnorm_rnaseqDT.Rdata")))) 
      de_DT <-  get(load(file.path(pipFolder, hicds, exprds, "1_runGeneDE", "DE_topTable.Rdata")))
      
      stopifnot(fam_dt$entrezID %in% gene2tadDT$entrezID)  ### NB: I took only genes from TADs !!!
      stopifnot(setequal(names(rna_geneList) ,rownames(norm_rnaseqDT)))
      stopifnot(sum(fam_dt$entrezID %in% names(rna_geneList)) <= sum(fam_dt$entrezID %in% rna_geneList))
      stopifnot((rna_geneList) %in% all_gene2tadDT$entrezID)
      
      stopifnot(de_DT$genes %in% names(rna_geneList) )
      de_DT$genes2 <- rna_geneList[de_DT$genes]
      stopifnot(de_DT$genes2 %in% rna_geneList)
      stopifnot(de_DT$genes2 %in% all_gene2tadDT$entrezID) # here I have genes from TADs in de_DT
      stopifnot(!is.na(de_DT$genes2))
      stopifnot(sum(fam_dt$entrezID %in% de_DT$genes2) >= sum(fam_dt$entrezID %in% de_DT$genes))
      

      #stopifnot(rna_geneList %in% de_DT$genes2 ) WRONG: the rna_geneList not filtered for minCount
      
      # reorder !!! very important for the indexing afterwards !!!
      norm_rnaseqDT <- norm_rnaseqDT[names(rna_geneList),]
      
      
      all_famCls <- unique(fam_dt$clique)
      all_meanCorr_ratioDown_famCls <- foreach(famCl=all_famCls) %dopar% {
        
        cl_genes <- fam_dt$entrezID[as.character(fam_dt$clique) == as.character(famCl)]
        stopifnot(length(cl_genes) >= minCmpntSize)
        
        corr_cl_gene2tad_dt <- gene2tadDT[gene2tadDT$entrezID %in% cl_genes &
                                       # corrected here 14.05 -> I need to check here that I have expression data for these genes
                                       gene2tadDT$entrezID %in% rna_geneList,
                                     ]
        
        fc_cl_gene2tad_dt <- gene2tadDT[gene2tadDT$entrezID %in% cl_genes &
                                       gene2tadDT$entrezID %in% de_DT$genes2 ,
                                     ] # need to subset here for then next if keptTADs !
        
        # the DE table is filtered for minCount, not the corr data
        
        stopifnot(nrow(fc_cl_gene2tad_dt) <= nrow(corr_cl_gene2tad_dt))
        stopifnot(fc_cl_gene2tad_dt$entrezID %in% corr_cl_gene2tad_dt$entrezID)
        
        stopifnot(corr_cl_gene2tad_dt$entrezID %in% rna_geneList) 
        stopifnot(fc_cl_gene2tad_dt$entrezID %in% rna_geneList) 
        
        if(max(table(corr_cl_gene2tad_dt$region)/nrow(corr_cl_gene2tad_dt)) > maxSameTAD){
          meanCorr <- paste0("sameTAD>", maxSameTAD)
          corr_keptTADs <- corr_cl_gene2tad_dt$region
          corr_keptGenes <- corr_cl_gene2tad_dt$entrezID
          
        } else {
          rowsToKeep <- which(rna_geneList %in% cl_genes)
          stopifnot(rowsToKeep == which(rna_geneList %in% corr_cl_gene2tad_dt$entrezID))
          stopifnot(rownames(norm_rnaseqDT) == names(rna_geneList))

          
          if(length(rowsToKeep) < minGenes) {
            meanCorr <- paste0("<", minGenes, "genes")
            corr_keptTADs <- corr_cl_gene2tad_dt$region
            corr_keptGenes <- corr_cl_gene2tad_dt$entrezID
            
            
          }else{
            # this is ok because I have reordered norm_rnaseqDT : norm_rnaseqDT <- norm_rnaseqDT[names(rna_geneList),]    
            subData <- as.data.frame(t(norm_rnaseqDT[rowsToKeep,,drop=F]))
            stopifnot(colnames(subData) == names(rna_geneList[rowsToKeep]))
            # columns => the genes
            # rows => the samples
            ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
            ### UPDATE: this should not happen in the latest version !!!
            # stopifnot(ncol(subData) == length(reg_genes))
            # stopifnot(ncol(subData) == length(geneList[which(geneList %in% reg_genes)]))
            stopifnot(ncol(subData) == length(rowsToKeep))
            #### CORRELATION
            corrMatrix_all <- cor(subData, method = corrMethod)
            # should be correlation of the genes
            ######################################################## BECAUSE I COULD HAVE DUPLICATED GENE ENTREZ ID FOR DIFFERENT NAMES !!!
            ### UPDATE: this should not happen in the latest version !!!
            # stopifnot(nrow(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
            # stopifnot(ncol(corrMatrix_all) == length(geneList[which(geneList %in% reg_genes)]))
            stopifnot(nrow(corrMatrix_all) == length(rowsToKeep))
            stopifnot(ncol(corrMatrix_all) == length(rowsToKeep))
            meanCorr <- mean(corrMatrix_all[lower.tri(corrMatrix_all, diag = withDiag)], na.rm=T)
            corr_keptGenes <- rna_geneList[rowsToKeep]
            stopifnot(length(rowsToKeep) == length(corr_keptGenes))
            corr_keptTADs <- corr_cl_gene2tad_dt$region
            }

        }
        
        
        if(max(table(fc_cl_gene2tad_dt$region)/nrow(fc_cl_gene2tad_dt)) > maxSameTAD){
          ratioDown <- paste0("sameTAD>", maxSameTAD)
          fc_keptTADs <- fc_cl_gene2tad_dt$region
          fc_keptGenes <- fc_cl_gene2tad_dt$entrezID
        } else {
          stopifnot(fc_cl_gene2tad_dt$entrezID %in% de_DT$genes2)
          cl_de_DT <- de_DT[de_DT$genes2 %in% fc_cl_gene2tad_dt$entrezID,]
          stopifnot(nrow(cl_de_DT) == nrow(fc_cl_gene2tad_dt))
          if(nrow(cl_de_DT) < minGenes) {
            ratioDown <- paste0("<", minGenes, "genes")
            fc_keptTADs <- fc_cl_gene2tad_dt$region
            fc_keptGenes <- fc_cl_gene2tad_dt$entrezID
            
          } else {
            ratioDown <- sum(sign(cl_de_DT$logFC) == -1)/nrow(cl_de_DT)
            fc_keptGenes <- cl_de_DT$genes2
            
            
            fc_keptTADs <- fc_cl_gene2tad_dt$region
          }

        }
        
        

        list(
          ratioDown=ratioDown,
          meanCorr=meanCorr,
          fc_keptGenes=fc_keptGenes,
          fc_keptTADs=fc_keptTADs,
          corr_keptGenes=corr_keptGenes,
          corr_keptTADs=corr_keptTADs
        )
      }
      
      cat(paste0("... end intra-cpt correlation\n"))
      
      names(all_meanCorr_ratioDown_famCls) <- all_famCls
      stopifnot(length(all_meanCorr_ratioDown_famCls) == length(all_famCls))
      
      outFile <- file.path(outFolder, hicds, exprds, "all_meanCorr_ratioDown_famCls.Rdata")
      dir.create(dirname(outFile), recursive = TRUE)
      save(all_meanCorr_ratioDown_famCls, file= outFile, version=2)
      cat(paste0("... written: ", outFile,  "\n"))
      
      # famCorr_data <- get(load("FAMILYMODULES_RUNMEANTADCORR/Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/all_meanCorr_ratioDown_famCls.Rdata"))
      toKeep <- unlist(lapply(all_meanCorr_ratioDown_famCls, function(x) is.numeric(x [["meanCorr"]])))
      stopifnot(length(toKeep) == length(all_meanCorr_ratioDown_famCls))
      famCorr_dataF <- all_meanCorr_ratioDown_famCls[toKeep]
      
      
      famCorr <- unlist(lapply(famCorr_dataF,  function(x) x[["meanCorr"]]))
      obsCorr <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata" )))
      outFile <- file.path(outFolder, hicds, exprds, paste0(hicds, "_", exprds, "_obs_famCl_meanCorr_density.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      plot_multiDens(
        list(famClique_meanCorr = famCorr,
             obsTAD_meanCorr = obsCorr),
        plotTit = paste0(hicds, " -  ", exprds)
      )
      mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize, "; minGenes=", minGenes,  "; maxSameTAD=", maxSameTAD), font=3)
      foo <- dev.off()
      cat(paste0("... written: ", outFile,  "\n"))
      
      
      
      # famRatioDown_data <- get(load("FAMILYMODULES_RUNMEANTADCORR/Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/all_meanCorr_ratioDown_famCls.Rdata"))
      
      toKeep <- unlist(lapply(all_meanCorr_ratioDown_famCls, function(x) is.numeric(x [["ratioDown"]])))
      stopifnot(length(toKeep) == length(all_meanCorr_ratioDown_famCls))
      famRatioDown_dataF <- all_meanCorr_ratioDown_famCls[toKeep]
      famRatioDown <- unlist(lapply(famRatioDown_dataF,  function(x) x[["ratioDown"]]))
      obsRatioDown <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown", "all_obs_ratioDown.Rdata" )))
      
      
      
      outFile <- file.path(outFolder, hicds, exprds, paste0(hicds, "_", exprds, "_obs_famCl_ratioDown_density.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      plot_multiDens(
        list(famClique_ratioDown = famRatioDown,
             obsTAD_ratioDown = obsRatioDown),
        plotTit = paste0(hicds, " -  ", exprds)
      )
      mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize, "; minGenes=", minGenes,  "; maxSameTAD=", maxSameTAD), font=3)
      foo <- dev.off()
      cat(paste0("... written: ", outFile,  "\n"))
      
      list(famClique_ratioDown = famRatioDown,
           obsTAD_ratioDown = obsRatioDown,
          famClique_meanCorr = famCorr,
           obsTAD_meanCorr = obsCorr
      )
    }
    names(exprds_corr) <- all_exprds[[paste0(hicds)]]
    exprds_corr
  }
  names(all_meanCorr_ratioDown) <- all_hicds
  
  outFile <- file.path(outFolder,  "all_meanCorr_ratioDown.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(all_meanCorr_ratioDown, file= outFile, version=2)
  cat(paste0("... written: ", outFile,  "\n"))
} else {
  outFile <- file.path(outFolder,  "all_meanCorr_ratioDown.Rdata")
  all_meanCorr_ratioDown <- get(load(outFile))
}


all_fam_corr <- lapply(all_meanCorr_ratioDown, function(sublist) lapply(sublist, function(x) x[["famClique_meanCorr"]]))
all_obs_corr <- lapply(all_meanCorr_ratioDown, function(sublist) lapply(sublist, function(x) x[["obsTAD_meanCorr"]]))
nDS <- length(unlist(all_fam_corr, recursive = FALSE))


outFile <- file.path(outFolder, paste0("allDS_obs_famCl_meanCorr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(famClique_meanCorr = unlist(all_fam_corr),
       obsTAD_meanCorr = unlist(all_obs_corr)),
  my_xlab = paste0("intra-TAD/component meanCorr"),
  plotTit = paste0( "famCliques - meanCorr - all datasets - n =", nDS )
)
mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize, "; minGenes=", minGenes,  "; maxSameTAD=", maxSameTAD), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))



all_fam_ratioDown <- lapply(all_meanCorr_ratioDown, function(sublist) lapply(sublist, function(x) x[["famClique_ratioDown"]]))
all_obs_ratioDown <- lapply(all_meanCorr_ratioDown, function(sublist) lapply(sublist, function(x) x[["obsTAD_ratioDown"]]))
nDS <- length(unlist(all_fam_ratioDown, recursive = FALSE))


outFile <- file.path(outFolder, paste0("allDS_obs_famCl_ratioDown_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(famClique_ratioDown = unlist(all_fam_ratioDown),
       obsTAD_ratioDown = unlist(all_obs_ratioDown)),
  my_xlab = paste0("intra-TAD/component ratioDown"),
  plotTit = paste0( "famCliques - ratioDown - all datasets - n =", nDS )
)
mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize, "; minGenes=", minGenes,  "; maxSameTAD=", maxSameTAD), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))


cat(paste0("*** DONE: ", script_name, "\n"))


# dt1 <- get(load("FAMILYMODULES_RUNMEANTADCORR/1/all_corr.Rdata"))
# dt2 <- get(load("FAMILYMODULES_RUNMEANTADCORR/all_corr.Rdata"))
# 
# 
# 

