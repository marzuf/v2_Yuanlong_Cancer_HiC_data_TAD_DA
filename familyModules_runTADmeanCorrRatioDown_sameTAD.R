#!/usr/bin/Rscript

startTime <- Sys.time()

# Rscript familyModules_runTADmeanCorrRatioDown_sameTAD.R

script_name <- "familyModules_runTADmeanCorrRatioDown_sameTAD.R"

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

nMaxSize <- 2

outFolder <- file.path("FAMILYMODULES_RUNTADMEANCORRRATIODOWN_SAMETAD", nMaxSize)

inFolder <- file.path("PREP_FAMILYMODULES", nMaxSize)

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
      
      fam_dt <- do.call(rbind, lapply(fam_data, function(x) x[["fam_cpt_dt"]]))
      
      fam_dt$entrezID <- as.character(fam_dt$entrezID)
      fam_dt$cpt <- as.character(fam_dt$cpt)
      
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
      
      
      all_famCpts <- unique(fam_dt$cpt)
      all_meanCorr_ratioDown_famCpts <- foreach(famCpt=all_famCpts) %dopar% {
        
        cpt_genes <- fam_dt$entrezID[as.character(fam_dt$cpt) == as.character(famCpt)]
        stopifnot(length(cpt_genes) >= minCmpntSize)
        
        corr_cpt_gene2tad_dt <- gene2tadDT[gene2tadDT$entrezID %in% cpt_genes &
                                       # corrected here 14.05 -> I need to check here that I have expression data for these genes
                                       gene2tadDT$entrezID %in% rna_geneList,
                                     ]
        
        fc_cpt_gene2tad_dt <- gene2tadDT[gene2tadDT$entrezID %in% cpt_genes &
                                       gene2tadDT$entrezID %in% de_DT$genes2 ,
                                     ] # need to subset here for then next if keptTADs !
        
        # the DE table is filtered for minCount, not the corr data
        
        stopifnot(nrow(fc_cpt_gene2tad_dt) <= nrow(corr_cpt_gene2tad_dt))
        stopifnot(fc_cpt_gene2tad_dt$entrezID %in% corr_cpt_gene2tad_dt$entrezID)
        
        stopifnot(corr_cpt_gene2tad_dt$entrezID %in% rna_geneList) 
        stopifnot(fc_cpt_gene2tad_dt$entrezID %in% rna_geneList) 
        
        if(max(table(corr_cpt_gene2tad_dt$region)/nrow(corr_cpt_gene2tad_dt)) <= maxSameTAD){
          meanCorr <- paste0("sameTAD<=", maxSameTAD)
          corr_keptTADs <- corr_cpt_gene2tad_dt$region
          corr_keptGenes <- corr_cpt_gene2tad_dt$entrezID
          
        } else {
          rowsToKeep <- which(rna_geneList %in% cpt_genes)
          stopifnot(rowsToKeep == which(rna_geneList %in% corr_cpt_gene2tad_dt$entrezID))
          stopifnot(rownames(norm_rnaseqDT) == names(rna_geneList))

          
          if(length(rowsToKeep) < minGenes) {
            meanCorr <- paste0("<", minGenes, "genes")
            corr_keptTADs <- corr_cpt_gene2tad_dt$region
            corr_keptGenes <- corr_cpt_gene2tad_dt$entrezID
            
            
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
            corr_keptTADs <- corr_cpt_gene2tad_dt$region
            }

        }
        
        
        if(max(table(fc_cpt_gene2tad_dt$region)/nrow(fc_cpt_gene2tad_dt)) <= maxSameTAD){
          ratioDown <- paste0("sameTAD<=", maxSameTAD)
          fc_keptTADs <- fc_cpt_gene2tad_dt$region
          fc_keptGenes <- fc_cpt_gene2tad_dt$entrezID
        } else {
          stopifnot(fc_cpt_gene2tad_dt$entrezID %in% de_DT$genes2)
          cpt_de_DT <- de_DT[de_DT$genes2 %in% fc_cpt_gene2tad_dt$entrezID,]
          stopifnot(nrow(cpt_de_DT) == nrow(fc_cpt_gene2tad_dt))
          if(nrow(cpt_de_DT) < minGenes) {
            ratioDown <- paste0("<", minGenes, "genes")
            fc_keptTADs <- fc_cpt_gene2tad_dt$region
            fc_keptGenes <- fc_cpt_gene2tad_dt$entrezID
            
          } else {
            ratioDown <- sum(sign(cpt_de_DT$logFC) == -1)/nrow(cpt_de_DT)
            fc_keptGenes <- cpt_de_DT$genes2
            
            
            fc_keptTADs <- fc_cpt_gene2tad_dt$region
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
      
      names(all_meanCorr_ratioDown_famCpts) <- all_famCpts
      stopifnot(length(all_meanCorr_ratioDown_famCpts) == length(all_famCpts))
      
      outFile <- file.path(outFolder, hicds, exprds, "all_meanCorr_ratioDown_famCpts.Rdata")
      dir.create(dirname(outFile), recursive = TRUE)
      save(all_meanCorr_ratioDown_famCpts, file= outFile, version=2)
      cat(paste0("... written: ", outFile,  "\n"))
      
      # famCorr_data <- get(load("FAMILYMODULES_RUNMEANTADCORR/Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/all_meanCorr_ratioDown_famCpts.Rdata"))
      toKeep <- unlist(lapply(all_meanCorr_ratioDown_famCpts, function(x) is.numeric(x [["meanCorr"]])))
      stopifnot(length(toKeep) == length(all_meanCorr_ratioDown_famCpts))
      famCorr_dataF <- all_meanCorr_ratioDown_famCpts[toKeep]
      
      
      famCorr <- unlist(lapply(famCorr_dataF,  function(x) x[["meanCorr"]]))
      obsCorr <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "4_runMeanTADCorr", "all_meanCorr_TAD.Rdata" )))
      outFile <- file.path(outFolder, hicds, exprds, paste0(hicds, "_", exprds, "_obs_famCpt_meanCorr_density.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      plot_multiDens(
        list(famCmpnt_meanCorr = famCorr,
             obsTAD_meanCorr = obsCorr),
        plotTit = paste0(hicds, " -  ", exprds)
      )
      mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize, "; minGenes=", minGenes,  "; maxSameTAD=", maxSameTAD), font=3)
      foo <- dev.off()
      cat(paste0("... written: ", outFile,  "\n"))
      
      
      
      # famRatioDown_data <- get(load("FAMILYMODULES_RUNMEANTADCORR/Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/all_meanCorr_ratioDown_famCpts.Rdata"))
      
      toKeep <- unlist(lapply(all_meanCorr_ratioDown_famCpts, function(x) is.numeric(x [["ratioDown"]])))
      stopifnot(length(toKeep) == length(all_meanCorr_ratioDown_famCpts))
      famRatioDown_dataF <- all_meanCorr_ratioDown_famCpts[toKeep]
      famRatioDown <- unlist(lapply(famRatioDown_dataF,  function(x) x[["ratioDown"]]))
      obsRatioDown <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown", "all_obs_ratioDown.Rdata" )))
      
      
      
      outFile <- file.path(outFolder, hicds, exprds, paste0(hicds, "_", exprds, "_obs_famCpt_ratioDown_density.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      plot_multiDens(
        list(famCmpnt_ratioDown = famRatioDown,
             obsTAD_ratioDown = obsRatioDown),
        plotTit = paste0(hicds, " -  ", exprds)
      )
      mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize, "; minGenes=", minGenes,  "; maxSameTAD=", maxSameTAD), font=3)
      foo <- dev.off()
      cat(paste0("... written: ", outFile,  "\n"))
      
      list(famCmpnt_ratioDown = famRatioDown,
           obsTAD_ratioDown = obsRatioDown,
          famCmpnt_meanCorr = famCorr,
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


all_fam_corr <- lapply(all_meanCorr_ratioDown, function(sublist) lapply(sublist, function(x) x[["famCmpnt_meanCorr"]]))
all_obs_corr <- lapply(all_meanCorr_ratioDown, function(sublist) lapply(sublist, function(x) x[["obsTAD_meanCorr"]]))
nDS <- length(unlist(all_fam_corr, recursive = FALSE))


outFile <- file.path(outFolder, paste0("allDS_obs_famCpt_meanCorr_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(famCmpnt_meanCorr = unlist(all_fam_corr),
       obsTAD_meanCorr = unlist(all_obs_corr)),
  my_xlab = paste0("intra-TAD/component meanCorr"),
  plotTit = paste0( "famCpts - meanCorr - all datasets - n =", nDS )
)
mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize, "; minGenes=", minGenes,  "; maxSameTAD=", maxSameTAD), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))



all_fam_ratioDown <- lapply(all_meanCorr_ratioDown, function(sublist) lapply(sublist, function(x) x[["famCmpnt_ratioDown"]]))
all_obs_ratioDown <- lapply(all_meanCorr_ratioDown, function(sublist) lapply(sublist, function(x) x[["obsTAD_ratioDown"]]))
nDS <- length(unlist(all_fam_ratioDown, recursive = FALSE))


outFile <- file.path(outFolder, paste0("allDS_obs_famCpt_ratioDown_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(famCmpnt_ratioDown = unlist(all_fam_ratioDown),
       obsTAD_ratioDown = unlist(all_obs_ratioDown)),
  my_xlab = paste0("intra-TAD/component ratioDown"),
  plotTit = paste0( "famCpts - ratioDown - all datasets - n =", nDS )
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

