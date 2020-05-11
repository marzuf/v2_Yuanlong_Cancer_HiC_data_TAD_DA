#!/usr/bin/Rscript

startTime <- Sys.time()

# Rscript familyModules_runMeanTADCorr.R

script_name <- "familyModules_runMeanTADCorr.R"

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


outFolder <- file.path("FAMILYMODULES_RUNMEANTADCORR", nMaxSize)

inFolder <- file.path("PREP_FAMILYMODULES", nMaxSize)

all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

hicds = "Barutcu_MCF-10A_40kb"

all_hicds=all_hicds[1]

all_corr <- foreach(hicds = all_hicds) %do%{
  cat(paste0("... start: ", hicds, "\n"))
  exprds_corr <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
    cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
    
    # retrieve file 

    famMod_file <- file.path(inFolder, hicds, exprds, "all_fams_dt.Rdata")
    
    ## !! UPDATE HERE !!!
    
    stopifnot(file.exists(famMod_file))
    fam_dt <- get(load(famMod_file))
    fam_dt$entrezID <- as.character(fam_dt$entrezID)
    fam_dt$cpt <- as.character(fam_dt$cpt)
    
    # INPUT DATA
    gene2tadDT_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(gene2tadDT_file))
    gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
    gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
    
    gene2tadDT <- gene2tadDT[grepl("_TAD", gene2tadDT$region),]
    
    
    pipeline_geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
    rna_geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_geneList.Rdata")))
    
    
    sum(fam_dt$entrezID %in% names(pipeline_geneList)) # 2493
    sum(fam_dt$entrezID %in% pipeline_geneList) # 2495
    sum(fam_dt$entrezID %in% names(rna_geneList)) # 7045
    sum(fam_dt$entrezID %in% rna_geneList) # 7059
    
    norm_rnaseqDT <- eval(parse(text = load(file.path(pipFolder, hicds, exprds,
                                                      "0_prepGeneData", "rna_qqnorm_rnaseqDT.Rdata")))) 
    
    # stopifnot(rna_geneList %in% rownames(norm_rnaseqDT)) # ! wrong
    stopifnot(names(rna_geneList) %in% rownames(norm_rnaseqDT))
    # reorder
    norm_rnaseqDT <- norm_rnaseqDT[names(rna_geneList),]
    
    stopifnot(fam_dt$entrezID %in% gene2tadDT$entrezID)  ### I took only genes from TADs !!!!
    # stopifnot(fam_dt$entrezID %in% names(pipeline_geneList))  ### NOT TRUE !!! I took only genes from TADs !!!!
    
    all_famCpts <- unique(fam_dt$cpt)
    famCpt = all_famCpts[1]
    all_meanCorr_famCpts <- foreach(famCpt=all_famCpts) %dopar% {
      
      cpt_genes <- fam_dt$entrezID[as.character(fam_dt$cpt) == as.character(famCpt)]
      stopifnot(length(cpt_genes) >= minCmpntSize)
      
      cpt_gene2tad_dt <- gene2tadDT[gene2tadDT$entrezID %in% cpt_genes,]
      stopifnot(nrow(cpt_gene2tad_dt) == length(cpt_genes))
      
      if(max(table(cpt_gene2tad_dt$region)/nrow(cpt_gene2tad_dt)) > maxSameTAD) return(paste0("sameTAD>", maxSameTAD))
      
      rowsToKeep <- which(rna_geneList %in% cpt_genes)

      keptGenes <- rna_geneList[rowsToKeep]
      keptTADs <- cpt_gene2tad_dt$region
      
      if(length(rowsToKeep) < minGenes) return(paste0("<", minGenes, "genes"))
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
      list(
        meanCorr=meanCorr,
        keptGenes=keptGenes,
        keptTADs=keptTADs
      )
    }
    
    cat(paste0("... end intra-cpt correlation\n"))
    
    names(all_meanCorr_famCpts) <- all_famCpts
    stopifnot(length(all_meanCorr_famCpts) == length(all_famCpts))
    
    outFile <- file.path(outFolder, hicds, exprds, "all_meanCorr_famCpts.Rdata")
    dir.create(dirname(outFile), recursive = TRUE)
    save(all_meanCorr_famCpts, file= outFile)
    cat(paste0("... written: ", outFile,  "\n"))
    
    # famCorr_data <- get(load("FAMILYMODULES_RUNMEANTADCORR/Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/all_meanCorr_famCpts.Rdata"))
    famCorr_data <- all_meanCorr_famCpts
    famCorr_dataF <- famCorr_data[lengths(famCorr_data) == 3]
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
    
    list(famCmpnt_meanCorr = famCorr,
         obsTAD_meanCorr = obsCorr
         )
  }
  names(exprds_corr) <- all_exprds[[paste0(hicds)]]
  exprds_corr
}
names(all_corr) <- all_hicds

outFile <- file.path(outFolder,  "all_corr.Rdata")
dir.create(dirname(outFile), recursive = TRUE)
save(all_corr, file= outFile)
cat(paste0("... written: ", outFile,  "\n"))


cat(paste0("*** DONE: ", script_name, "\n"))




    
