#!/usr/bin/Rscript

startTime <- Sys.time()

# Rscript familyModules_runTADratioDown.R

script_name <- "familyModules_runTADratioDown.R"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


plotType <- "svg"
myHeight <- 5
myWidth <- 7

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


nMaxSize <- 1


outFolder <- file.path("FAMILYMODULES_RUNTADRATIODOWN", nMaxSize)
inFolder <- file.path("PREP_FAMILYMODULES", nMaxSize)


all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
# all_hicds=all_hicds[1]
# all_hicds=all_hicds[2:length(all_hicds)]
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

hicds = "Barutcu_MCF-10A_40kb"

all_hicds=all_hicds[1:2]
exprds="TCGAbrca_lum_bas"

buildData <- TRUE

if(buildData){
  
  all_ratioDown <- foreach(hicds = all_hicds) %do%{
    cat(paste0("... start: ", hicds, "\n"))
    exprds_ratioDown <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
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
      
      gene2tadDT <- gene2tadDT[grepl("_TAD", gene2tadDT$region),]
      
      
      pipeline_geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
      rna_geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_geneList.Rdata")))
      
      de_DT <-  get(load(file.path(pipFolder, hicds, exprds, "1_runGeneDE", "DE_topTable.Rdata")))
      # stopifnot(names(rna_geneList) %in% de_DT$genes) FALSE
      # stopifnot(de_DT$genes %in% rna_geneList ) # FALSE
      stopifnot(de_DT$genes %in% names(rna_geneList) )
      
      sum(fam_dt$entrezID %in% names(pipeline_geneList)) # 2493
      sum(fam_dt$entrezID %in% pipeline_geneList) # 2495
      sum(fam_dt$entrezID %in% names(rna_geneList)) # 7045
      sum(fam_dt$entrezID %in% rna_geneList) # 7059
      
      sum(names(rna_geneList) %in% de_DT$genes)
      sum((rna_geneList) %in% de_DT$genes)
      
      de_DT <- de_DT[de_DT$genes %in% names(rna_geneList),]
      nrow(de_DT)
      rna_geneList <- rna_geneList[names(rna_geneList) %in% de_DT$genes]
      
      # stopifnot(rna_geneList %in% rownames(norm_rnaseqDT)) # ! wrong
      # stopifnot(names(rna_geneList) %in% rownames(norm_rnaseqDT))
      # reorder
      # norm_rnaseqDT <- norm_rnaseqDT[names(rna_geneList),]
      
      stopifnot(fam_dt$entrezID %in% gene2tadDT$entrezID)  ### I took only genes from TADs !!!!
      # stopifnot(fam_dt$entrezID %in% names(pipeline_geneList))  ### NOT TRUE !!! I took only genes from TADs !!!!
      
      all_famCpts <- unique(fam_dt$cpt)
      famCpt = all_famCpts[1]
      all_ratioDown_famCpts <- foreach(famCpt=all_famCpts) %dopar% {
        
        cpt_genes <- fam_dt$entrezID[as.character(fam_dt$cpt) == as.character(famCpt)]
        stopifnot(length(cpt_genes) >= minCmpntSize)
        
        cpt_gene2tad_dt <- gene2tadDT[gene2tadDT$entrezID %in% cpt_genes,]
        stopifnot(nrow(cpt_gene2tad_dt) == length(cpt_genes))
        
        keptTADs <- cpt_gene2tad_dt$region
        
        if(max(table(cpt_gene2tad_dt$region)/nrow(cpt_gene2tad_dt)) > maxSameTAD) return(paste0("sameTAD>", maxSameTAD))
        
        stopifnot(de_DT$genes %in% names(rna_geneList) )
        
        cpt_rna_geneList <- rna_geneList[rna_geneList %in% cpt_genes]
        cpt_de_DT <- de_DT[de_DT$genes %in% names(cpt_rna_geneList),]
        
        stopifnot(cpt_de_DT$genes %in% names(cpt_rna_geneList) )
        stopifnot(length(cpt_rna_geneList) == nrow(cpt_de_DT))
        
        if(nrow(cpt_de_DT) < minGenes) return(paste0("<", minGenes, "genes"))
        
        cpt_ratioDown <- sum(sign(cpt_de_DT$logFC) == -1)/nrow(cpt_de_DT)
        
        
        list(
          ratioDown=cpt_ratioDown,
          keptGenes=cpt_rna_geneList,
          keptTADs=keptTADs
        )
      }
      
      cat(paste0("... end intra-cpt ratioDown\n"))
      
      names(all_ratioDown_famCpts) <- all_famCpts
      stopifnot(length(all_ratioDown_famCpts) == length(all_famCpts))
      
      outFile <- file.path(outFolder, hicds, exprds, "all_ratioDown_famCpts.Rdata")
      dir.create(dirname(outFile), recursive = TRUE)
      save(all_ratioDown_famCpts, file= outFile)
      cat(paste0("... written: ", outFile,  "\n"))
      
      # famRatioDown_data <- get(load("FAMILYMODULES_RUNMEANTADCORR/Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/all_ratioDown_famCpts.Rdata"))
      famRatioDown_data <- all_ratioDown_famCpts
      famRatioDown_dataF <- famRatioDown_data[lengths(famRatioDown_data) == 3]
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
           obsTAD_ratioDown = obsRatioDown
      )
    }
    names(exprds_ratioDown) <- all_exprds[[paste0(hicds)]]
    exprds_ratioDown
  }
  names(all_ratioDown) <- all_hicds
  
  outFile <- file.path(outFolder,  "all_ratioDown.Rdata")
  dir.create(dirname(outFile), recursive = TRUE)
  save(all_ratioDown, file= outFile, version=2)
  cat(paste0("... written: ", outFile,  "\n"))
  

  
  
} else {
  
  outFile <- file.path(outFolder,  "all_ratioDown.Rdata")
  all_ratioDown <- get(load(outFile))
  
}

all_fam_ratioDown <- lapply(all_ratioDown, function(sublist) lapply(sublist, function(x) x[["famCmpnt_ratioDown"]]))
all_obs_ratioDown <- lapply(all_ratioDown, function(sublist) lapply(sublist, function(x) x[["obsTAD_ratioDown"]]))
nDS <- length(unlist(all_fam_ratioDown, recursive = FALSE))


outFile <- file.path(outFolder, paste0("allDS_obs_famCpt_ratioDown_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  list(famCmpnt_ratioDown = unlist(all_fam_ratioDown),
       obsTAD_ratioDown = unlist(all_obs_ratioDown)),
  my_xlab = paste0("intra-TAD/component ratioDown"),
  plotTit = paste0( "all datasets - n =", nDS )
)
mtext(side=3, text = paste0("minCmpntSize=", minCmpntSize, "; minGenes=", minGenes,  "; maxSameTAD=", maxSameTAD), font=3)
foo <- dev.off()
cat(paste0("... written: ", outFile,  "\n"))



cat(paste0("*** DONE: ", script_name, "\n"))

# 
# 
# 
# #!/usr/bin/Rscript
# 
# startTime <- Sys.time()
# 
# ################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# # - script0: pipeline_regionList.Rdata
# # - script0: rna_geneList.Rdata
# # - script0: pipeline_geneList.Rdata
# # - script0: rna_madnorm_rnaseqDT.Rdata
# # - script1: DE_topTable.Rdata
# # - script1: DE_geneList.Rdata
# ################################################################################
# 
# ################  OUTPUT
# # - /all_meanLogFC_TAD.Rdata
# ################################################################################
# 
# SSHFS <- F
# setDir <- ifelse(SSHFS, "/media/electron", "")
# 
# args <- commandArgs(trailingOnly = TRUE)
# stopifnot(length(args) == 1)
# settingF <- args[1]
# stopifnot(file.exists(settingF))
# 
# pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
# 
# script0_name <- "0_prepGeneData"
# script1_name <- "1_runGeneDE"
# script_name <- "3_runMeanTADLogFC"
# stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
# cat(paste0("> START ", script_name,  "\n"))
# 
# source("main_settings.R")
# source(settingF)
# source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
# suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
# suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
# suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
# 
# # create the directories
# curr_outFold <- paste0(pipOutFold, "/", script_name)
# system(paste0("mkdir -p ", curr_outFold))
# 
# pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
# system(paste0("rm -f ", pipLogFile))
# 
# registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R
# 
# # ADDED 16.11.2018 to check using other files
# txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
# printAndLog(txt, pipLogFile)
# txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
# printAndLog(txt, pipLogFile)
# txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
# printAndLog(txt, pipLogFile)
# txt <- paste0("settingF\t=\t", settingF, "\n")
# printAndLog(txt, pipLogFile)
# 
# ################################***********************************************************************************
# ############ LOAD INPUT DATA
# ################################***********************************************************************************
# gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
# gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
# 
# DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
# DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
# 
# pipeline_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
# pipeline_regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))
# 
# if(useTADonly) {
#   if(any(grepl("_BOUND", pipeline_regionList))) {
#     stop("! data were not prepared for \"useTADonly\" !")
#   }
# }
# 
# stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
# stopifnot(!any(duplicated(names(DE_geneList))))
# 
# entrezList <- unlist(sapply(DE_topTable$genes, function(x) DE_geneList[x]))
# names(entrezList) <- DE_topTable$genes
# stopifnot(length(entrezList) == length(DE_topTable$genes))
# 
# # replace the gene symbol rownames by ensemblID rownames
# logFC_DT <- data.frame(entrezID =  entrezList,
#                        logFC = DE_topTable$logFC, stringsAsFactors = F)
# 
# rownames(logFC_DT) <- NULL
# initNrow <- nrow(logFC_DT)
# logFC_DT <- logFC_DT[logFC_DT$entrezID %in% pipeline_geneList,]
# txt <- paste0(toupper(script_name), "> Take only filtered genes: ", nrow(logFC_DT), "/", initNrow, "\n")
# printAndLog(txt, pipLogFile)
# 
# ### take only the filtered data according to initial settings
# gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(pipeline_geneList),]
# initLen <- length(unique(gene2tadDT$region))
# gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]
# txt <- paste0(toupper(script_name), "> Take only filtered regions: ", length(unique(gene2tadDT$region)), "/", initLen, "\n")
# printAndLog(txt, pipLogFile)
# 
# ################################***********************************************************************************
# ################################********************************************* get observed logFC for all regions
# ################################***********************************************************************************
# 
# cat(paste0("... start computing mean logFC by TAD \n"))
# 
# head(logFC_DT)
# 
# mergedDT <- left_join(logFC_DT, gene2tadDT[,c("entrezID", "region")], by="entrezID")
# 
# 
# save(mergedDT, file="mergedDT.Rdata")
# save(logFC_DT, file="logFC_DT.Rdata")
# save(gene2tadDT, file="gene2tadDT.Rdata")
# 
# stopifnot(nrow(mergedDT) == nrow(na.omit(mergedDT)))
# 
# mean_DT <- aggregate(logFC ~ region, data=mergedDT, FUN=mean)
# all_meanLogFC_TAD <- setNames(mean_DT$logFC, mean_DT$region)
# stopifnot(length(all_meanLogFC_TAD) == length(unique(gene2tadDT$region)))
# txt <- paste0(toupper(script_name), "> Number of regions for which mean logFC computed: ", length(all_meanLogFC_TAD), "\n")
# printAndLog(txt, pipLogFile)
# 
# if(useTADonly) {
#     initLen <- length(all_meanLogFC_TAD)
#     all_meanLogFC_TAD <- all_meanLogFC_TAD[grep("_TAD", names(all_meanLogFC_TAD))]
#     txt <- paste0(toupper(script_name), "> Take only the TAD regions: ", length(all_meanLogFC_TAD),"/", initLen, "\n")
#     printAndLog(txt, pipLogFile)
# }
# 
# save(all_meanLogFC_TAD, file= paste0(curr_outFold, "/all_meanLogFC_TAD.Rdata"))
# cat(paste0("... written: ", curr_outFold, "/all_meanLogFC_TAD.Rdata", "\n"))
# 
# txt <- paste0(startTime, "\n", Sys.time(), "\n")
# printAndLog(txt, pipLogFile)
# 
# cat(paste0("*** DONE: ", script_name, "\n"))
# 
