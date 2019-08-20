suppressPackageStartupMessages(library(grid, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(tools, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


setDir = ifelse(SSHFS, "/media/electron", "")


mainDir = file.path(setDir, "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA")


#source(file.path(mainDir, "utils_fct.R"))
source(file.path("..", "Cancer_HiC_data_TAD_DA", "utils_fct.R"))


pipelineDir = file.path(mainDir, "PIPELINE")
settingDir = file.path(pipelineDir, "INPUT_FILES")
pipOutDir =  file.path(pipelineDir, "OUTPUT_FOLDER")
pipScriptDir = paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
mainSettingsFile = file.path(paste0(pipScriptDir, "_",  "TopDom"), "main_settings.R")

plot_lolliTAD_ds <- function(exprds, hicds, all_TADs, orderByLolli = "startPos", titTADonly=TRUE){
                          # mainDir = file.path(setDir, "/mnt/etemp/marie/Cancer_HiC_data_TAD_DA"),
                          # pipelineDir = file.path(mainDir, "PIPELINE"),
                          # settingDir = file.path(pipelineDir, "INPUT_FILES"),
                          # #pipOutDir =  file.path(pipelineDir, "OUTPUT_FOLDER"),
                          # pipOutDir = file.path(pipOutDir, hicds, exprds),
                          # pipScriptDir = paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2"),
                          # mainSettingsFile = file.path(paste0(pipScriptDir, "_",  "TopDom"), "main_settings.R") # needed for entrezDT_file
                          # ) {
  
  stopifnot(length(exprds) == 1)
  stopifnot(length(hicds) == 1)

  pipOutDir = file.path(pipOutDir, hicds, exprds)
  
  stopifnot(dir.exists(file.path(mainDir, hicds)))
  stopifnot(dir.exists(settingDir))
  stopifnot(dir.exists(pipOutDir))
  stopifnot(dir.exists(pipScriptDir))
  stopifnot(file.exists(mainSettingsFile))
  
  # PIPELINE/INPUT_FILES/NCI-H460_40kb/run_settings_TCGAluad_mutKRAS_mutEGFR.R
  settingF <- file.path(settingDir, hicds, paste0("run_settings_", exprds, ".R"))
  stopifnot(file.exists(settingF))
  
  script0_name <- "0_prepGeneData"
  script1_name <- "1_runGeneDE"
  script3_name <- "3_runMeanTADLogFC"
  
  cat("... source files\n")
  
  source(mainSettingsFile)
  source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
  source(settingF) # overwrite the variables loaded from mainSettingsFile
  
  ################################****************************************************************************************
  ####################################################### PREPARE INPUT
  ################################****************************************************************************************
  
  # INPUT DATA
  # txt <- paste0("... gene2tadDT_file = ", gene2tadDT_file, "\n")
  # printAndLog(txt, logFile = logfile)
  
  cat(paste0("... read gene2tadDT from: ", gene2tadDT_file, "\n"))
  
  gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
  gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
  
  cat(paste0("... read entrez2symbDT from: ", entrezDT_file, "\n"))
  
  # txt <- paste0("... entrezDT_file = ", entrezDT_file, "\n")
  # printAndLog(txt, logFile = logfile)
  entrez2symbDT <- read.delim(entrezDT_file, header=T, stringsAsFactors=F)
  entrez2symbDT <- entrez2symbDT[,c("entrezID", "symbol")]
  colnames(entrez2symbDT) <- c("entrezID", "geneName")
  entrez2symbDT$entrezID <- as.character(entrez2symbDT$entrezID)
  
  meanTADlogFC <- eval(parse(text = load(file.path(pipOutDir,script3_name,"all_meanLogFC_TAD.Rdata"))))
  
  ########################################################################################
  ########### PREPARE THE GENES TO PLOT THE logFC 
  ########################################################################################
  samp1File <- file.path(setDir, sample1_file)
  cat(paste0("... load ", samp1File, "\n"))
  samp1 <- eval(parse(text = load(samp1File)))
  samp2File <- file.path(setDir, sample2_file)
  cat(paste0("... load ", samp2File, "\n"))
  samp2 <- eval(parse(text = load(samp2File)))
  
  # UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
  rnaseqDTFile <- file.path(pipOutDir,script0_name, "rna_rnaseqDT.Rdata")
  initListFile <- file.path(pipOutDir,script0_name, "rna_geneList.Rdata")
  geneListFile <- file.path(pipOutDir,script0_name, "pipeline_geneList.Rdata")

  cat(paste0("... load ", rnaseqDTFile, "\n"))
  rnaseqDT <- eval(parse(text = load(rnaseqDTFile)))
  cat(paste0("... load ", initListFile, "\n"))
  initList <- eval(parse(text = load(initListFile)))
  cat(paste0("... load ", geneListFile, "\n"))
  geneList <- eval(parse(text = load(geneListFile)))
  
  # txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
  # printAndLog(txt, logfile)
  
  rnaseqDT <- rnaseqDT[names(geneList),]    
  stopifnot(all(rownames(rnaseqDT) == names(geneList)))
  stopifnot(!any(duplicated(names(geneList))))
  stopifnot(!any(duplicated(geneList)))
  
  DE_topTable <- eval(parse(text = load(file.path(pipOutDir, script1_name, "DE_topTable.Rdata"))))
  stopifnot(!any(duplicated(names(geneList))))
  DE_topTable <- DE_topTable[DE_topTable$genes %in% names(geneList),]
  stopifnot(nrow(DE_topTable) > 0)
  
  gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
  
  DE_topTable$genes <- unlist(sapply(DE_topTable$genes, function(x) geneList[x]))
  rownames(DE_topTable) <- NULL
  
  # ! duplicated row.names are not allowed !
  dupEntrez <- unique(geneList[duplicated(geneList)])
  geneList <- geneList[! geneList %in% dupEntrez]
  stopifnot(!any(duplicated(geneList)))
  
  DE_topTable <- DE_topTable[!DE_topTable$genes %in% dupEntrez,]
  stopifnot(!any(duplicated(DE_topTable$genes)))
  
  initNrow <- nrow(rnaseqDT)
  rnaseqDT <- rnaseqDT[which(rownames(rnaseqDT) %in% names(geneList)),]
  # txt <- paste0(toupper(script_name), "> Discard duplicated symbol, retain: ", nrow(rnaseqDT), "/", initNrow , " genes\n")
  # printAndLog(txt, logfile)
  
  stopifnot(all(rownames(rnaseqDT) == names(geneList)))
  stopifnot(is.numeric(rnaseqDT[1,1]))
  #if(applyVoomAndCPM) {
  if(inputDataType == "raw" | inputDataType == "RSEM") {
    log2_rnaseqDT <- log2(rnaseqDT + 0.0001) 
  } else{
    log2_rnaseqDT <- rnaseqDT
  }
  meanExpr <- rowMeans(log2_rnaseqDT, na.rm=T)
  names(meanExpr) <- unlist(sapply(names(meanExpr), function(x) geneList[x]))
  
  ################################****************************************************************************************
  ####################################################### DO THE PLOTS 
  ################################****************************************************************************************
  ###################################################################### PLOT LOLLI TAD
  
  # retrieve which condition is cond1, i.e. the one that is more expressed when logFC is positive
  geneHighestLogFC <- names(geneList[geneList == DE_topTable$genes[which.max(DE_topTable$logFC)] ])
  
  samp1_vect <- rnaseqDT[geneHighestLogFC,samp1, drop=F]
  stopifnot(dim(samp1_vect) == c(1,length(samp1)))
  samp2_vect <- rnaseqDT[geneHighestLogFC,samp2, drop=F]
  stopifnot(dim(samp2_vect) == c(1,length(samp2)))
  
  plot_cond1 <- ifelse(as.numeric(rowMeans(samp1_vect, na.rm = T)) > as.numeric(rowMeans(samp2_vect, na.rm = T)), cond1, cond2)
  plot_cond2 <- ifelse(as.numeric(rowMeans(samp1_vect, na.rm = T)) > as.numeric(rowMeans(samp2_vect, na.rm = T)), cond2, cond1)
  
  n_topTAD_toplot <- length(all_TADs)
  stopifnot(length(n_topTAD_toplot) > 0)
  
  cat("... start plotting\n")
  
  vect_plot <- list()
  for(i_plot in 1:n_topTAD_toplot)  {
    tad_to_plot <- all_TADs[i_plot]
    vect_plot[[i_plot]] <- plot_lolliTAD(TAD_to_plot = tad_to_plot,
                                         meanExpr_vect = meanExpr, 
                                         DE_table = DE_topTable,
                                         g2t_table = gene2tadDT, 
                                         id2name_table=entrez2symbDT, 
                                         geneList = geneList,
                                         textLeft =  meanTADlogFC[tad_to_plot] > 0,
                                         orderBy = orderByLolli,
                                         cond1=plot_cond1, cond2=plot_cond2)
  }
  if(titTADonly) {
    myTit <- ""
  } else {
    myTit <-   paste(hicds," - " , exprds)
  }
  all_plots <- do.call(grid.arrange, c(vect_plot,  list(ncol=ifelse(n_topTAD_toplot == 1, 1, 2), top=textGrob(myTit,
                                                                             gp=gpar(fontsize=20,font=2)))))
  
  return(all_plots)
  
  
}





