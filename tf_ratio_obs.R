startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)
require(dplyr)

# Rscript tf_ratio_obs.R crisp
# Rscript tf_ratio_obs.R chea3_all
# Rscript tf_ratio_obs.R chea3_lung
# Rscript tf_ratio_obs.R trrust
# Rscript tf_ratio_obs.R motifmap
# 
#

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- ifelse(plotType=="png", 400, 7)
plotCex <- 1.4

hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAluad_norm_luad"
dsIn <- "chea3_lung"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1 | length(args) == 3)
dsIn <- args[1]
if(length(args) == 3) {
  all_hicds <- args[2]
  all_exprds <- args[3]
  
} else {
  all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))
}

stopifnot(dsIn %in% c("crisp",  "trrust", "motifmap", "chea3_all", "chea3_lung"))

outFolder <- file.path(paste0("TF_RATIO_OBS_", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)

stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)


if(buildData){
  
  
  if(dsIn == "crisp") {
    reg_file <- file.path("gene_set_library_crisp_processed.txt")
    reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
    cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
    reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    # regSymbol targetSymbol targetEntrezID
    # 2    ARID3A         A1BG              1
    # 3    ARID3A     A1BG-AS1         503538
    # 4    ARID3A         A1CF          29974
    
  } else if(dsIn == "chea3_all") {
    reg_file <- file.path("chea3_all_tissues_TFs_processed.txt")
    reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
    cat(paste0("with target Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
    reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
  } else if(dsIn == "chea3_lung") {
    reg_file <- file.path("chea3_lung_TFs_processed.txt")
    reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
    cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
    cat(paste0("with target Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
    reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
  } else if(dsIn == "trrust"){
    reg_file <- file.path("trrust_rawdata.human.tsv")
    reg_dt <- read.delim(reg_file, sep="\t", header=FALSE, stringsAsFactors = FALSE,
                         col.names = c("regSymbol", "targetSymbol", "direction", "ID"))
    cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
    cat(paste0("with target Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
    reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    
    
  }  else if(dsIn == "motifmap"){
    reg_file <- file.path("MOTIFMAP_ALLGENES/overlapDT_bp.Rdata")
    reg_dt <- get(load(reg_file))
    colnames(reg_dt)[colnames(reg_dt)=="entrezID"] <- "targetEntrezID"
    cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    
    reg_dt$regSymbol_init <- reg_dt$regSymbol
    
    reg_dt$regSymbol <- toupper(gsub(".+=(.+)", "\\1", reg_dt$regSymbol_init) )
    
  }
  
  stopifnot(!is.na(reg_dt$targetEntrezID))
  reg_dt <- reg_dt[reg_dt$regSymbol %in% names(symb2entrez),]
  cat(paste0("with Reg Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
  reg_dt$regEntrezID <- symb2entrez[reg_dt$regSymbol]
  reg_dt$regEntrezID <- as.character(reg_dt$regEntrezID)
  stopifnot(!is.na(reg_dt$regEntrezID))
  
  
  init_reg_dt <- reg_dt
  rm("reg_dt")

  
  
  all_ratio_data <- foreach(hicds = all_hicds) %dopar%{
    # all_ratio_data <- foreach(hicds = all_hicds[7]) %dopar%{
    
    reg_dt <- init_reg_dt
    
    cat(paste0("> START - ", hicds,"\n"))
    
    g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
    g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    stopifnot(!duplicated(g2t_dt$entrezID))
    g2t_vect <- setNames(g2t_dt$region, g2t_dt$entrezID)


    #here keep only the reg elements that are assigned to a TAD
    hicds_reg_dt <- left_join(reg_dt, g2t_dt[,c("entrezID", "region")], by=c("targetEntrezID" = "entrezID"), all.x=FALSE, all.y=FALSE)  
    colnames(hicds_reg_dt)[colnames(hicds_reg_dt) == "region"] <- "targetRegion"
    hicds_reg_dt <- left_join(hicds_reg_dt, g2t_dt[,c("entrezID", "region")], by=c("regEntrezID" = "entrezID"), all.x=FALSE, all.y=FALSE)  # force here to keep only regulatory elements that are in  TAD
    colnames(hicds_reg_dt)[colnames(hicds_reg_dt) == "region"] <- "regRegion"
    
    
    #here keep also target elements that are regulated by a reg element not necessarily in a TAD
    # (but the target is inside a TAD)
    hicds_target_withBD_dt <- left_join(reg_dt, g2t_dt[,c("entrezID", "region")], by=c("targetEntrezID" = "entrezID"), all.x=FALSE, all.y=FALSE)
    colnames(hicds_target_withBD_dt)[colnames(hicds_target_withBD_dt) == "region"] <- "targetRegion"
    

    hicds_targetOnly_dt <- hicds_reg_dt[,c("targetEntrezID", "targetRegion")]
    hicds_targetOnly_dt <- unique(hicds_targetOnly_dt)
    hicds_regOnly_dt <- hicds_reg_dt[,c("regEntrezID", "regRegion")]
    hicds_regOnly_dt <- unique(hicds_regOnly_dt)
    
    
    hicds_target_withBD_dt <- hicds_target_withBD_dt[,c("targetEntrezID", "targetRegion")]
    hicds_target_withBD_dt <- unique(hicds_target_withBD_dt)
    
    # # of regulated genes by TAD (regulated by a regulatory element that is in  TAD)
    hicds_nTarget_dt <- aggregate(targetEntrezID ~ targetRegion, data = hicds_targetOnly_dt, FUN=length)
    hicds_nTarget_dt_vect <- setNames(hicds_nTarget_dt$targetEntrezID, hicds_nTarget_dt$targetRegion)
    
    # # of regulatory elements by TAD
    hicds_nReg_dt <- aggregate(regEntrezID ~ regRegion, data = hicds_regOnly_dt, FUN=length)
    hicds_nReg_dt_vect <- setNames(hicds_nReg_dt$regEntrezID, hicds_nReg_dt$regRegion)
    
    
    # # of regulated genes by TAD (regulated by a regulatory element that is not necessarily in  TAD)
    hicds_nTargetWithBD_dt <- aggregate(targetEntrezID ~ targetRegion, data = hicds_target_withBD_dt, FUN=length)
    hicds_nTargetWithBD_dt_vect <- setNames(hicds_nTargetWithBD_dt$targetEntrezID, hicds_nTargetWithBD_dt$targetRegion)
    
    # # of genes by TAD
    hicds_nGenes_dt <- aggregate(entrezID ~ region, data = g2t_dt, FUN=length)
    hicds_nGenes_dt_vect <- setNames(hicds_nGenes_dt$entrezID, hicds_nGenes_dt$region)
    
    
    all_tads <- unique(g2t_dt$region)
    
    hicds_data <- list(
      hicds_nGenes = hicds_nGenes_dt_vect[paste0(all_tads)],
      hicds_nTarget = hicds_nTarget_dt_vect[paste0(all_tads)],
      hicds_nTargetWithBD = hicds_nTargetWithBD_dt_vect[paste0(all_tads)],
      hicds_nReg = hicds_nReg_dt_vect[paste0(all_tads)]
    )
    

    exprds_data <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      # exprds_data <- foreach(exprds = all_exprds[[paste0(hicds)]][1]) %do% {
      
      if(dsIn == "chea3_lung") {
        if(! (grepl("lusc", exprds) | grepl("luad", exprds))) return(NULL)
      }
      
      
      geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata") ))
      
      de_dt <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "1_runGeneDE", "DE_topTable.Rdata") ))
      
      stopifnot(geneList %in% g2t_dt$entrezID)
      exprds_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      exprds_all_tads <- unique(exprds_g2t_dt$region)
      
      
      # target gene is a pipeline gene - no constraint on regulatory element
      target_allReg_exprds_reg_dt <- reg_dt[reg_dt$targetEntrezID %in% geneList,]
      nrow(target_allReg_exprds_reg_dt)
      target_allReg_exprds_reg_dt <- left_join(target_allReg_exprds_reg_dt, g2t_dt[,c("entrezID", "region")], by=c("targetEntrezID"="entrezID"), all.x=FALSE, all.y=FALSE)
      colnames(target_allReg_exprds_reg_dt)[colnames(target_allReg_exprds_reg_dt) == "region"] <- "targetRegion"
      
      # target gene is a pipeline gene - regulatory element is in a TAD
      target_tadReg_exprds_reg_dt <- reg_dt[reg_dt$targetEntrezID %in% geneList & reg_dt$regEntrezID %in% exprds_g2t_dt$entrezID,]
      nrow(target_tadReg_exprds_reg_dt)
      target_tadReg_exprds_reg_dt <- left_join(target_tadReg_exprds_reg_dt, g2t_dt[,c("entrezID", "region")], by=c("targetEntrezID"="entrezID"), all.x=FALSE, all.y=FALSE)
      colnames(target_tadReg_exprds_reg_dt)[colnames(target_tadReg_exprds_reg_dt) == "region"] <- "targetRegion"
      target_tadReg_exprds_reg_dt <- left_join(target_tadReg_exprds_reg_dt, g2t_dt[,c("entrezID", "region")], by=c("regEntrezID"="entrezID"), all.x=FALSE, all.y=FALSE)
      colnames(target_tadReg_exprds_reg_dt)[colnames(target_tadReg_exprds_reg_dt) == "region"] <- "regRegion"
      
      
      # target gene is a pipeline gene - regulatory element is expressed
      target_exprReg_exprds_reg_dt <- reg_dt[reg_dt$targetEntrezID %in% geneList & reg_dt$regEntrezID %in% de_dt$genes,]
      nrow(target_exprReg_exprds_reg_dt)
      target_exprReg_exprds_reg_dt <- left_join(target_exprReg_exprds_reg_dt, g2t_dt[,c("entrezID", "region")], by=c("targetEntrezID"="entrezID"), all.x=FALSE, all.y=FALSE)
      colnames(target_exprReg_exprds_reg_dt)[colnames(target_exprReg_exprds_reg_dt) == "region"] <- "targetRegion"

      
      
      # target gene is a pipeline gene - regulatory element is in a TAD + expressed
      target_tadAndExprReg_exprds_reg_dt <- reg_dt[reg_dt$targetEntrezID %in% geneList  & reg_dt$regEntrezID %in% exprds_g2t_dt$entrezID & reg_dt$regEntrezID %in% de_dt$genes,]
      nrow(target_tadAndExprReg_exprds_reg_dt)
      target_tadAndExprReg_exprds_reg_dt <- left_join(target_tadAndExprReg_exprds_reg_dt, g2t_dt[,c("entrezID", "region")], by=c("targetEntrezID"="entrezID"), all.x=FALSE, all.y=FALSE)
      colnames(target_tadAndExprReg_exprds_reg_dt)[colnames(target_tadAndExprReg_exprds_reg_dt) == "region"] <- "targetRegion"
      target_tadAndExprReg_exprds_reg_dt <- left_join(target_tadAndExprReg_exprds_reg_dt, g2t_dt[,c("entrezID", "region")], by=c("regEntrezID"="entrezID"), all.x=FALSE, all.y=FALSE)
      colnames(target_tadAndExprReg_exprds_reg_dt)[colnames(target_tadAndExprReg_exprds_reg_dt) == "region"] <- "regRegion"
    
      
      
      exprds_nGenes_dt <- aggregate(entrezID ~region , data=exprds_g2t_dt, FUN=length)
      exprds_nGenes_dt_vect <- setNames(exprds_nGenes_dt$entrezID, exprds_nGenes_dt$region)
      
      
      # target gene is a pipeline gene - no constraint on regulatory element
      countTarget_target_allReg_exprds_reg_dt <- aggregate(targetEntrezID ~ targetRegion, data=unique(target_allReg_exprds_reg_dt[,c("targetEntrezID", "targetRegion")]), FUN=length)
      countTarget_target_allReg_exprds_reg_vect <- setNames(countTarget_target_allReg_exprds_reg_dt$targetEntrezID, countTarget_target_allReg_exprds_reg_dt$targetRegion)
      
      # target gene is a pipeline gene - regulatory element is in a TAD
      countTarget_target_tadReg_exprds_reg_dt <- aggregate(targetEntrezID ~ targetRegion, data=unique(target_tadReg_exprds_reg_dt[,c("targetEntrezID", "targetRegion")]), FUN=length)
      countTarget_target_tadReg_exprds_reg_vect <- setNames(countTarget_target_tadReg_exprds_reg_dt$targetEntrezID, countTarget_target_tadReg_exprds_reg_dt$targetRegion)
      
      countReg_target_tadReg_exprds_reg_dt <- aggregate(regEntrezID ~ regRegion, data=unique(target_tadReg_exprds_reg_dt[,c("regEntrezID", "regRegion")]), FUN=length)
      countReg_target_tadReg_exprds_reg_vect <- setNames(countReg_target_tadReg_exprds_reg_dt$regEntrezID, countReg_target_tadReg_exprds_reg_dt$regRegion)
      
      # target gene is a pipeline gene - regulatory element is expressed
      countTarget_target_exprReg_exprds_reg_dt <- aggregate(targetEntrezID ~ targetRegion, data=unique(target_exprReg_exprds_reg_dt[,c("targetEntrezID", "targetRegion")]), FUN=length)
      countTarget_target_exprReg_exprds_reg_vect <- setNames(countTarget_target_exprReg_exprds_reg_dt$targetEntrezID, countTarget_target_exprReg_exprds_reg_dt$targetRegion)
      
      # target gene is a pipeline gene - regulatory element is in a TAD + expressed
      countTarget_target_tadAndExprReg_exprds_reg_dt <- aggregate(targetEntrezID ~ targetRegion, data=unique(target_tadAndExprReg_exprds_reg_dt[,c("targetEntrezID", "targetRegion")]), FUN=length)
      countTarget_target_tadAndExprReg_exprds_reg_vect <- setNames(countTarget_target_tadAndExprReg_exprds_reg_dt$targetEntrezID, countTarget_target_tadAndExprReg_exprds_reg_dt$targetRegion)
      
      countReg_target_tadAndExprReg_exprds_reg_dt <- aggregate(regEntrezID ~ regRegion, data=unique(target_tadAndExprReg_exprds_reg_dt[,c("regEntrezID", "regRegion")]), FUN=length)
      countReg_target_tadAndExprReg_exprds_reg_vect <- setNames(countReg_target_tadAndExprReg_exprds_reg_dt$regEntrezID, countReg_target_tadAndExprReg_exprds_reg_dt$regRegion)
      
      list(
        exprds_nGenes = exprds_nGenes_dt_vect[paste0(exprds_all_tads)],
        # target gene is a pipeline gene - no constraint on regulatory element
        exprds_nTargetWithBD =countTarget_target_allReg_exprds_reg_vect[paste0(exprds_all_tads)],
        # target gene is a pipeline gene - regulatory element is in a TAD -> # of target
        exprds_nTarget = countTarget_target_tadReg_exprds_reg_vect[paste0(exprds_all_tads)],
        # target gene is a pipeline gene - regulatory element is in a TAD -> # of reg element
        exprds_nReg = countReg_target_tadReg_exprds_reg_vect[paste0(exprds_all_tads)],
        # target gene is a pipeline gene - regulatory element is expressed
        exprds_nTargetWithBDExpr = countTarget_target_exprReg_exprds_reg_vect[paste0(exprds_all_tads)],
        # target gene is a pipeline gene - regulatory element is in TAD +  expressed -> # of target
        exprds_nTargetExpr = countTarget_target_tadAndExprReg_exprds_reg_vect[paste0(exprds_all_tads)],
        # target gene is a pipeline gene - regulatory element is in TAD +  expressed -> # of reg element
        exprds_nRegExpr = countReg_target_tadAndExprReg_exprds_reg_vect[paste0(exprds_all_tads)]
      )

    } # end-foreach iterating exprds
    
    if(!is.null(exprds_data))
      names(exprds_data) <- all_exprds[[paste0(hicds)]]
    
    list(
      hicds_data, 
      exprds_data
    )
  } # end-for iterating over hicds
  names(all_ratio_data) <- all_hicds
  outFile <- file.path(outFolder, paste0("all_ratio_data.Rdata"))
  save(all_ratio_data, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile="TF_RATIO_OBS_CHEA3_LUNG/all_ratio_data.Rdata"
  outFile <- file.path(outFolder, paste0("all_ratio_data.Rdata"))
  all_ratio_data <- get(load(outFile))
}  

all_hicds_nTargetByGenes <- unlist(lapply(all_ratio_data, function(sl) sl[[1]][[paste0("hicds_nTarget")]]/sl[[1]][[paste0("hicds_nGenes")]]))
all_hicds_nTargetWithBDByGenes <- unlist(lapply(all_ratio_data, function(sl) sl[[1]][[paste0("hicds_nTargetWithBD")]]/sl[[1]][[paste0("hicds_nGenes")]]))
all_hicds_nRegByGenes <- unlist(lapply(all_ratio_data, function(sl) sl[[1]][[paste0("hicds_nReg")]]/sl[[1]][[paste0("hicds_nGenes")]]))

stopifnot(length(all_hicds_nTargetByGenes) == length(all_hicds_nTargetWithBDByGenes))
stopifnot(length(all_hicds_nTargetByGenes) == length(all_hicds_nRegByGenes))


all_exprds_nGenes <- unlist(lapply(all_ratio_data, function(sl) lapply(sl[[2]], function(x) x[[paste0("exprds_nGenes")]])))






#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))
