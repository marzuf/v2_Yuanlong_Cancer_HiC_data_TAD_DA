
# Rscript check_diffTADs_randommidpos.R

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "check_diffTADs_randommidpos.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"

require(foreach)
require(doMC)
registerDoMC(40)

pipFolder<- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CHECK_DIFFTADS_RANDOMMIDPOS" 
dir.create(outFolder, recursive = TRUE)


all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[grepl("_RANDOMMIDPOS", all_hicds) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

  
hicds = all_hicds[1]
all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    cat("... start ", hicds, " - ", exprds, "\n")
    
    ### RETRIEVE THE GENE2TAD ASSIGNMENT
    g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    
    ### RETRIEVE REGIONLIST
    stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
    regionListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
    stopifnot(file.exists(regionListFile))
    pipeline_regionList <- eval(parse(text = load(regionListFile))) # not adjusted
    
    ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
    stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
    geneListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
    stopifnot(pipeline_geneList %in% g2t_DT$entrezID)    
    g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList,]
    
    ### RETRIEVE THE GENE2TAD ASSIGNMENT - INIT
    init_hicds <- gsub("_RANDOMMIDPOS_40kb", "_40kb", hicds)
    init_g2tFile <- file.path(pipFolder, init_hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(init_g2tFile))
    init_g2t_DT <- read.delim(init_g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "init_region"), stringsAsFactors = FALSE)
    init_g2t_DT$entrezID <- as.character(init_g2t_DT$entrezID)
    
    cmp_init_g2t_dt <- merge(g2t_DT[,c("entrezID", "region")], init_g2t_DT[,c("entrezID", "init_region")], all.x = TRUE, all.y=FALSE, by="entrezID")
    stopifnot(!is.na(cmp_init_g2t_dt))
    stopifnot(nrow(cmp_init_g2t_dt) == nrow(g2t_DT))
    
    cmp_init_g2t_dt_noBD <- cmp_init_g2t_dt[!grepl("_BOUND", cmp_init_g2t_dt$init_region),]
    initRegByReg_dt_noBD <- aggregate(init_region ~ region, data =cmp_init_g2t_dt_noBD, FUN = function(x) length(unique(x)))
    
    initRegByReg_dt <- aggregate(init_region ~ region, data =cmp_init_g2t_dt, FUN = function(x) length(unique(x)))
    mean(initRegByReg_dt$init_region == 1)
    mean(initRegByReg_dt$init_region == 2)
    mean(initRegByReg_dt$init_region == 3)
    
    # RETAIN ONLY THE GENES THAT ENCOMPASS GENES FROM DIFFERENT TADs:
    keepTADs <- initRegByReg_dt$region[initRegByReg_dt$init_region > 1]
    stopifnot(length(keepTADs) == sum(initRegByReg_dt$init_region == 2) + sum(initRegByReg_dt$init_region == 3))
    stopifnot(keepTADs %in% g2t_DT$region)
    stopifnot(!duplicated(keepTADs))
    stopifnot(initRegByReg_dt$init_region == 1 | initRegByReg_dt$init_region == 2  | initRegByReg_dt$init_region == 3)  # might be 3 because of the BD
    stopifnot(initRegByReg_dt_noBD$init_region == 1 | initRegByReg_dt_noBD$init_region == 2  | initRegByReg_dt_noBD$init_region == 3)  # might be 3 because of rounding BD pos
    
    keep_g2t_dt <- g2t_DT[g2t_DT$region %in% keepTADs,]
    stopifnot(length(unique(keep_g2t_dt$region)) == length(keepTADs))
    
    stopifnot(keep_g2t_dt$entrezID %in% pipeline_geneList)
    stopifnot(keep_g2t_dt$region %in% pipeline_regionList)
    
    new_pipeline_geneList <- pipeline_geneList[ pipeline_geneList %in% keep_g2t_dt$entrezID]
    new_pipeline_regionList <- pipeline_regionList[ pipeline_regionList %in% keep_g2t_dt$region]
    
    stopifnot(length(new_pipeline_regionList) > 0)
    stopifnot(length(new_pipeline_geneList) > 0)
    
    outFile <- file.path(outFolder,hicds,exprds,"pipeline_geneList.Rdata")
    dir.create(dirname(outFile), recursive = TRUE)
    save(new_pipeline_geneList, file = outFile, version=2)
    cat(paste0("... written: ", outFile, "\n"))
    
    outFile <- file.path(outFolder,hicds,exprds,"pipeline_regionList.Rdata")
    dir.create(dirname(outFile), recursive = TRUE)
    save(new_pipeline_regionList, file = outFile, version=2)
    cat(paste0("... written: ", outFile, "\n"))
    
    data.frame(
      hicds =hicds,
      exprds=exprds,
      ratio_onlyOneTAD = mean(initRegByReg_dt$init_region == 1),
      nKeepTADs = length(keepTADs),
      stringsAsFactors = FALSE
    )
    
    
  }
}

outFile <- file.path(outFolder, "all_result_dt.txt")
write.table(all_result_dt, file = outFile, sep="\t", quote=F, append=F, col.names=T, row.names=F)
cat(paste0("... written: ", outFile, "\n"))


# 
# ### RETRIEVE THE TAD ASSIGNMENT
# tadFile <- file.path(pipFolder, hicds, "genes2tad", "all_assigned_regions.txt")
# stopifnot(file.exists(tadFile))
# tad_DT <- read.delim(tadFile, header=F, col.names = c("chromo","region", "start", "end"), stringsAsFactors = FALSE)
# tad_DT <- tad_DT[grepl("_TAD", tad_DT$region),]
# stopifnot(nrow(tad_DT) > 0)
# 
# 
# init_tadFile <- file.path(pipFolder, init_hicds, "genes2tad", "all_assigned_regions.txt")
# stopifnot(file.exists(tadFile))
# init_tad_DT <- read.delim(init_tadFile, header=F, col.names = c("chromo","region", "start", "end"), stringsAsFactors = FALSE)
# init_tad_DT <- init_tad_DT[grepl("_TAD", init_tad_DT$region),]
# stopifnot(nrow(tad_DT) > 0)


# > init_tad_DT[8944:8948,]
# chromo      region     start       end
# 9074   chr6 chr6_TAD474 119080001 119160000 # midpos 119120000
# 9075   chr6 chr6_TAD475 119160001 119240000 # midpos 119200000
# 9077   chr6 chr6_TAD476 119280001 119720000 # midpos 119500000
# 9078   chr6 chr6_TAD477 119720001 120000000
# 
# 
# > tad_DT[8927:8930,]
# chromo      region     start       end
# 8965   chr6 chr6_TAD474 119120001 119200000
# 8966   chr6 chr6_TAD475 119200001 119520000
# 8967   chr6 chr6_TAD476 119520001 119880000
# > 
#   
# 
#   entrezID      region init_region
# 2481   254394 chr6_TAD475 chr6_TAD474
# 2562    25842 chr6_TAD475 chr6_TAD475
# 7074    79632 chr6_TAD475 chr6_TAD476


# 
logfc <- get(load("PIPELINE/OUTPUT_FOLDER/Barutcu_MCF-10A_RANDOMMIDPOSDISC_40kb/TCGAbrca_lum_bas/3_runMeanTADLogFC/all_meanLogFC_TAD.Rdata"))
reglist <- get(load("CHECK_DIFFTADS_RANDOMMIDPOS/Barutcu_MCF-10A_RANDOMMIDPOS_40kb/TCGAbrca_lum_bas/pipeline_regionList.Rdata"))
stopifnot(setequal(names(logfc), reglist))


