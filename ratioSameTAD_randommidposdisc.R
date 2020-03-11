
# Rscript ratioSameTAD_randommidposdisc.R

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "ratioSameTAD_randommidposdisc.R"
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

outFolder <- file.path("RATIOSAMETAD_RANDOMMIDPOSDISC")
dir.create(outFolder, recursive = TRUE)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")


############################################ RANDOMMIDPOSDISC

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[grepl("_RANDOMMIDPOSDISC_", all_hicds) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


hicds = all_hicds[1]
all_result_randommidposdisc <- foreach(hicds = all_hicds, .combine='c') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_data <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
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
    init_hicds <- gsub("_RANDOMMIDPOSDISC_40kb", "_40kb", hicds)
    init_g2tFile <- file.path(pipFolder, init_hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(init_g2tFile))
    init_g2t_DT <- read.delim(init_g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "init_region"), stringsAsFactors = FALSE)
    init_g2t_DT$entrezID <- as.character(init_g2t_DT$entrezID)
    
    cmp_init_g2t_dt <- merge(g2t_DT[,c("entrezID", "region")], init_g2t_DT[,c("entrezID", "init_region")], all.x = TRUE, all.y=FALSE, by="entrezID")
    stopifnot(!is.na(cmp_init_g2t_dt))
    stopifnot(nrow(cmp_init_g2t_dt) == nrow(g2t_DT))
    
    cmp_init_g2t_dt_noBD <- cmp_init_g2t_dt[!grepl("_BOUND", cmp_init_g2t_dt$init_region),]
    initRegByReg_dt_noBD <- aggregate(init_region ~ region, data =cmp_init_g2t_dt_noBD, FUN = function(x) length(unique(x)))
    
    
    save(cmp_init_g2t_dt, file="cmp_init_g2t_dt.Rdata", version=2)
    
    # change here for strict version
    allRatioRegByReg <- foreach(reg = unique(cmp_init_g2t_dt$region), .combine='c') %do% {
      x <- cmp_init_g2t_dt$init_region[cmp_init_g2t_dt$region == reg]
      regRatio <- as.numeric(table(x))/length(x) # 1 = all the genes come from the same TAD; # 0.5 => max. 50% come from a same TAD
      stopifnot(sum(regRatio) == 1)
      regRatio
    }
    allRatioRegByReg
    
    save(allRatioRegByReg, file="allRatioRegByReg.Rdata", version=2)
    
    allRatioRegByReg
  }
  hicds_data
}

############################################ RANDOMMIDPOS

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[grepl("_RANDOMMIDPOS_", all_hicds) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


hicds = all_hicds[1]
all_result_randommidpos <- foreach(hicds = all_hicds, .combine='c') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_data <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
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
    
    
    save(cmp_init_g2t_dt, file="cmp_init_g2t_dt.Rdata", version=2)
    
    # change here for strict version
    allRatioRegByReg <- foreach(reg = unique(cmp_init_g2t_dt$region), .combine='c') %do% {
      x <- cmp_init_g2t_dt$init_region[cmp_init_g2t_dt$region == reg]
      regRatio <- as.numeric(table(x))/length(x) # 1 = all the genes come from the same TAD; # 0.5 => max. 50% come from a same TAD
      stopifnot(sum(regRatio) == 1)
      regRatio
    }
    allRatioRegByReg
    
    save(allRatioRegByReg, file="allRatioRegByReg.Rdata", version=2)
    
    allRatioRegByReg
  }
  hicds_data
}



############################################ RANDOMMIDPOSSTRICT
all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[grepl("_RANDOMMIDPOSSTRICT_", all_hicds) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


hicds = all_hicds[1]
all_result_randommidposstrict <- foreach(hicds = all_hicds, .combine='c') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_data <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
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
    init_hicds <- gsub("_RANDOMMIDPOSSTRICT_40kb", "_40kb", hicds)
    init_g2tFile <- file.path(pipFolder, init_hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(init_g2tFile))
    init_g2t_DT <- read.delim(init_g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "init_region"), stringsAsFactors = FALSE)
    init_g2t_DT$entrezID <- as.character(init_g2t_DT$entrezID)
    
    cmp_init_g2t_dt <- merge(g2t_DT[,c("entrezID", "region")], init_g2t_DT[,c("entrezID", "init_region")], all.x = TRUE, all.y=FALSE, by="entrezID")
    stopifnot(!is.na(cmp_init_g2t_dt))
    stopifnot(nrow(cmp_init_g2t_dt) == nrow(g2t_DT))
    
    cmp_init_g2t_dt_noBD <- cmp_init_g2t_dt[!grepl("_BOUND", cmp_init_g2t_dt$init_region),]
    initRegByReg_dt_noBD <- aggregate(init_region ~ region, data =cmp_init_g2t_dt_noBD, FUN = function(x) length(unique(x)))
    
    
    save(cmp_init_g2t_dt, file="cmp_init_g2t_dt.Rdata", version=2)
    
    # change here for strict version
    allRatioRegByReg <- foreach(reg = unique(cmp_init_g2t_dt$region), .combine='c') %do% {
      x <- cmp_init_g2t_dt$init_region[cmp_init_g2t_dt$region == reg]
      regRatio <- as.numeric(table(x))/length(x) # 1 = all the genes come from the same TAD; # 0.5 => max. 50% come from a same TAD
      stopifnot(sum(regRatio) == 1)
      regRatio
    }
    allRatioRegByReg
    
    save(allRatioRegByReg, file="allRatioRegByReg.Rdata", version=2)
    
    allRatioRegByReg
  }
  hicds_data
}


############################################ RANDOMMIDPOSDISC

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[grepl("_RANDOMMIDPOSDISC_", all_hicds) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


hicds = all_hicds[1]
all_result_randommidposdisc <- foreach(hicds = all_hicds, .combine='c') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_data <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
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
    init_hicds <- gsub("_RANDOMMIDPOSDISC_40kb", "_40kb", hicds)
    init_g2tFile <- file.path(pipFolder, init_hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(init_g2tFile))
    init_g2t_DT <- read.delim(init_g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "init_region"), stringsAsFactors = FALSE)
    init_g2t_DT$entrezID <- as.character(init_g2t_DT$entrezID)
    
    cmp_init_g2t_dt <- merge(g2t_DT[,c("entrezID", "region")], init_g2t_DT[,c("entrezID", "init_region")], all.x = TRUE, all.y=FALSE, by="entrezID")
    stopifnot(!is.na(cmp_init_g2t_dt))
    stopifnot(nrow(cmp_init_g2t_dt) == nrow(g2t_DT))
    
    cmp_init_g2t_dt_noBD <- cmp_init_g2t_dt[!grepl("_BOUND", cmp_init_g2t_dt$init_region),]
    initRegByReg_dt_noBD <- aggregate(init_region ~ region, data =cmp_init_g2t_dt_noBD, FUN = function(x) length(unique(x)))
    
    
    save(cmp_init_g2t_dt, file="cmp_init_g2t_dt.Rdata", version=2)
    
    # change here for strict version
    allRatioRegByReg <- foreach(reg = unique(cmp_init_g2t_dt$region), .combine='c') %do% {
      x <- cmp_init_g2t_dt$init_region[cmp_init_g2t_dt$region == reg]
      regRatio <- as.numeric(table(x))/length(x) # 1 = all the genes come from the same TAD; # 0.5 => max. 50% come from a same TAD
      stopifnot(sum(regRatio) == 1)
      regRatio
    }
    allRatioRegByReg
    
    save(allRatioRegByReg, file="allRatioRegByReg.Rdata", version=2)
    
    allRatioRegByReg
  }
  hicds_data
}
# outFile <- file.path(outFolder, "all_result.Rdata")
# save(all_result, file = outFile, version=2)
# cat(paste0("... written: ", outFile, "\n"))



outFile <- file.path(outFolder, "ratioTAD_within_TAD_all_randommidpos_density.svg")
do.call(svg, list(outFile, height=7, width=9))
plot_multiDens(
  list(randommidpos=all_result_randommidpos,
       randommidposdisc = all_result_randommidposdisc,
       randommidposstrict = all_result_randommidposstrict), legPos = "topright",
  plotTit = paste0("Distribution ratio init. TADs")
)
mtext(side=3, paste0("all DS"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

