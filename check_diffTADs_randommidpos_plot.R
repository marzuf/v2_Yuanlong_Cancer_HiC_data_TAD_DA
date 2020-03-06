
# Rscript check_diffTADs_randommidpos_plot.R

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

outFolder <- "CHECK_DIFFTADS_RANDOMMIDPOS_PLOT" 
dir.create(outFolder, recursive = TRUE)


all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[! (grepl("_RANDOMMIDPOS", all_hicds) | grepl("PERMUT", all_hicds)) ]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

  
all_hicds = "ENCSR489OCU_NCI-H460_40kb"

rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T", "RANDOMMIDPOSDISC" )


buildData <- TRUE

if(buildData){
  
  

hicds = all_hicds[1]
all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    
    hicds_perm_dt <-foreach(rd_patt = rd_patterns, .combine='rbind') %do% {
      
    
      cat("... start ", hicds, " - ", exprds, " - ", rd_patt, "\n")
      
      rd_hicds <- file.path(dirname(hicds), gsub("_40kb", paste0("_", rd_patt, "_40kb"), basename(hicds)))
      
      
      ### RETRIEVE THE GENE2TAD ASSIGNMENT - RANDOM
      rd_g2tFile <- file.path(pipFolder, rd_hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(rd_g2tFile))
      rd_g2t_DT <- read.delim(rd_g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      rd_g2t_DT$entrezID <- as.character(rd_g2t_DT$entrezID)
      
            
      ### RETRIEVE THE GENES USED IN THE PIPELINE - script0 - RANDOM
      stopifnot(dir.exists(file.path(pipOutFolder, rd_hicds, exprds, script0_name)))
      rd_geneListFile <- file.path(pipOutFolder, rd_hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(rd_geneListFile))
      rd_pipeline_geneList <- eval(parse(text = load(rd_geneListFile))) # not adjusted
      stopifnot(rd_pipeline_geneList %in% rd_g2t_DT$entrezID)    
      rd_g2t_DT <- rd_g2t_DT[rd_g2t_DT$entrezID %in% pipeline_geneList,]
      
      
    
    
      ### RETRIEVE THE GENE2TAD ASSIGNMENT
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      
      ### RETRIEVE THE GENES USED IN THE PIPELINE - script0
      stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
      geneListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneListFile))
      pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
      stopifnot(pipeline_geneList %in% g2t_DT$entrezID)    
      g2t_DT <- g2t_DT[g2t_DT$entrezID %in% pipeline_geneList,]
      
      
      rd_g2t <- by(rd_g2t_DT, rd_g2t_DT$region, function(x) sort(x$entrezID))
      
      obs_g2t <- by(g2t_DT, g2t_DT$region, function(x) sort(x$entrezID))
      
      
      # EXACT MATCH
      rd_obs_match <- rd_g2t  %in% obs_g2t
      
      # # NESTED MATCH:
      rd_obs_nestedMatch <- unlist(lapply(rd_g2t, function(x) any(unlist(lapply(obs_g2t, function(obs_tad) all(x %in% obs_tad))))))
      
      if(hicds == "" & exprds == "TCGAluad_norm_luad" & rd_patt =="RANDOMMIDPOSDISC") save(rd_obs_nestedMatch, file = "rd_obs_nestedMatch.Rdata", version=2)
      
      stopifnot(sum(rd_obs_match) == sum(!is.na(match(rd_g2t, obs_g2t))))
      
      data.frame(
        hicds = hicds,
        randomData= rd_patt,
        nObsTADs = length(obs_g2t),
        nRdTADs = length(rd_g2t),
        nSameTADs = sum(rd_obs_match),
        nNestedTADs = sum(rd_obs_nestedMatch),
        stringsAsFactors = FALSE
      )
    }
    hicds_perm_dt
  }
  hicds_dt
}
}

      
outFile <- file.path(outFolder, "all_result_dt.Rdata")
save(all_result_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))




