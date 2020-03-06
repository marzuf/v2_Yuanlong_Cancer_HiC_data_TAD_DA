
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
require(ggpubr)
registerDoMC(40)

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 10

pipFolder<- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "CHECK_DIFFTADS_RANDOMMIDPOS_PLOT" 
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files(pipOutFolder)
all_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds)) ]
# all_hicds <- all_hicds[grepl("ENCSR489OCU_NCI-H460_40kb", all_hicds)]
# all_hicds = "ENCSR489OCU_NCI-H460_40kb"

all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))

  
rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T", "RANDOMMIDPOSDISC" )
rd_patt = rd_patterns[1]

buildData <- FALSE

if(buildData){
  
  

hicds = all_hicds[1]
hicds = all_hicds[2]
all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    rd_patt = rd_patterns[1]
    hicds_perm_dt <-foreach(rd_patt = rd_patterns, .combine='rbind') %do% {
      
    
      cat("... start ", hicds, " - ", exprds, " - ", rd_patt, "\n")
      
      rd_hicds <- file.path(dirname(hicds), gsub("_40kb", paste0("_", rd_patt, "_40kb"), basename(hicds)))
      
      
      ### RETRIEVE THE GENE2TAD ASSIGNMENT - RANDOM
      rd_g2tFile <- file.path(pipFolder, rd_hicds, "genes2tad", "all_genes_positions.txt")
      
      if(!file.exists(rd_g2tFile)) save(rd_g2tFile, file ="rd_g2tFile.Rdata", version=2)
      
      stopifnot(file.exists(rd_g2tFile))
      rd_g2t_DT <- read.delim(rd_g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      rd_g2t_DT$entrezID <- as.character(rd_g2t_DT$entrezID)
      
            
      ### RETRIEVE THE GENES USED IN THE PIPELINE - script0 - RANDOM
      stopifnot(dir.exists(file.path(pipOutFolder, rd_hicds, exprds, script0_name)))
      rd_geneListFile <- file.path(pipOutFolder, rd_hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(rd_geneListFile))
      rd_pipeline_geneList <- eval(parse(text = load(rd_geneListFile))) # not adjusted
      stopifnot(rd_pipeline_geneList %in% rd_g2t_DT$entrezID)    
      rd_g2t_DT <- rd_g2t_DT[rd_g2t_DT$entrezID %in% rd_pipeline_geneList,]
      
    
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
      
      
      # rd_g2t <- list(
      #   c("a", "b", "c", "d"),
      #   c("e", "f"),
      #   c("g", "h", "i", "j")
      # )
      # 
      # obs_g2t <- list(
      #   c("a", "d", "i"),
      #   c("e", "f", "g")
      # )
      # 
      # unlist(lapply(rd_g2t, function(x) any(unlist(lapply(obs_g2t, function(obs_tad) all(x %in% obs_tad))))))
      
      
      if(hicds == "ENCSR489OCU_NCI-H460_40kb" & exprds == "TCGAluad_norm_luad" & rd_patt =="RANDOMMIDPOSDISC") 
          save(rd_obs_nestedMatch, file = "rd_obs_nestedMatch.Rdata", version=2)
      
      stopifnot(sum(rd_obs_match) == sum(!is.na(match(rd_g2t, obs_g2t))))
      
      data.frame(
        hicds = hicds,
        exprds =exprds,
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
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  save(all_result_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  all_result_dt <- get(load(outFile))
}

nDS <- length(unique(file.path(all_result_dt$hicds, all_result_dt$exprds)))

mT = c("Nested","Same")[1]
for(mT in c("Nested", "Same")) {
  
  all_result_dt$ratio <- all_result_dt[,paste0("n", mT, "TADs")]/all_result_dt$nRdTADs
  
  box_ratio <- ggboxplot(all_result_dt, x="randomData", y="ratio",
            xlab = "", 
            ylab = paste0("ratio ", mT)) +
    labs(title = paste0("ratio ", mT ," with obs. data"), subtitle = paste0("all DS - n=", nDS))+
    scale_y_continuous(
    breaks = scales::pretty_breaks(n = 10))+
  theme( # Increase size of axis lines
    strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"))

outFile <- file.path(outFolder, paste0("allDS_ratio", mT, "TAD_boxplot.", plotType))
ggsave(box_ratio, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


  
  
}




