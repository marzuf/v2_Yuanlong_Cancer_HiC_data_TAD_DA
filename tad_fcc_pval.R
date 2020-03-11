
# Rscript tad_fcc_pval.R

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "tad_fcc_pval.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"

require(foreach)
require(doMC)
require(ggpubr)
registerDoMC(40)

plotType <- "png"
myHeight <- 400
myWidth <- 400
plotCex <- 1.4

pipFolder<- file.path(".")
stopifnot(dir.exists(pipFolder))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "TAD_FCC_PVAL" 
dir.create(outFolder, recursive = TRUE)

all_hicds <- list.files(pipOutFolder)
# all_hicds <- all_hicds[! (grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds)) ]
# all_hicds <- all_hicds[grepl("ENCSR489OCU_NCI-H460_40kb", all_hicds)]
# all_hicds = "ENCSR489OCU_NCI-H460_40kb"
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))


rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMMIDPOSDISC" , "RANDOMMIDPOSSTRICT", "RANDOMSHIFT", "PERMUTG2T")
rd_patt = rd_patterns[1]

buildData <- TRUE

if(buildData){
  hicds = all_hicds[1]
  hicds = all_hicds[2]
  all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
        cat("... start ", hicds, " - ", exprds, " \n")
        
        rd_hicds <- file.path(dirname(hicds), gsub("_40kb", paste0("_", rd_patt, "_40kb"), basename(hicds)))
        
        
        fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCC_runAllDown", "all_obs_prodSignedRatio.Rdata")
        if(!file.exists(fcc_file)) return(NULL)
        
        pval_file <- file.path(pipOutFolder, hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")
        if(!file.exists(pval_file)) return(NULL)
        
        all_fcc <- get(load(fcc_file))
        all_pval <- get(load(pval_file))
        all_adjPval <- p.adjust(all_pval, method = "BH")
        stopifnot(names(all_adjPval) == names(all_fcc))
        
        data.frame(
          hicds = hicds,
          exprds=exprds,
          FCC = all_fcc,
          adjCombPval = all_adjPval,
          stringsAsFactors = FALSE
        )
        
        
        
    }
    hicds_dt
  }
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  save(all_result_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  all_result_dt <- get(load(outFile))
}  
all_result_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                                   gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                        gsub(".+RANDOMMIDPOSSTRICT_40kb", "RANDOMMIDPOSSTRICT",
                                             gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                                  gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                                       gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_result_dt$hicds))))))
all_result_dt$hicds_lab[! all_result_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"

all_types <- unique(all_result_dt$hicds_lab)

for(a_t in all_types) {
  
  plot_dt <- all_result_dt[all_result_dt$hicds_lab == a_t,]
  stopifnot(nrow(plot_dt) > 0)
  
  nDS <- length(unique(file.path(plot_dt$hicds, plot_dt$exprds)))
  
  my_x <- plot_dt$adjCombPval
  my_y <- plot_dt$FCC
  
  x_lab <- "TAD adj. comb. p-val [log10]"
  y_lab <- "TAD FCC"
  
  outFile <- file.path(outFolder, paste0("FCC_vs_adjCombPval_", a_t, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  densplot(
    x=my_x,
    y=my_y,
    xlab=x_lab,
    ylab=y_lab,
    main = paste0("FCC vs. p-val - ", a_t),
    cex.main=plotCex,
    cex.lab = plotCex,
    cex.axis = plotCex,
    cex=0.7
    
  )
  mtext(side=3, text = paste0("all DS (n=", nDS, ")"))
  addCorr(x=my_x, y=my_y, bty="n")
  
}

  
  
        