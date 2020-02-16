# Rscript  prep_enhancer2tad.R

require(foreach)
require(doMC)
registerDoMC(40)

outFolder <- file.path("PREP_ENHANCER2TAD")
dir.create(outFolder, recursive = TRUE)

enhancer_dt <- read.delim("enhanceratlas_enhancer_Lung.bed", col.names = c("chromo", "start", "end", "score"), stringsAsFactors = FALSE)

runFolder <- ".."

all_hicds <- list.files(file.path(runFolder, "PIPELINE/OUTPUT_FOLDER"))

myhicds <- "ENCSR489OCU_NCI-H460"

all_hicds <- all_hicds[grep(myhicds, all_hicds)]

all_exprds <- sapply(all_hicds, function(x) list.files(file.path("../PIPELINE/OUTPUT_FOLDER", x)))

foo <- foreach(hicds = all_hicds, .combine='rbind') %do%{
  
  cat(paste0("... start - ", hicds,  "\n"))
  
  
  # exprds = all_exprds[[paste0(hicds)]][1]
  # foo <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    # if(! (grepl("TCGAluad", exprds) | grepl("TCGAlusc", exprds))) return(NULL)
    
    # cat(paste0("... start - ", hicds, " - ", exprds, "\n"))
    
    tad_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    tad_dt <- read.delim(tad_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("chromo","region", "start", "end"))
    
    enhancer2tad_dt <- foreach(i = 1:nrow(enhancer_dt), .combine='rbind') %dopar% {
      # assign enhancer to TAD/BD according to start position
      
      i_match <- which(
      tad_dt$chromo == enhancer_dt$chromo[i] &
      tad_dt$start <= enhancer_dt$start[i]  &
      tad_dt$end >= enhancer_dt$start[i] 
      )
      if(length(i_match) == 0) {
        stopifnot(max(tad_dt$end[tad_dt$chromo == enhancer_dt$chromo[i]]) < enhancer_dt$start[i])
        match_region <- NA
      }else{
        match_region <- tad_dt$region[i_match]
      }
    
      data.frame(
        enhancer_chromo = enhancer_dt$chromo[i],
        enhancer_start = enhancer_dt$start[i],
        enhancer_end = enhancer_dt$end[i],
        region = match_region,
        stringsAsFactors = FALSE
      )
      
    }
      # outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_enhancer2tad_dt.Rdata"))
    outFile <- file.path(outFolder, paste0(hicds, "_enhancer2tad_dt.Rdata"))
    save(enhancer2tad_dt, file = outFile, version=2)
    cat(paste0("... written: ", outFile, "\n"))
    
  # } # end-foreach exprds
} # end-foreach hicds