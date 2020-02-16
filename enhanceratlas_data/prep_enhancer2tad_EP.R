# Rscript  prep_enhancer2tad_EP.R

require(foreach)
require(doMC)
registerDoMC(40)

outFolder <- file.path("PREP_ENHANCER2TAD_EP")
dir.create(outFolder, recursive = TRUE)

ep_dt <- read.delim("enhanceratlas_AllEPs_Lung_EP.txt", sep="$", stringsAsFactors = FALSE)
ep_dt <- ep_dt[,1:2]
colnames(ep_dt) <- c("enhancer", "target_symbol")
ep_dt$enhancer_chromo <- gsub("(.+):.+-.+_ENSG.+","\\1", ep_dt$enhancer )
ep_dt$enhancer_start <- as.numeric(gsub(".+:(.+)-.+_ENSG.+","\\1", ep_dt$enhancer ))
stopifnot(!is.na(ep_dt$enhancer_start))
ep_dt$enhancer_end <- gsub(".+:.+-(.+)_ENSG.+","\\1", ep_dt$enhancer )
stopifnot(!is.na(ep_dt$enhancer_end))
ep_dt$target_ensembl <- gsub(".+:.+-.+_(ENSG.+)","\\1", ep_dt$enhancer )

enhancer_dt <- ep_dt

cat(paste0(nrow(ep_dt), " \n"))

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
      tad_dt$chromo == enhancer_dt$enhancer_chromo[i] &
      tad_dt$start <= enhancer_dt$enhancer_start[i]  &
      tad_dt$end >= enhancer_dt$enhancer_start[i] 
      )
      if(length(i_match) == 0) {
        stopifnot(max(tad_dt$end[tad_dt$chromo == enhancer_dt$enhancer_chromo[i]]) < enhancer_dt$enhancer_start[i])
        match_region <- NA
      }else{
        match_region <- tad_dt$region[i_match]
      }
    
      data.frame(
        enhancer_chromo = enhancer_dt$enhancer_chromo[i],
        enhancer_start = enhancer_dt$enhancer_start[i],
        enhancer_end = enhancer_dt$enhancer_end[i],
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