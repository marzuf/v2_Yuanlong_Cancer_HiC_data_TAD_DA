# Rscript look_fithic_results.R ENCSR444WCZ_A549 chr10 20000

cl <- "ENCSR444WCZ_A549"
chromo <- "chr10"
binResol <- 10000

args <- commandArgs(trailingOnly = TRUE)
cl <- args[1]
chromo <- args[2]
binResol <- as.numeric(args[3])

stopifnot(!is.na(binResol))

tad_dt <- read.delim("../lung_hist_data/g38/LIFTOVER_COORD/AKR1C1_all_histOverlap_dt.txt", sep="\t", stringsAsFactors = FALSE, header=T)
tad_dt <- tad_dt[grepl(cl, tad_dt$hicds),]
tad_dt <- unique(tad_dt)
stopifnot(nrow(tad_dt) == 1)
tadStart <- tad_dt$hg19_tadStart_oi
tadEnd <- tad_dt$hg19_tadEnd_oi
geneName <- tad_dt$gene_oi
tadName <- tad_dt$tad_oi

fitoutFolder <- "FITHIC_OUTPUT"

outFolder <- file.path("LOOK_FITHIC_RESULTS",  paste0(cl, "_mat_", chromo, "_",binResol/1000, "kb"))

dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- 400
myWidth <- 400

sigThresh <- 0.05

for(i in 1:2){
  
  count_dist <- read.delim(file.path(fitoutFolder, paste0(cl, "_mat_", chromo, "_",binResol/1000, "kb"), paste0(".fithic_pass", i, ".txt")),header=T)
  count_dist$avgGenomicDist_log10 <- log10(count_dist$avgGenomicDist)
  # plot(contactProbability ~ avgGenomicDist, data=count_dist)
  
  outFile <- file.path(outFolder, paste0("proba_dist_fithic_pass", i, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(contactProbability ~ avgGenomicDist_log10, data=count_dist, pch=16, cex =0.7, main=paste0("fithic_pass", i))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  sig_file <- file.path(fitoutFolder,  paste0(cl, "_mat_", chromo, "_",binResol/1000, "kb"), paste0(".spline_pass", i, ".significances.txt.gz"))
  sig_dt <- read.csv(gzfile(sig_file), sep="\t", header=TRUE)
  
  tad_sig_dt <- sig_dt[(sig_dt$fragmentMid1 >= tadStart & sig_dt$fragmentMid1 <= tadEnd) &
            (sig_dt$fragmentMid2 >= tadStart & sig_dt$fragmentMid2 <= tadEnd),]
  tad_sig_dt <- na.omit(tad_sig_dt)
  nrow(tad_sig_dt)
  
  if(nrow(tad_sig_dt) > 0) {
    
    save(tad_sig_dt, file="tad_sig_dt.Rdata", version=2)
    
    stopifnot(is.numeric(tad_sig_dt$q_value))
    stopifnot(is.numeric(tad_sig_dt$p_value))
    
    tad_sig_dt <- tad_sig_dt[tad_sig_dt$p_value <= sigThresh | tad_sig_dt$q_value <= sigThresh,] 
    
    if(nrow(tad_sig_dt) > 0) {
    outFile <- file.path(outFolder, paste0("spline_pass", i, "_", geneName, "_", tadName, "_signif.txt"))
    write.table(tad_sig_dt, file =outFile, sep="\t", quote=F, col.names=TRUE, row.names=FALSE)
    cat(paste0("... written: ", outFile, "\n"))
    }
  }
  
  
}
