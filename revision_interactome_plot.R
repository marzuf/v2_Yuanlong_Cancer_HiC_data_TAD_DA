# hicds             exprds       region meanLogFC  adjPvalComb
# 1 GSE118514_RWPE1_40kb TCGAprad_norm_prad chr12_TAD194  3.354541 1.415930e-06
# 2 GSE118514_RWPE1_40kb TCGAprad_norm_prad  chr7_TAD424 -1.552305 3.467131e-05
# 3 GSE118514_RWPE1_40kb TCGAprad_norm_prad chr17_TAD174  1.960796 5.311100e-05

# hicds             exprds       region     start       end
# 1 GSE118514_RWPE1_40kb TCGAprad_norm_prad chr12_TAD194  54160001  54440000
# 2 GSE118514_RWPE1_40kb TCGAprad_norm_prad  chr7_TAD424 116080001 116320000
# 3 GSE118514_RWPE1_40kb TCGAprad_norm_prad chr17_TAD174  46720001  46880000

# tads_to_plot <- list(
#   list(chromo="chr12", start=54160001, end= 54440000),
#   list(chromo="chr7", start =116080001 , end=116320000),
#   list(chromo="chr17", start = 46720001, end=46880000)
# )
# 

# Rscript revision_interactome_plot.R

setDir <- "/media/electron"
setDir <-  ""

outFolder <- file.path("REVISION_INTERACTOME_PLOT")
dir.create(outFolder, recursive = TRUE)

myHeight <- myWidth <- 400
plotType <- "png"

options(scipen = 100)

runFolder <- "."
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

hicds_of_interest <- c( "GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "GSE118514_22Rv1_40kb")
# hicds_of_interest <- c( "GSE118514_RWPE1_40kb")
exprds_of_interest <- c("TCGAprad_norm_prad")

dt <- resultData[resultData$hicds %in% hicds_of_interest & resultData$exprds %in% exprds_of_interest ,]
dt <- dt[order(dt$adjPvalComb),]
stopifnot(nrow(dt) > 0)

out_dt = do.call(rbind, by(dt, dt$dataset, function(x) {tmp=x[1:3,];tmp$tad_rank <- rank(tmp$adjPvalComb); tmp}))
out_dt$chromo <- gsub("(chr.+)_.+", "\\1", out_dt$region)
rownames(out_dt) <- NULL

all_ref_hicds <- unique(out_dt$hicds)
ref_hicds <- "GSE118514_RWPE1_40kb"
ref_hicds_lab <- "RWPE1"
all_matches <- c("LNCaP", "22Rv1")

all_ref_hicds_labs <- setNames( c("RWPE1", "LNCaP", "22Rv1"), 
                                c("GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "GSE118514_22Rv1_40kb"))
  
for(ref_hicds in all_ref_hicds) {
  
  cat(paste0("... start ", ref_hicds, "\n"))
  
  stopifnot(ref_hicds  %in% names(all_ref_hicds_labs))
  ref_hicds_lab <- as.character(all_ref_hicds_labs[ref_hicds])
  
  curr_out_dt <- out_dt[out_dt$hicds == ref_hicds,]
  
  print(curr_out_dt[,c("hicds", "exprds", "region", "adjPvalComb", "meanLogFC")])
  
  stopifnot(nrow(curr_out_dt) > 0)
  
  
  all_matches <- as.character(all_ref_hicds_labs[!all_ref_hicds_labs==ref_hicds_lab])
  stopifnot(length(all_matches) == length(all_ref_hicds_labs)-1)
  
  
  ref_data <- get(load(file.path(setDir,
                                 "/mnt/ndata/Yuanlong/1.Projects/19.With_Marie/1.Data",
                                 paste0(ref_hicds, "_TADs_for_hicdc.output.Rdata"))))
  
  
  for(i in 1:nrow(curr_out_dt)) {
    
    chromo <- curr_out_dt$chromo[i]
    tad_start <- curr_out_dt$start[i]
    tad_end <- curr_out_dt$end[i]
    region_genes <- curr_out_dt$region_genes[i]
    
    ref_tad <- paste0(chromo, ":", tad_start, "-", tad_end)
    
    cat(paste0(ref_hicds_lab, " - ", ref_tad, " \n"))
    
    tad_lab <- paste0(gsub("chr", "", chromo), ":", tad_start, ":", tad_end)
    stopifnot(tad_lab %in% names(ref_data))
    
    exp_dt <- ref_data[[tad_lab]][["summary_tab"]]
    exp_dt <- exp_dt[,c(which(colnames(exp_dt) == "gene"), grep("exp_", colnames(exp_dt)))]
    exp_dt <- exp_dt[exp_dt$gene != ".",]
    stopifnot(ncol(exp_dt) == 4)
    exp_dt[,2:4] <- round(exp_dt[,2:4],4)
    print(exp_dt)
    print(colSums(exp_dt[,2:4]))
    
    print_dt <- ref_data[[tad_lab]][["summary_tab"]]
    print_dt <- unique(print_dt[,grepl("mean_p", colnames(print_dt)) | grepl("sig_", colnames(print_dt))])
    print_dt <- round(print_dt,4)
    curr_mats <- ref_data[[tad_lab]][["hic_dc_mat"]]
    
    print(round(print_dt,4))
    
    ref_mat <- as.matrix(curr_mats[[ref_hicds_lab]])
    
    
    match_ds = all_matches[1]
    for(match_ds in all_matches) {
      stopifnot(ref_hicds_lab %in% names(curr_mats))
      stopifnot(match_ds %in% names(curr_mats))
      match_mat <- as.matrix(curr_mats[[match_ds]])
      
      outFile <- file.path(outFolder, paste0(ref_hicds_lab, "_", ref_tad, "_vs_", match_ds, ".", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
      par(mfrow=c(1,2))
      image(ref_mat, main=paste0(ref_hicds_lab, " ", ref_tad))
      mtext(side=3, text=paste0(ref_hicds_lab, " - ", region_genes))
      image(match_mat)
      mtext(side=3, text=paste0(match_ds))
      
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
    }
    
    
  }
  
  
  
  
}


