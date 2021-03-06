setDir <- "/media/electron"
setDir <- ""

stop("-- sure do not want to use corrected version ?\n")

# v2 normalization: divide by the median [because of outliers] of diago instead of zscore
# so then I can take the median -> no negative values -> can take ratio

# Rscript revision_inter_intra_proba2_v2.R
script_name="revision_inter_intra_proba2_v2.R"
startTime <- Sys.time()
cat("> START ", script_name, "\n")


require(foreach)
require(doMC)
registerDoMC(60)

binSize <- 40000
all_chrs <- paste0("chr", 1:22)
# all_chrs=all_chrs[1]

all_ds <-  c(
#  "Barutcu_MCF-10A_40kb"="AWS_Barutcu_MCF-10A",
#  "Barutcu_MCF-7_40kb"="AWS_Barutcu_MCF-7",
#  "ENCSR079VIJ_G401_40kb" ="mega_ENCSR079VIJ_G401",
#  "ENCSR312KHQ_SK-MEL-5_40kb"="mega_ENCSR312KHQ_SK-MEL-5",
#  "ENCSR401TBQ_Caki2_40kb"="mega_ENCSR401TBQ_Caki2",
#  "ENCSR504OTV_transverse_colon_40kb"="ENCSR504OTV_transverse_colon",
  "ENCSR549MGQ_T47D_40kb"="mega_ENCSR549MGQ_T47D",
  "ENCSR862OGI_RPMI-7951_40kb"="mega_ENCSR862OGI_RPMI-7951",
  "GSE105194_cerebellum_40kb"="mega_GSE105194_cerebellum",
  "GSE105194_spinal_cord_40kb"="mega_GSE105194_spinal_cord",
  "GSE105318_DLD1_40kb"="mega_GSE105318_DLD1", 
  "GSE109229_BT474_40kb"="GSE109229_BT474", 
  "GSE109229_SKBR3_40kb"="GSE109229_SKBR3",
  "GSE118588_Panc_beta_40kb"="GSE118588_Panc_beta",
  "GSE99051_786_O_40kb" = "GSE99051_786_O",
  "PA2_40kb"="Compendium_PA2",
  "PA3_40kb"="Compendium_PA3",
  "Panc1_rep12_40kb"="mega_Panc1_rep12",
  "Rao_HCT-116_2017_40kb"="AWS_Rao_HCT-116_2017",
  "K562_40kb"="AWS_K562",
  "HMEC_40kb"="AWS_HMEC"
)





# all_ds = all_ds[1]

buildTable <- TRUE

outFolder <- "REVISION_INTER_INTRA_PROBA2_V2"
dir.create(outFolder, recursive = TRUE)

i=1
if(buildTable) {
  all_inter_intra_dt <- foreach(i_ds = seq_along(all_ds), .combine='rbind') %do% {
    
    cell_line <- all_ds[i_ds]
    hicds <- names(all_ds)[i_ds]
    
    all_tads_dt <- read.delim(file.path(hicds, "genes2tad/all_assigned_regions.txt"), col.names=c("chromo", "region", "start", "end"), 
                              header=FALSE, stringsAsFactors = FALSE)
    all_tads_dt <- all_tads_dt[grepl("_TAD", all_tads_dt$region),]
    stopifnot(nrow(all_tads_dt) > 0)
    
    stopifnot(all_chrs %in% all_tads_dt$chromo)
    
    chromo = "chr21"
    all_chromo_dt <- foreach(chromo = all_chrs,.combine='rbind') %do% {
      
      
      cat(paste0("START - ", hicds, " - ", chromo, "\n"))
      
      ### PREPARE TAD DATA
      tad_dt <- all_tads_dt[all_tads_dt$chromo == chromo,]
      stopifnot(nrow(tad_dt) > 0)
      
      # convert to 0-based bin
      tad_dt$startBin <- (tad_dt$start-1)/binSize
      tad_dt$endBin <- (tad_dt$end)/binSize-1
      stopifnot(tad_dt$startBin %% 1 == 0)
      stopifnot(tad_dt$endBin %% 1 == 0)
      stopifnot(tad_dt$startBin <= tad_dt$endBin)
      stopifnot(nrow(tad_dt) > 0)
      
      ### PREPARE HIC DATA
      matfile <- file.path(setDir, "/mnt/ndata/Yuanlong/2.Results/1.Juicer",
                           cell_line, "contact_mat", paste0("mat_", chromo, "_", binSize/1000, "kb_ob.txt.gz"))
      
      
      mat_dt <- read.csv(gzfile(matfile,'rt')  ,header=FALSE, sep="\t", col.names=c("coordA", "coordB", "count")) 
      stopifnot(is.numeric(mat_dt$coordA))
      stopifnot(is.numeric(mat_dt$coordB))
      stopifnot(is.numeric(mat_dt$count))
      mat_dt$binA <- mat_dt$coordA/binSize
      mat_dt$binB <- mat_dt$coordB/binSize
      stopifnot(mat_dt$binA %% 1 == 0)
      stopifnot(mat_dt$binB %% 1 == 0)
      stopifnot(mat_dt$binA <= mat_dt$binB) ### check upper right stored
      mat_dt$diagoDist <- mat_dt$binB-mat_dt$binA
      init_nrow <- nrow(mat_dt)
      mat_dt <- na.omit(mat_dt)
      cat(paste0("... after discarding NA values: ", nrow(mat_dt), "/", init_nrow, "\n"))
      cat(paste0("... performing median normalization...\n"))
      matNorm_dt <- do.call(rbind, by(mat_dt, mat_dt$diagoDist, function(x) {x$normCount <- x$count/median(x$count); x})) ## CHANGE HERE v2 NORM
      stopifnot(!is.na(matNorm_dt$count))
      stopifnot(!is.na(matNorm_dt$normCount))
      # can produce Na if not enough value at one diagodist # not true in v2
      matNorm_dt <- na.omit(matNorm_dt)
      cat(paste0("... after discarding NA values: ", nrow(matNorm_dt), "/", nrow(mat_dt), "\n"))
      rm("mat_dt")
      
      i=1
      overtads_dt <- foreach(i = 1:nrow(tad_dt), .combine='rbind') %dopar% {
        
        t_region <- tad_dt$region[i]
        
        t_start <- tad_dt$startBin[i]
        t_end <- tad_dt$endBin[i]
        
        next_start <- tad_dt$startBin[i+1]
        next_end <- tad_dt$endBin[i+1]
        
        prev_start <- tad_dt$startBin[i-1]
        prev_end <- tad_dt$endBin[i-1]
        
        intra_values <- matNorm_dt$count[(t_start <= matNorm_dt$binA  & t_end >= matNorm_dt$binA) &
                                           (t_start <= matNorm_dt$binB  & t_end >= matNorm_dt$binB)]
        intra_normValues <- matNorm_dt$normCount[(t_start <= matNorm_dt$binA  & t_end >= matNorm_dt$binA) &
                                                   (t_start <= matNorm_dt$binB  & t_end >= matNorm_dt$binB)]
        stopifnot(!is.na(intra_values))
        stopifnot(!is.na(intra_normValues))
        
        if(i == nrow(tad_dt)) { # if last tad -> no next tad
          # store upper right corner, so binA in current, binB in next
          next_inter_values <- NA
          next_inter_normValues <- NA
        } else {
          next_inter_values <- matNorm_dt$count[(t_start <= matNorm_dt$binA  & t_end >= matNorm_dt$binA) &
                                                  (next_start <= matNorm_dt$binB  & next_end >= matNorm_dt$binB)]
          next_inter_normValues <- matNorm_dt$normCount[(t_start <= matNorm_dt$binA  & t_end >= matNorm_dt$binA) &
                                                          (next_start <= matNorm_dt$binB  & next_end >= matNorm_dt$binB)]
          stopifnot(!is.na(next_inter_values))
          stopifnot(!is.na(next_inter_normValues))
        }
        
        if(i == 1) { # if 1st tad -> no previous tad 
          prev_inter_values <- NA
          prev_inter_normValues <- NA
        } else {
          # store upper right corner, so binA in previous, binB in current
          prev_inter_values <- matNorm_dt$count[(prev_start <= matNorm_dt$binA  & prev_end >= matNorm_dt$binA) &
                                                  (t_start <= matNorm_dt$binB  & t_end >= matNorm_dt$binB)]
          prev_inter_normValues <- matNorm_dt$normCount[(prev_start <= matNorm_dt$binA  & prev_end >= matNorm_dt$binA) &
                                                          (t_start <= matNorm_dt$binB  & t_end >= matNorm_dt$binB)]
          
          stopifnot(!is.na(prev_inter_values))
          stopifnot(!is.na(prev_inter_normValues))
          
        }
        data.frame(
          hicds = hicds,
          region = t_region,
          binStart = t_start,
          binEnd = t_end,
          start = tad_dt$start[i],
          end = tad_dt$end[i],
          mean_intra = mean(intra_values),
          mean_intraNorm = mean(intra_normValues),
          mean_inter_prev = mean(prev_inter_values),
          mean_inter_prevNorm = mean(prev_inter_normValues),
          mean_inter_next = mean(next_inter_values),
          mean_inter_nextNorm = mean(next_inter_normValues),
          stringsAsFactors = FALSE
          
        )
      } # end iterate over TADs
      overtads_dt
    } # end iterate over chromo
    all_chromo_dt
    outFile <- file.path(outFolder, paste0(hicds, "_all_chromo_dt.Rdata"))
    save(all_chromo_dt, file=outFile, version=2)
    cat(paste0("... written: ", outFile, "\n"))
    all_chromo_dt
  } # end iterate over DS
  outFile <- file.path(outFolder, "all_inter_intra_dt.Rdata")
  save(all_inter_intra_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_inter_intra_dt.Rdata")
  all_inter_intra_dt <- get(load(outFile))
}






######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




# 
# # data from Rao et al. indicate bins by genome coordinates, they need to be turned into indeces
#   chr.data$binsA = chr.data$binsA/bin.size
#   chr.data$binsB = chr.data$binsB/bin.size
# }
# 
# chr.data$binsA = chr.data$binsA+1
# chr.data$binsB = chr.data$binsB+1
