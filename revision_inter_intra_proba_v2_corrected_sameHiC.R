setDir <- "/media/electron"
setDir <- ""

# v2 normalization: divide by the median [because of outliers] of diago instead of zscore
# so then I can take the mean -> no negative values -> can take ratio

# CORRECTED: before normalization and mean, add the 0 values in the count vector !

# Rscript revision_inter_intra_proba_v2_corrected_sameHiC.R
script_name="revision_inter_intra_proba_v2_corrected_sameHiC.R"
startTime <- Sys.time()
cat("> START ", script_name, "\n")


require(foreach)
require(doMC)
registerDoMC(40)

binSize <- 40000
all_chrs <- paste0("chr", 1:22)
# all_chrs=all_chrs[1]

tmp_pairs <- c(
  file.path("LG1_40kb" ,"LG2_40kb"),
  file.path("PA2_40kb" ,"PA3_40kb" ),
  file.path("Barutcu_MCF-7_40kb", "GSE109229_BT474_40kb"), 
  file.path("GSE118514_22Rv1_40kb", "ENCSR346DCU_LNCaP_40kb")
  )
rev_pairs <- as.character(sapply(tmp_pairs, function(x) file.path(basename(x), dirname(x))))
all_pairs <- c(tmp_pairs, rev_pairs)
stopifnot(!duplicated(all_pairs))

# all_pairs=all_pairs[1]

all_ds <-  c(
  "Barutcu_MCF-7_40kb"="AWS_Barutcu_MCF-7",
  "GSE109229_BT474_40kb"="GSE109229_BT474",
  "LG1_40kb" ="Compendium_LG1",
  "LG2_40kb"="Compendium_LG2",
  "PA2_40kb" ="Compendium_PA2",
  "PA3_40kb"="Compendium_PA3",
  "ENCSR346DCU_LNCaP_40kb"="mega_ENCSR346DCU_LNCaP",
  "GSE118514_22Rv1_40kb"="GSE118514_22Rv1"
)


buildTable <- TRUE

outFolder <- "REVISION_INTER_INTRA_PROBA_V2_CORRECTED_SAMEHIC"
dir.create(outFolder, recursive = TRUE)

i=1
if(buildTable) {
  all_inter_intra_dt <- foreach(i_ds = seq_along(all_pairs), .combine='rbind') %do% {
    
    hicds <- dirname(all_pairs[i_ds])
    matched_hicds <- basename(all_pairs[i_ds])
    
    stopifnot(hicds %in% names(all_ds))               
    stopifnot(matched_hicds %in% names(all_ds))               
    
    cell_line <- all_ds[hicds]
    matched_cell_line <- all_ds[matched_hicds]
    
    all_tads_dt <- read.delim(file.path(hicds, "genes2tad/all_assigned_regions.txt"), col.names=c("chromo", "region", "start", "end"), 
                              header=FALSE, stringsAsFactors = FALSE)
    all_tads_dt <- all_tads_dt[grepl("_TAD", all_tads_dt$region),]
    stopifnot(nrow(all_tads_dt) > 0)
    
    stopifnot(all_chrs %in% all_tads_dt$chromo)
    
    chromo = "chr21"
    # all_chrs=all_chrs[1]
    all_chromo_dt <- foreach(chromo = all_chrs,.combine='rbind') %dopar% {
      
      
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
      
      zz <- gzfile(matfile,'rt')
      mat_dt <- read.csv(zz  ,header=FALSE, sep="\t", col.names=c("coordA", "coordB", "count")) 
      stopifnot(is.numeric(mat_dt$coordA))
      stopifnot(is.numeric(mat_dt$coordB))
      stopifnot(is.numeric(mat_dt$count))
      mat_dt$binA <- mat_dt$coordA/binSize  # 0-based bin
      mat_dt$binB <- mat_dt$coordB/binSize
      stopifnot(mat_dt$binA %% 1 == 0)
      stopifnot(mat_dt$binB %% 1 == 0)
      stopifnot(mat_dt$binA <= mat_dt$binB) ### check upper right stored
      mat_dt$diagoDist <- mat_dt$binB-mat_dt$binA
      
      mainDiagoSize <- max(mat_dt$binB)+1
      # mat_dt$diagoSize <- mainDiagoSize - mat_dt$diagoDist 
      
      init_nrow <- nrow(mat_dt)
      mat_dt <- na.omit(mat_dt)
      cat(paste0("... after discarding NA values: ", nrow(mat_dt), "/", init_nrow, "\n"))
      cat(paste0("... performing median normalization...\n"))
      # CORRECTED: CHANGED HERE -> add the 0s to have the full vector - 26.02.21
      matNorm_dt <- do.call(rbind, by(mat_dt, mat_dt$diagoDist, function(x) {
        diagDist <- as.numeric(unique(as.character(x$diagoDist)))
        stopifnot(length(diagDist) ==1)
        stopifnot(!is.na(diagDist))
        
        sparseVect <- x$count
        sparseSize <- length(sparseVect)
        fullSize <- mainDiagoSize - diagDist # how long should be the diago vect; to fill with 0s
        
        stopifnot(sparseSize <= fullSize)
        
        fullVect <- sparseVect
        if(fullSize > sparseSize) fullVect[(sparseSize+1):fullSize] <- 0
        
        stopifnot(length(fullVect) == fullSize)
        
        x$normCount <- x$count/median(fullVect)  ### here corrected -> compute median based on fullVect
        x
        })) ## CHANGE HERE v2 NORM
      stopifnot(!is.na(matNorm_dt$count))
      stopifnot(!is.na(matNorm_dt$normCount))
      # can produce Na if not enough value at one diagodist # not true for v2
      matNorm_dt <- na.omit(matNorm_dt)
      cat(paste0("... after discarding NA values: ", nrow(matNorm_dt), "/", nrow(mat_dt), "\n"))
      rm("mat_dt")
      
      
      ### PREPARE HIC DATA - for the matched data
      matched_matfile <- file.path(setDir, "/mnt/ndata/Yuanlong/2.Results/1.Juicer",
                                matched_cell_line, "contact_mat", paste0("mat_", chromo, "_", binSize/1000, "kb_ob.txt.gz"))
      
      stopifnot(matched_matfile != matfile)
      
      zz2 <- gzfile(matched_matfile,'rt')
      matched_mat_dt <- read.csv(zz2 ,header=FALSE, sep="\t", col.names=c("coordA", "coordB", "count")) 
      stopifnot(is.numeric(matched_mat_dt$coordA))
      stopifnot(is.numeric(matched_mat_dt$coordB))
      stopifnot(is.numeric(matched_mat_dt$count))
      matched_mat_dt$binA <- matched_mat_dt$coordA/binSize  # 0-based bin
      matched_mat_dt$binB <- matched_mat_dt$coordB/binSize
      stopifnot(matched_mat_dt$binA %% 1 == 0)
      stopifnot(matched_mat_dt$binB %% 1 == 0)
      stopifnot(matched_mat_dt$binA <= matched_mat_dt$binB) ### check upper right stored
      matched_mat_dt$diagoDist <- matched_mat_dt$binB-matched_mat_dt$binA
      
      mainDiagoSize <- max(matched_mat_dt$binB)+1
      # matched_mat_dt$diagoSize <- mainDiagoSize - matched_mat_dt$diagoDist 
      
      init_nrow <- nrow(matched_mat_dt)
      matched_mat_dt <- na.omit(matched_mat_dt)
      cat(paste0("... after discarding NA values: ", nrow(matched_mat_dt), "/", init_nrow, "\n"))
      cat(paste0("... performing median normalization...\n"))
      # CORRECTED: CHANGED HERE -> add the 0s to have the full vector - 26.02.21
      matched_matNorm_dt <- do.call(rbind, by(matched_mat_dt, matched_mat_dt$diagoDist, function(x) {
        diagDist <- as.numeric(unique(as.character(x$diagoDist)))
        stopifnot(length(diagDist) ==1)
        stopifnot(!is.na(diagDist))
        
        sparseVect <- x$count
        sparseSize <- length(sparseVect)
        fullSize <- mainDiagoSize - diagDist # how long should be the diago vect; to fill with 0s
        
        stopifnot(sparseSize <= fullSize)
        
        fullVect <- sparseVect
        if(fullSize > sparseSize) fullVect[(sparseSize+1):fullSize] <- 0
        
        stopifnot(length(fullVect) == fullSize)
        
        x$normCount <- x$count/median(fullVect)  ### here corrected -> compute median based on fullVect
        x
      })) ## CHANGE HERE v2 NORM
      stopifnot(!is.na(matched_matNorm_dt$count))
      stopifnot(!is.na(matched_matNorm_dt$normCount))
      # can produce Na if not enough value at one diagodist # not true for v2
      matched_matNorm_dt <- na.omit(matched_matNorm_dt)
      cat(paste0("... after discarding NA values: ", nrow(matched_matNorm_dt), "/", nrow(matched_mat_dt), "\n"))
      rm("matched_mat_dt")
      
      
      
      
      
      # if tad longer than matched matrix, will put 0

      
      i=1 # CORRECTED VERSION: DO NOT DIRECTLY TAKE THE MEAN -> BUT SUM AND DIVIDE BY TOTAL POSSIBLE # OF VALUES
      overtads_dt <- foreach(i = 1:nrow(tad_dt), .combine='rbind') %do% {
        
        t_region <- tad_dt$region[i]
        
        t_start <- tad_dt$startBin[i]
        t_end <- tad_dt$endBin[i]
        
        # if bin start = 0 and bin end = 3 -> submatrix of 4x4 => should have 10 values / for 5 and 7 -> 3x3, 6 values
        t_size <- t_end - t_start + 1
        nvalues_current <- floor(t_size * (t_size-1) * 0.5) + t_size  # need the floor for the case of size 1 = 1 value
        
        next_start <- tad_dt$startBin[i+1]
        next_end <- tad_dt$endBin[i+1]
        next_size <- next_end - next_start + 1
        # nvalues_next <- floor(next_size * (next_size-1) * 0.5) + next_size  # need the floor for the case of size 1 = 1 value
        
        
        prev_start <- tad_dt$startBin[i-1]
        prev_end <- tad_dt$endBin[i-1]
        prev_size <- prev_end - prev_start + 1
        # nvalues_prev <- floor(prev_size * (prev_size-1) * 0.5) + prev_size  # need the floor for the case of size 1 = 1 value
        
        
        # nbr values in inter: if a first TAD goes 0-3 and the next TAD goes 5-7 -> max # values = 12
        nvalues_interNext <- t_size * next_size
        nvalues_interPrev <- t_size * prev_size
        
        intra_values <- matNorm_dt$count[(t_start <= matNorm_dt$binA  & t_end >= matNorm_dt$binA) &
                                           (t_start <= matNorm_dt$binB  & t_end >= matNorm_dt$binB)]
        intra_normValues <- matNorm_dt$normCount[(t_start <= matNorm_dt$binA  & t_end >= matNorm_dt$binA) &
                                                   (t_start <= matNorm_dt$binB  & t_end >= matNorm_dt$binB)]
        stopifnot(!is.na(intra_values))
        stopifnot(!is.na(intra_normValues))
        
        stopifnot(length(intra_values) <= nvalues_current)
        stopifnot(length(intra_values) == length(intra_normValues))
        
        
        ##### do the same for matched hi-c
        matched_intra_values <- matched_matNorm_dt$count[(t_start <= matched_matNorm_dt$binA  & t_end >= matched_matNorm_dt$binA) &
                                                           (t_start <= matched_matNorm_dt$binB  & t_end >= matched_matNorm_dt$binB)]
        matched_intra_normValues <- matched_matNorm_dt$normCount[(t_start <= matched_matNorm_dt$binA  & t_end >= matched_matNorm_dt$binA) &
                                                                   (t_start <= matched_matNorm_dt$binB  & t_end >= matched_matNorm_dt$binB)]
        stopifnot(!is.na(matched_intra_values))
        stopifnot(!is.na(matched_intra_normValues))
        
        stopifnot(length(matched_intra_values) <= nvalues_current)
        stopifnot(length(matched_intra_values) == length(matched_intra_normValues))
        
        
        
        if(i == nrow(tad_dt)) { # if last tad -> no next tad
          # store upper right corner, so binA in current, binB in next
          # next_inter_values <- NA
          # next_inter_normValues <- NA
          mean_inter_next <-  NA
          mean_inter_nextNorm <-  NA
          
          matched_mean_inter_next <-  NA
          matched_mean_inter_nextNorm <-  NA
          
        
        } else {
          next_inter_values <- matNorm_dt$count[(t_start <= matNorm_dt$binA  & t_end >= matNorm_dt$binA) &
                                                  (next_start <= matNorm_dt$binB  & next_end >= matNorm_dt$binB)]
          next_inter_normValues <- matNorm_dt$normCount[(t_start <= matNorm_dt$binA  & t_end >= matNorm_dt$binA) &
                                                          (next_start <= matNorm_dt$binB  & next_end >= matNorm_dt$binB)]
          stopifnot(!is.na(next_inter_values))
          stopifnot(!is.na(next_inter_normValues))
          stopifnot(length(next_inter_values) <= nvalues_interNext)
          stopifnot(length(next_inter_values) == length(next_inter_normValues))
          
          mean_inter_next <-  sum(next_inter_values)/nvalues_interNext
          mean_inter_nextNorm <-  sum(next_inter_normValues)/nvalues_interNext
          
          #### do the same for the matched hi-c
          
          
          matched_next_inter_values <- matched_matNorm_dt$count[(t_start <= matched_matNorm_dt$binA  & t_end >= matched_matNorm_dt$binA) &
                                                                  (next_start <= matched_matNorm_dt$binB  & next_end >= matched_matNorm_dt$binB)]
          matched_next_inter_normValues <- matched_matNorm_dt$normCount[(t_start <= matched_matNorm_dt$binA  & t_end >= matched_matNorm_dt$binA) &
                                                                          (next_start <= matched_matNorm_dt$binB  & next_end >= matched_matNorm_dt$binB)]
          stopifnot(!is.na(matched_next_inter_values))
          stopifnot(!is.na(matched_next_inter_normValues))
          stopifnot(length(matched_next_inter_values) <= nvalues_interNext)
          stopifnot(length(matched_next_inter_values) == length(matched_next_inter_normValues))
          
          matched_mean_inter_next <-  sum(matched_next_inter_values)/nvalues_interNext
          matched_mean_inter_nextNorm <-  sum(matched_next_inter_normValues)/nvalues_interNext
          
        }
        
        if(i == 1) { # if 1st tad -> no previous tad 
          # prev_inter_values <- NA
          # prev_inter_normValues <- NA
          # 
          mean_inter_prev <- NA
          mean_inter_prevNorm <- NA
          
          matched_mean_inter_prev <- NA
          matched_mean_inter_prevNorm <- NA
          
          
        } else {
          # store upper right corner, so binA in previous, binB in current
          prev_inter_values <- matNorm_dt$count[(prev_start <= matNorm_dt$binA  & prev_end >= matNorm_dt$binA) &
                                                  (t_start <= matNorm_dt$binB  & t_end >= matNorm_dt$binB)]
          prev_inter_normValues <- matNorm_dt$normCount[(prev_start <= matNorm_dt$binA  & prev_end >= matNorm_dt$binA) &
                                                          (t_start <= matNorm_dt$binB  & t_end >= matNorm_dt$binB)]
          
          stopifnot(!is.na(prev_inter_values))
          stopifnot(!is.na(prev_inter_normValues))
          stopifnot(length(prev_inter_values) <= nvalues_interPrev)
          stopifnot(length(prev_inter_values) == length(prev_inter_normValues))
          
          
          mean_inter_prev <- sum(prev_inter_values)/nvalues_interPrev
          mean_inter_prevNorm <- sum(prev_inter_normValues)/nvalues_interPrev
          
          
          ### do the same for matched hi-c
          matched_prev_inter_values <- matched_matNorm_dt$count[(t_start <= matched_matNorm_dt$binA  & t_end >= matched_matNorm_dt$binA) &
                                                                  (prev_start <= matched_matNorm_dt$binB  & prev_end >= matched_matNorm_dt$binB)]
          matched_prev_inter_normValues <- matched_matNorm_dt$normCount[(t_start <= matched_matNorm_dt$binA  & t_end >= matched_matNorm_dt$binA) &
                                                                          (prev_start <= matched_matNorm_dt$binB  & prev_end >= matched_matNorm_dt$binB)]
          stopifnot(!is.na(matched_prev_inter_values))
          stopifnot(!is.na(matched_prev_inter_normValues))
          stopifnot(length(matched_prev_inter_values) <= nvalues_interPrev)
          stopifnot(length(matched_prev_inter_values) == length(matched_prev_inter_normValues))
          
          matched_mean_inter_prev <-  sum(matched_prev_inter_values)/nvalues_interPrev
          matched_mean_inter_prevNorm <-  sum(matched_prev_inter_normValues)/nvalues_interPrev
          
          
          
        }
        
      
        data.frame(
          hicds = hicds,
          matched_hicds = matched_hicds,
          region = t_region,
          binStart = t_start,
          binEnd = t_end,
          start = tad_dt$start[i],
          end = tad_dt$end[i],
          #### CHANGE HERE CORRECTED VERSION
          # mean_intra = mean(intra_values),
          # mean_intraNorm = mean(intra_normValues),
          # mean_inter_prev = mean(prev_inter_values),
          # mean_inter_prevNorm = mean(prev_inter_normValues),
          # mean_inter_next = mean(next_inter_values),
          # mean_inter_nextNorm = mean(next_inter_normValues),
          mean_intra = sum(intra_values)/nvalues_current,
          mean_intraNorm = sum(intra_normValues)/nvalues_current,
          
          mean_inter_prev = mean_inter_prev,
          mean_inter_prevNorm = mean_inter_prevNorm,
          
          mean_inter_next = mean_inter_next,
          mean_inter_nextNorm = mean_inter_nextNorm,
          
          # add the matched values
          matched_mean_intra = sum(matched_intra_values)/nvalues_current,
          matched_mean_intraNorm = sum(matched_intra_normValues)/nvalues_current,
          
          matched_mean_inter_prev = matched_mean_inter_prev,
          matched_mean_inter_prevNorm = matched_mean_inter_prevNorm,
          
          matched_mean_inter_next = matched_mean_inter_next,
          matched_mean_inter_nextNorm = matched_mean_inter_nextNorm,
          
          
          
          stringsAsFactors = FALSE
          
        )
      } # end iterate over TADs
      cat(paste0("... done for chromo ", chromo,"\n"))
      close(zz)
      close(zz2)
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
