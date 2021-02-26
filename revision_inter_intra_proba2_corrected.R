setDir <- "/media/electron"
setDir <- ""

# v1 normalization: zscore

# CORRECTED: before normalization and mean, add the 0 values in the count vector !

# Rscript revision_inter_intra_proba2_corrected.R
script_name="revision_inter_intra_proba2_corrected.R"
startTime <- Sys.time()
cat("> START ", script_name, "\n")


require(foreach)
require(doMC)
registerDoMC(50)

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

outFolder <- "REVISION_INTER_INTRA_PROBA2_CORRECTED"
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
    # all_chrs=all_chrs[1]
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
        
        fullNormCount <- as.numeric(scale(fullVect)) ### use fullVect instead of  x$count !! CORRECTED

        stopifnot(length(fullNormCount) == fullSize)
        
        normCount <- fullNormCount[1:sparseSize]
        stopifnot(length(normCount) == nrow(x))
        x$normCount <- normCount
        x
      }))
      stopifnot(!is.na(matNorm_dt$count))
      stopifnot(!is.na(matNorm_dt$normCount))
      # can produce Na if not enough value at one diagodist # not true for v2
      matNorm_dt <- na.omit(matNorm_dt)
      cat(paste0("... after discarding NA values: ", nrow(matNorm_dt), "/", nrow(mat_dt), "\n"))
      rm("mat_dt")
      
      i=1 # CORRECTED VERSION: DO NOT DIRECTLY TAKE THE MEAN -> BUT SUM AND DIVIDE BY TOTAL POSSIBLE # OF VALUES
      overtads_dt <- foreach(i = 1:nrow(tad_dt), .combine='rbind') %dopar% {
        
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
        
        if(i == nrow(tad_dt)) { # if last tad -> no next tad
          # store upper right corner, so binA in current, binB in next
          # next_inter_values <- NA
          # next_inter_normValues <- NA
          mean_inter_next <-  NA
          mean_inter_nextNorm <-  NA
        
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
        }
        
        if(i == 1) { # if 1st tad -> no previous tad 
          # prev_inter_values <- NA
          # prev_inter_normValues <- NA
          # 
          mean_inter_prev <- NA
          mean_inter_prevNorm <- NA
          
          
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
          
          
        }
        
      
        data.frame(
          hicds = hicds,
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
          
          stringsAsFactors = FALSE
          
        )
      } # end iterate over TADs
      cat(paste0("... done for chromo ", chromo,"\n"))
      close(zz)
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
