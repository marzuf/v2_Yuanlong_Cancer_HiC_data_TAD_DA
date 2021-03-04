dt=get(load("hic_dc_comp_exp_info.Rdata"))

source("revision_settings.R")

all_chrs <- paste0("chr", 1:22)

all_interactome_dt <- do.call(rbind, lapply(seq_along(dt), function(i) {
  
  i_dt <- dt[[i]]
  
  chromo <- gsub("(.+):.+:.+", "chr\\1", names(dt)[i])
  stopifnot(chromo %in% all_chrs)
  chromo_nbr <- as.numeric(gsub("(.+):.+:.+", "\\1", names(dt)[i]))
  stopifnot(!is.na(chromo_nbr))
  
  tad_start <- as.numeric(gsub(".+:(.+):.+", "\\1", names(dt)[i]))
  stopifnot(!is.na(tad_start))
  tad_end <- as.numeric(gsub(".+:.+:(.+)", "\\1", names(dt)[i]))
  stopifnot(!is.na(tad_end))
  
  stopifnot(tad_start == min(i_dt[["pos_start"]]))
  stopifnot(tad_end == max(i_dt[["pos_end"]]))
  stopifnot(chromo_nbr == i_dt[["chr_num"]])
  
  sub_dt <- unique(i_dt[,grepl("mean_p", colnames(i_dt)) | grepl("sig_ratio", colnames(i_dt))])
  stopifnot(nrow(sub_dt) == 1)                                        
  sub_dt$chromo <- chromo
  sub_dt$tad_start <- tad_start
  sub_dt$tad_end <- tad_end
  sub_dt
}))
  
all_pairs_cols <- c("ENCSR346DCU_LNCaP_40kb" = "LNCaP",
                    "GSE118514_RWPE1_40kb" = "RWPE1",
                    "GSE118514_22Rv1_40kb"="22RV1")

myds <- all_pairs[grepl("prad", basename(all_pairs))]
all_cmps <- c(myds,
              as.character(sapply(myds , function(x) file.path(basename(dirname(x)), dirname(dirname(x)), basename(x)))))
### PLOT: 

final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_dt <- get(load(final_table_file))
final_table_DT <- final_dt
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)
colnames(final_table_DT)[colnames(final_table_DT) == "start"] <- "tad_start"
colnames(final_table_DT)[colnames(final_table_DT) == "end"] <- "tad_end"
final_table_DT$chromo <- gsub("(chr.+)_TAD.+", "\\1", final_table_DT$region )

sub_final_table_DT <- final_table_DT[final_table_DT$hicds %in% dirname(dirname(all_cmps)) &
                                   final_table_DT$exprds %in% basename(all_cmps),
                                   ] 
stopifnot(nrow(sub_final_table_DT) > 0)

merged_dt <- merge(sub_final_table_DT, all_interactome_dt, by=c("chromo", "tad_start", "tad_end"))

stopifnot(!duplicated(merged_dt$regionID))

cmp=all_cmps[1]
for(cmp in all_cmps) {

  plot_dt <- merged_dt
  
  ref_hicds <- dirname(dirname(cmp))
  match_hicds <- basename(dirname(cmp))
  exprds <- basename(cmp)

  sub_plot_dt <- plot_dt[plot_dt$hicds == ref_hicds &
                           plot_dt$exprds == exprds, ]
  stopifnot(nrow(sub_plot_dt) > 0)
  
  stopifnot(ref_hicds %in% names(all_pairs_cols))
  stopifnot(match_hicds %in% names(all_pairs_cols))
  ref_hicds_col <- as.character(all_pairs_cols[paste0(ref_hicds)])
  match_hicds_col <- as.character(all_pairs_cols[paste0(match_hicds)])

  icol ="sig_ratio"
  for(icol in c("sig_ratio", "mean_p")) {
    
    sub_plot_dt$match_ref_diff <- sub_plot_dt[,paste0(icol, "_", match_hicds_col)] - 
                                    sub_plot_dt[,paste0(icol, "_", ref_hicds_col)] 
    
    myxlab <- paste0("match - ref ", icol)
    myylab <- "ref. TAD adj. Pval [-log10]"
      
    mySub <- paste0("ref = ", ref_hicds, "; match = ", match_hicds, "; ", exprds, "; # TADs = ", nrow(sub_plot_dt))
    
    myy <- -log10(sub_plot_dt$adjPvalComb)
    myx <- sub_plot_dt$match_ref_diff
    
    plot(x=myx,
         y=myy,
         cex.main=plotCex,
         cex.lab=plotCex,
         cex.axis=plotCex)
    
    
  }
  
}

# cmp=all_cmps[1]
# for(cmp in all_cmps) {
#   
#   ref_hicds <- dirname(dirname(cmp))
#   match_hicds <- basename(dirname(cmp))
#   exprds <- basename(cmp)
#   
#   stopifnot(ref_hicds %in% names(all_pairs_cols))
#   stopifnot(match_hicds %in% names(all_pairs_cols))
#   ref_hicds_col <- all_pairs_cols[paste0(ref_hicds)]
#   match_hicds_col <- all_pairs_cols[paste0(match_hicds)]
#   
#   sub_final_dt <- final_table_DT[final_table_DT$hicds == ref_hicds &
#                                   final_table_DT$exprds == exprds, ]
#   stopifnot(nrow(sub_final_dt) > 0)
#   
#   # find the matching TAD
#   foreach(i = 1:nrow(sub_final_dt), .combine='rbind') %dopar% {
#     
#     curr_start <- sub_final_dt$start[i]
#     curr_end <- sub_final_dt$end[i]
#     curr_chromo <- sub_final_dt$chromo[i]
#     
#     tad_interact_dt <- all_interactome_dt[all_interactome_dt$chromo == chromo & 
#                                             all_interactome_dt$tad_start == curr_start & 
#                                             all_interactome_dt$tad_end == curr_end ,
#                                             ]
#     tad_interact_dt <- unique(tad_interact_dt)
#     stopifnot(nrow(sub_interact_dt) == 1)
#     
#     
#     
#     
#   }
#   
# }

# # load the interactome data
# # assign the signif. interactions to the TADs
# 
# binSize <- 40*10^3
# 
# all_chrs <- paste0("chr", 1:22)
# 
# all_interactMatch_dt <- foreach(hicds = all_hicds, .combine= 'rbind') %dopar% {
#   
#   
#   interactome_dt <- paste0(">>> CHANGE HERE")
#   stopifnot(interactome_dt$binA %% 1 == 0)
#   stopifnot(interactome_dt$binB %% 1 == 0)
#   
#   # ensure only intra
#   stopifnot(interactome_dt$chromoA == interactome_dt$chromoB)
#   
#   
#   all_domains_dt <- read.delim(file.path(hicds, "genes2tad/all_assigned_regions.txt"), col.names=c("chromo", "region", "start", "end"), 
#                                header=FALSE, stringsAsFactors = FALSE)
#   ### keep only the TADs !!!
#   all_tads_dt <- all_domains_dt[grepl("_TAD", all_domains_dt$region),]
#   stopifnot(nrow(all_tads_dt) > 0)
#   
#   # convert to 0-based bin
#   all_tads_dt$startBin <- (all_tads_dt$start-1)/binSize
#   stopifnot(all_tads_dt$startBin %% 1 == 0)
#   all_tads_dt$endBin <- (all_tads_dt$end)/binSize-1
#   stopifnot(all_tads_dt$endBin %% 1 == 0)
#   
#   stopifnot(all_chrs %in% all_tads_dt$chromo)
#   
#   interactome_dt$hicds <- hicds
#   interactome_dt$binA_tadMatch <- NA
#   interactome_dt$binB_tadMatch <- NA
#   
#   hicds_interactMatch_dt  <- foreach(i = 1:nrow(interactome_dt), .combine='rbind') %dopar% {
#     
#     # find the region of the 1st bin
#     binA_tadMatch <- all_tads_dt$region[all_tads_dt$chromo <= interactome_dt$chromoA[i] & 
#                                           all_tads_dt$startBin <= interactome_dt$binA[i] & 
#                                           all_tads_dt$endBin >= interactome_dt$binA[i]]
#     
#     binB_tadMatch <- all_tads_dt$region[all_tads_dt$chromo <= interactome_dt$chromoB[i] & 
#                                           all_tads_dt$startBin <= interactome_dt$binB[i] & 
#                                           all_tads_dt$endBin >= interactome_dt$binB[i]]
#     
#     if(length(binA_tadMatch) == 0) binA_tadMatch <- NA
#     if(length(binB_tadMatch) == 0) binB_tadMatch <- NA
#     
#     interactome_dt$binA_tadMatch <- binA_tadMatch
#     interactome_dt$binB_tadMatch <- binB_tadMatch
#     
#   }
#   
#   
#   
# }
