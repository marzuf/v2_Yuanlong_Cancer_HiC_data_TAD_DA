
# Rscript revision_interactome.R

require(foreach)
require(doMC)
registerDoMC(50)
require(ggplot2)
require(ggpubr)
require(ggsci)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

initbinsize <- 40000
hicdcbinsize <- 20000
signifThresh <- 0.05

plotType <- "png"
myHeight <- 400
myWidth <- 400
myWidthGG <- 7
myHeightGG <- 5
plotCex <- 1.2

buildTable <- F

tad_signifThresh <- 0.01


setDir <- "/media/electron"
setDir <- ""

all_pairs_cols <- c("ENCSR346DCU_LNCaP_40kb" = "LNCaP",
                    "GSE118514_RWPE1_40kb" = "RWPE1",
                    "GSE118514_22Rv1_40kb"="22Rv1")

outFolder <- file.path("REVISION_INTERACTOME")
dir.create(outFolder, recursive = TRUE)

all_chrs <- paste0("chr", 1:22)

ds=names(all_pairs_cols)[1]

if(buildTable){
  # all_hicdc_dt <- foreach(ds = c("GSE118514_22Rv1_40kb", .combine='rbind') %do% {
  
  # all_hicdc_dt <- foreach(ds = names(all_pairs_cols), .combine='rbind') %do% {
  for(ds in names(all_pairs_cols)) {
    
    cat(paste0("start ", ds, "\n"))
    
    dt <- get(load(file.path(setDir,
                             "/mnt/ndata/Yuanlong/1.Projects/19.With_Marie/1.Data",
                             paste0(ds, "_TADs_for_hicdc.output.Rdata"))))
    # level 1 - list of TADs,
    # level 2 - for a given TAD
    # - 2a "hic_dc_mat"  = list of 3 matrices ("LNCaP" "22Rv1" "RWPE1")
    # - 2b"summary_tab" = data frame with colnames:
    # [1] "gene"                          "LNCaP_prostate_cancer"        
    # [3] "22Rv1_prostate_cancer"         "RWPE1_prostate_epithelium"    
    # [5] "exp_LNCaP_prostate_cancer"     "exp_22Rv1_prostate_cancer"    
    # [7] "exp_RWPE1_prostate_epithelium" "sig_ratio_LNCaP"              
    # [9] "sig_ratio_22Rv1"               "sig_ratio_RWPE1"              
    # [11] "mean_p_LNCaP"                  "mean_p_22Rv1"                 
    # [13] "mean_p_RWPE1"                  "pos_start"                    
    # [15] "pos_end"                       "chr_num"                      
    
    # iterate over the TADs
    all_interactome_dt <- foreach(i= seq_along(dt), .combine='rbind') %dopar% {
      
      i_dt <- dt[[i]]
      
      if(i_dt[["hic_dc_mat"]] == "query_out_of_range") return(NULL)
      
      summary_dt <- i_dt[["summary_tab"]]
      
      chromo <- gsub("(.+):.+:.+", "chr\\1", names(dt)[i])
      stopifnot(chromo %in% all_chrs)
      chromo_nbr <- as.numeric(gsub("(.+):.+:.+", "\\1", names(dt)[i]))
      stopifnot(!is.na(chromo_nbr))
      
      tad_start <- as.numeric(gsub(".+:(.+):.+", "\\1", names(dt)[i]))
      stopifnot(!is.na(tad_start))
      tad_end <- as.numeric(gsub(".+:.+:(.+)", "\\1", names(dt)[i]))
      stopifnot(!is.na(tad_end))
      
      
      # stopifnot(tad_start == min(summary_dt[["pos_start"]]))# not always true, si e.g.
      # > dt=get(load("/mnt/ndata/Yuanlong/1.Projects/19.With_Marie/1.Data/GSE118514_22Rv1_40kb_TADs_for_hicdc.output.Rdata"))                                                                              
      # > names(dt)[9649]
      # [1] "20:29400001:29840000"
      # > min(dt[[9649]][["summary_tab"]][["pos_start"]])
      # [1] 29440001
      stopifnot(tad_start <= min(summary_dt[["pos_start"]]))# not always true, si e.g.
      stopifnot(tad_end >= max(summary_dt[["pos_end"]]))
      # stopifnot(tad_end == max(summary_dt[["pos_end"]]))
      stopifnot(chromo_nbr == summary_dt[["chr_num"]])
      
      sub_dt <- unique(summary_dt[,grepl("mean_p", colnames(summary_dt)) | grepl("sig_ratio", colnames(summary_dt))])
      
      # if(nrow(sub_dt) != 1) {
      #   outFile <- file.path(paste0(ds, "_",i, "_", "i_dt_debug.Rdata"))
      #   save(i_dt, file=outFile, version=2)
      #   outFile <- file.path(paste0(ds, "_",i, "_", "sub_dt_debug.Rdata"))
      #   save(sub_dt, file=outFile, version=2)
      # }
      
      # stopifnot(nrow(sub_dt) == 1)   # can be 0 e.g.                                     
      # > dt[[5122]][["summary_tab"]]
      # [1] gene                          LNCaP_prostate_cancer        
      # [3] 22Rv1_prostate_cancer         RWPE1_prostate_epithelium    
      # [5] exp_LNCaP_prostate_cancer     exp_22Rv1_prostate_cancer    
      # [7] exp_RWPE1_prostate_epithelium sig_ratio_LNCaP              
      # [9] sig_ratio_22Rv1               sig_ratio_RWPE1              
      # [11] mean_p_LNCaP                  mean_p_22Rv1                 
      # [13] mean_p_RWPE1                  pos_start                    
      # [15] pos_end                       chr_num                      
      # 
      if(nrow(sub_dt) == 0) {
        sub_dt[1,] <- NA
      }
      # return(NULL)
      stopifnot(nrow(sub_dt) == 1)   # 
      sub_dt$chromo <- chromo
      sub_dt$tad_start <- tad_start
      sub_dt$tad_end <- tad_end

      # if(ds== "GSE118514_22Rv1_40kb" & i == 5710) {
      #   outFile <- file.path(paste0(ds, "_",i, "_", "i_dt_debug.Rdata"))
      #   save(i_dt, file=outFile, version=2)
      #   outFile <- file.path(paste0(ds, "_",i, "_", "sub_dt_debug.Rdata"))
      #   save(sub_dt, file=outFile, version=2)
      # }
      # 
      # 
      # sur chacune des submats, avec hicds comme ref. compter le nombre d'interactions qui sont différentes
      all_mats <- i_dt[["hic_dc_mat"]]
      ref_ds <- as.character(all_pairs_cols[ds])
      i_ref <- which(names(all_mats) == ref_ds)
      stopifnot(length(i_ref) == 1)
      i_others <- which(names(all_mats) != ref_ds)
      stopifnot(length(i_others) == 2)
      ref_mat <- as.matrix(all_mats[[i_ref]])
      stopifnot(nrow(ref_mat) == ncol(ref_mat))

      stopifnot(nrow(ref_mat) == (tad_end-tad_start+1)/initbinsize*(initbinsize/hicdcbinsize))
      ref_values <- ref_mat[lower.tri(ref_mat, diag=TRUE)]
      stopifnot(length(ref_values) == nrow(ref_mat) * (nrow(ref_mat)+1) * 0.5)
      i_mat=2
      for(i_mat in i_others) {
        match_mat <- all_mats[[i_mat]]
        stopifnot(dim(match_mat) == dim(ref_mat))
        match_values <- match_mat[lower.tri(match_mat, diag=TRUE)]
        stopifnot(length(match_values) == nrow(match_mat) * (nrow(match_mat)+1) * 0.5)
        # matching_vals <- na.omit((ref_values <= signifThresh) == (match_values <= signifThresh))
        ### THEY ARE ALREADY IN -LOG10 !!!
        matching_vals <- na.omit((ref_values >= -log10(signifThresh)) == (match_values >= -log10(signifThresh)))

        nMatchingEntries <- length(matching_vals)
        nSameEntries <- sum(matching_vals)
        stopifnot(nSameEntries/nMatchingEntries >= 0 &  nSameEntries/nMatchingEntries <= 1)
        sub_dt[,paste0("nMatchingEntries_", ref_ds, "_", names(all_mats)[i_mat])] <- nMatchingEntries
        sub_dt[,paste0("nSameEntries_", ref_ds, "_", names(all_mats)[i_mat])] <- nSameEntries
        sub_dt[,paste0("nRatioMatchEntries_", ref_ds, "_", names(all_mats)[i_mat])] <- nSameEntries/nMatchingEntries
      } # end iterating over hicdc matrices
      
      
      sub_dt
    } # end iterating over TAD
    
    all_interactome_dt$hicds_hicdc <- ds
    # all_interactome_dt
    outFile <- file.path(outFolder, paste0(ds, "_all_interactome_dt.Rdata"))
    save(all_interactome_dt, file=outFile, version=2)
    cat(paste0("... written: ", outFile, "\n"))
    # return(NULL)
    
  } # end iterating over ds
  ### NEED TO BE STORED SEPARATELY BECAUSE OF NON MATCHING COLUMN NAMES
  # outFile <- file.path(outFolder, "all_hicdc_dt.Rdata")
  # save(all_hicdc_dt, file=outFile, version=2)
  # cat(paste0("... written: ", outFile, "\n"))
} 
# else {
#   outFile <- file.path(outFolder, "all_hicdc_dt.Rdata")
#   all_hicdc_dt <- get(load(outFile))
# }
# load("REVISION_INTERACTOME/all_hicdc_dt.Rdata")
# stop("--ok \n")


source("revision_settings.R")

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

stopifnot(length(unique(sub_final_table_DT$exprds)) == 1)

densplotcols <- c("nRatioMatchEntries", "sig_ratio", "mean_p")
densitycols <- c("signif_lab", "signif_lab2", "signif_lab3")

cmp=all_cmps[1]
all_plot_dt <- foreach(cmp = all_cmps, .combine='rbind') %dopar% {
  
  cat(paste0("> start plotting - ", cmp, "\n"))
  
  ref_hicds <- dirname(dirname(cmp))
  match_hicds <- basename(dirname(cmp))
  exprds <- basename(cmp)
  
  ref_hicds_col <- as.character(all_pairs_cols[paste0(ref_hicds)])
  match_hicds_col <- as.character(all_pairs_cols[paste0(match_hicds)])
  
  # ref_lab <- gsub(".+_(.+)_40kb", "\\1", ref_hicds)
  # match_lab <- gsub(".+_(.+)_40kb", "\\1", match_hicds)
  
  ref_lab <- ref_hicds_col
  match_lab <- match_hicds_col
  
  ds_hicdc_dt <- get(load(file.path(outFolder,paste0(ref_hicds, "_all_interactome_dt.Rdata"))))
  ds_hicdc_dt$hicds <- ds_hicdc_dt$hicds_hicdc
  
  merged_dt <- merge(sub_final_table_DT, ds_hicdc_dt, by=c("hicds", "chromo", "tad_start", "tad_end"))
  
  stopifnot(!duplicated(merged_dt$regionID))
  
  merged_dt$tad_signif <- merged_dt$adjPvalComb <= tad_signifThresh
  merged_dt$signif_lab <- ifelse(merged_dt$tad_signif, "signif.", "not signif.")
  merged_dt$signif_lab2 <- ifelse(merged_dt$meanLogFC > 0, "up", "down")
  merged_dt$signif_lab3 <- paste0(merged_dt$signif_lab, merged_dt$signif_lab2)
  merged_dt$signif_lab2[merged_dt$signif_lab == "not signif."] <- "not signif."
  
  plot_dt <- merged_dt
  

  sub_plot_dt <- plot_dt[plot_dt$hicds == ref_hicds &
                           plot_dt$exprds == exprds, ]
  stopifnot(nrow(sub_plot_dt) == nrow(plot_dt))
  stopifnot(nrow(sub_plot_dt) > 0)
  
  stopifnot(ref_hicds %in% names(all_pairs_cols))
  stopifnot(match_hicds %in% names(all_pairs_cols))

  icol ="sig_ratio"
  for(icol in densplotcols) {
    
    icol_lab <- gsub("_", "", icol)
    
    if(icol == "nRatioMatchEntries") {
      # cat(paste0(paste0(icol,"_",ref_hicds_col, "_", match_hicds_col  ), "\n"))
      # cat(paste0("..\n",colnames(sub_plot_dt), "\n"))
      sub_plot_dt$match_ref_diff <- sub_plot_dt[,paste0(icol,"_",ref_hicds_col, "_", match_hicds_col  )]
    } else {
      

      
      
      myxlab <-paste0( ref_lab, " ", icol_lab)
      myylab <- paste0(ref_lab, " TAD adj. p-val[-log10]")
      
      plotTit <- paste0(icol_lab, " ", ref_lab)
      
      mySub <- paste0(exprds, "; # signif. TADs = ", sum(sub_plot_dt$tad_signif), "/", nrow(sub_plot_dt))
      
      myy <- -log10(sub_plot_dt$adjPvalComb)
      myx <- sub_plot_dt[,paste0(icol, "_", ref_hicds_col)] 

      
      outFile <- file.path(outFolder, paste0(icol, "_", ref_hicds_col, "_vs_adjPval.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      densplot(x=myx,
               y=myy,
               xlab = myxlab,
               ylab =myylab,
               main =plotTit,
               cex.main=plotCex,
               cex.lab=plotCex,
               cex.axis=plotCex)
      mtext(side=3, text = plotTit)
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      for(plotvar in densitycols) {
        legTitle <- paste0("")
        # my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], sort(unique(sub_plot_dt[,paste0("tad_signif")])))
        my_cols <- setNames(pal_jama()(5)[c(3, 2,4, 1)], sort(unique(sub_plot_dt[,paste0(plotvar)])))
        # plotTit <- paste0(icol)
        plotTit <- paste0(icol_lab, " ", ref_lab)
        
        mySub <- paste0(names(table(sub_plot_dt[,plotvar])),"=", as.numeric(table(sub_plot_dt[,plotvar])), collapse="; ") 
        
        p3 <- ggdensity(sub_plot_dt,
                        x = paste0(paste0(icol, "_", ref_hicds_col)),
                        y = "..density..",
                        # combine = TRUE,                  # Combine the 3 plots
                        xlab = paste0(myxlab),
                        # add = "median",                  # Add median line.
                        rug = FALSE,                      # Add marginal rug
                        color = paste0(plotvar),
                        fill = paste0(plotvar),
                        # color = paste0("tad_signif"),
                        # fill = paste0("tad_signif"),
                        palette = "jco"
        ) +
          ggtitle(plotTit, subtitle = mySub)+
          scale_color_manual(values=my_cols)+
          scale_fill_manual(values=my_cols)  +
          labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
          guides(color=FALSE)+
          scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
          scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
          mytheme
        
        outFile <- file.path(outFolder, paste0(icol, "_", ref_hicds_col, "_", plotvar, "_density.", plotType))
        ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
        cat(paste0("... written: ", outFile, "\n"))
        
        
        
        
      } # end iterating over signif_lab
      
      
      
      
      
      sub_plot_dt$match_ref_diff <- sub_plot_dt[,paste0(icol, "_", match_hicds_col)] - 
        sub_plot_dt[,paste0(icol, "_", ref_hicds_col)] 
    }
    
    
    myxlab <-paste0(match_lab, " - ", ref_lab, " ", icol_lab)
    myylab <- paste0(ref_lab, " TAD adj. p-val[-log10]")
    
    plotTit <- paste0(icol_lab, " ", ref_lab, " (ref) vs. ", match_lab)
      
    mySub <- paste0(exprds, "; # signif. TADs = ", sum(sub_plot_dt$tad_signif), "/", nrow(sub_plot_dt))
    
    myy <- -log10(sub_plot_dt$adjPvalComb)
    myx <- sub_plot_dt$match_ref_diff
    
    outFile <- file.path(outFolder, paste0(icol, "_", ref_hicds_col, "_vs_", match_hicds_col, "_vs_adjPval.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(x=myx,
         y=myy,
         xlab = myxlab,
         ylab =myylab,
         main =plotTit,
         cex.main=plotCex,
         cex.lab=plotCex,
         cex.axis=plotCex)
    mtext(side=3, text = plotTit)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    
    for(plotvar in densitycols) {
      legTitle <- paste0("")
      # my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], sort(unique(sub_plot_dt[,paste0("tad_signif")])))
      my_cols <- setNames(pal_jama()(5)[c(3, 2,4, 1)], sort(unique(sub_plot_dt[,paste0(plotvar)])))
      # plotTit <- paste0(icol)
      plotTit <- paste0(icol_lab, " ", ref_lab, " (ref) vs. ", match_lab)
      
      mySub <- paste0(names(table(sub_plot_dt[,plotvar])),"=", as.numeric(table(sub_plot_dt[,plotvar])), collapse="; ") 
      
      p3 <- ggdensity(sub_plot_dt,
                      x = paste0("match_ref_diff"),
                      y = "..density..",
                      # combine = TRUE,                  # Combine the 3 plots
                      xlab = paste0(myxlab),
                      # add = "median",                  # Add median line.
                      rug = FALSE,                      # Add marginal rug
                      color = paste0(plotvar),
                      fill = paste0(plotvar),
                      # color = paste0("tad_signif"),
                      # fill = paste0("tad_signif"),
                      palette = "jco"
      ) +
        ggtitle(plotTit, subtitle = mySub)+
        scale_color_manual(values=my_cols)+
        scale_fill_manual(values=my_cols)  +
        labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
        guides(color=FALSE)+
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
        mytheme
      
      outFile <- file.path(outFolder, paste0(icol, "_", ref_hicds_col, "_vs_", match_hicds_col, "_", plotvar, "_density.", plotType))
      ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
      cat(paste0("... written: ", outFile, "\n"))
      
      
      
      
    } # end iterating over signif_lab
    
    
    sub_plot_dt[,paste0(icol, "_match_ref_diff")] <- sub_plot_dt$match_ref_diff
  } # end iterating over columns to plot
  sub_plot_dt$ref_hicds <- ref_hicds
  sub_plot_dt$match_hicds <- match_hicds
  sub_plot_dt[,c("ref_hicds","exprds","adjPvalComb", "match_hicds","tad_signif",paste0(densplotcols, "_match_ref_diff"), densitycols)]
} # end iterating over ds

outFile <- file.path(outFolder, "all_plot_dt.Rdata")
save(all_plot_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

exprds <- unique(all_plot_dt$exprds)
stopifnot(length(exprds) == 1)

icol ="sig_ratio"
for(icol in densplotcols) {
  
  plotcol <- paste0(icol, "_match_ref_diff")
  
  nCmps <- length(unique(file.path(all_plot_dt$ref_hicds, all_plot_dt$match_hicds)))
  
  icol_lab <- gsub("_", "", icol)
  
  myxlab <-paste0( "match - ref", " ", icol_lab)
  myylab <- paste0( "ref TAD adj. p-val[-log10]")
  
  plotTit <- paste0(icol_lab, " all DS (# cmps = ", nCmps, ")")
  
  mySub <- paste0(exprds, "; # signif. TADs = ", sum(all_plot_dt$tad_signif), "/", nrow(all_plot_dt))
  
  myy <- -log10(all_plot_dt$adjPvalComb)
  myx <- all_plot_dt[,plotcol]
  
  outFile <- file.path(outFolder, paste0(icol, "", "ref", "_vs_", "match", "_vs_adjPval_allDS.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(x=myx,
           y=myy,
           xlab = myxlab,
           ylab =myylab,
           main =plotTit,
           cex.main=plotCex,
           cex.lab=plotCex,
           cex.axis=plotCex)
  mtext(side=3, text = plotTit)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  for(plotvar in densitycols) {
    legTitle <- paste0("")
    # my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], sort(unique(sub_plot_dt[,paste0("tad_signif")])))
    my_cols <- setNames(pal_jama()(5)[c(3, 2,4, 1)], sort(unique(all_plot_dt[,paste0(plotvar)])))
    # plotTit <- paste0(icol)
    plotTit <- paste0(icol_lab, " allDS - ", "ref", " vs. ", "match",  " all DS (# cmps = ", nCmps, ")")
    
    mySub <- paste0(names(table(all_plot_dt[,plotvar])),"=", as.numeric(table(all_plot_dt[,plotvar])), collapse="; ") 
    
    p3 <- ggdensity(all_plot_dt,
                    # x = paste0("match_ref_diff"),
                    x = paste0(plotcol),
                    y = "..density..",
                    # combine = TRUE,                  # Combine the 3 plots
                    xlab = paste0(myxlab),
                    # add = "median",                  # Add median line.
                    rug = FALSE,                      # Add marginal rug
                    color = paste0(plotvar),
                    fill = paste0(plotvar),
                    # color = paste0("tad_signif"),
                    # fill = paste0("tad_signif"),
                    palette = "jco"
    ) +
      ggtitle(plotTit, subtitle = mySub)+
      scale_color_manual(values=my_cols)+
      scale_fill_manual(values=my_cols)  +
      labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
      guides(color=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
      mytheme
    
    outFile <- file.path(outFolder, paste0(icol, "ref_vs_match", "_", plotvar, "_density_allDS.", plotType))
    ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
  } # end iterating over signif_lab
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
