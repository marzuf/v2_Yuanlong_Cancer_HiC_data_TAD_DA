startTime <- Sys.time()

require(ggplot2)
require(ggpubr)

options(scipen=100)

source("revision_settings.R")
# to get all_pairs

source("../MANUSCRIPT_FIGURES/full_dataset_names.R")

source("revision_sim_metrics.R")

# Rscript revision_cellLines_sim.R

plotType <- "png"
myHeightGG <- 5
myWidthGG <- 6
myHeight <- 400
myWidth <- 400

buildTable <- FALSE

outFolder <- file.path("REVISION_CELLLINES_SIM")

require(foreach)
require(doMC)
registerDoMC(20)

all_hicds <- names(hicds_names)

runFold <- "."

bin_size <- 40*10^3
coverTADmatchRatio <-  0.8
boundariesJI_tolRad <- 2*bin_size


onlyPipTADs <- FALSE  # only implemented for normal vs. tumor !!! not working luad and lusc!!!

if(onlyPipTADs) stop("--not correctly implemented")

if(onlyPipTADs) outFolder <- file.path("REVISION_CELLLINES_SIM_ONLYPIPTADS")
dir.create(outFolder, recursive = TRUE)

# all_hicds <- all_hicds[1:3]
all_ds_pairs <- combn(x=all_hicds, m=2)
stopifnot(nrow(all_ds_pairs) == 2)
stopifnot(ncol(all_ds_pairs) == length(all_hicds) * (length(all_hicds) -1) * 0.5)


hicds = all_hicds[1]

if(buildTable) {
  all_tads <- foreach(hicds = all_hicds) %dopar% {
    tad_file <- file.path(runFold, hicds, "genes2tad", "all_assigned_regions.txt" )
    genome_dt <- read.delim(tad_file, sep="\t", header=FALSE, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
    tad_dt <- genome_dt[grepl("_TAD", genome_dt$region),]
    
    if(onlyPipTADs) {
      hicds_exprds <- list.files(file.path("PIPELINE", "OUTPUT_FOLDER", hicds))
      hicds_exprds <- hicds_exprds[grepl("_norm_", hicds_exprds)]
      if(length(hicds_exprds) == 0) return(NULL)
      stopifnot(length(hicds_exprds) == 1)
      pip_tads <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, hicds_exprds, "0_prepGeneData", "pipeline_regionList.Rdata")))
      stopifnot(pip_tads %in% tad_dt$region)
      tad_dt <- tad_dt[tad_dt$region %in% pip_tads,]
    }
    
    stopifnot(nrow(tad_dt) > 0)
    tad_dt
  }
  names(all_tads) <- all_hicds
  outFile <- file.path(outFolder, "all_tads.Rdata")
  save(all_tads, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  all_pairs_dt <- foreach(i = 1:ncol(all_ds_pairs), .combine='rbind') %dopar% {
    
    
    hicds1 <- all_ds_pairs[1,i]
    hicds2 <- all_ds_pairs[2,i]
    
    stopifnot(hicds1 %in% names(all_tads))
    stopifnot(hicds2 %in% names(all_tads))
    
    
    hicds1_tads <- all_tads[[paste0(hicds1)]]
    hicds2_tads <- all_tads[[paste0(hicds2)]]
    
    stopifnot(grepl("_TAD", hicds1_tads$region))
    hicds1_tads$region <- NULL
    stopifnot(grepl("_TAD", hicds2_tads$region))
    hicds2_tads$region <- NULL
    
    all_chromo <- unique(c(hicds1_tads$chromo, hicds2_tads$chromo))
    stopifnot(all_chromo %in% hicds1_tads$chromo)
    stopifnot(all_chromo %in% hicds2_tads$chromo)
    
    # all_chromo=all_chromo[1]
    
    chr_dt <- foreach(chr = all_chromo, .combine='rbind') %dopar% {
      
      cat(paste0("... start:\t", hicds1, " vs. ", hicds2, " - ", chr, "\n"))
      
      hicds1_chr_tads <- hicds1_tads[hicds1_tads$chromo == chr,]
      stopifnot(nrow(hicds1_chr_tads) > 0)
      
      hicds2_chr_tads <- hicds2_tads[hicds2_tads$chromo == chr,]
      stopifnot(nrow(hicds2_chr_tads) > 0)
      
      chrsize <- max(c(hicds1_chr_tads$end, hicds2_chr_tads$end))
      stopifnot(is.numeric(chrsize))
      
      cat(paste0("chrsize= ", chrsize, "\n"))
      

      ds1_ds2_moc <- get_MoC(hicds1_chr_tads, hicds2_chr_tads, chrSize=chrsize)
      ds1_ds2_binJI <- get_bin_JaccardIndex(hicds1_chr_tads, hicds2_chr_tads, binSize=bin_size)
      ds1_ds2_boundJI <- get_boundaries_JaccardIndex(hicds1_chr_tads, hicds2_chr_tads, tolRad=boundariesJI_tolRad,matchFor="set1" )
      ds1_ds2_matchingTADs <- get_ratioMatchingTADs(hicds1_chr_tads, hicds2_chr_tads, coverMatchRatioThresh=coverTADmatchRatio, matchFor="set1" )
      ds1_ds2_VI <- get_variationInformation(hicds1_chr_tads, hicds2_chr_tads)
      
      ds2_ds1_boundJI <- get_boundaries_JaccardIndex(hicds2_chr_tads, hicds1_chr_tads, 
                                                     tolRad=boundariesJI_tolRad,matchFor="set1" )
      ds2_ds1_matchingTADs <- get_ratioMatchingTADs(hicds2_chr_tads, hicds1_chr_tads, 
                                                    coverMatchRatioThresh=coverTADmatchRatio, matchFor="set1" )
      ds2_ds1_VI <- get_variationInformation(hicds2_chr_tads, hicds1_chr_tads)
      
      stopifnot(ds1_ds2_VI == ds2_ds1_VI)
      
      # rbind(
      data.frame(
        ds1 = hicds1,
        ds2 = hicds2, 
        chromo = chr,
        moc =ds1_ds2_moc,
        binJI =ds1_ds2_binJI,
        VI = ds1_ds2_VI,
        ds1_ds2_boundJI =ds1_ds2_boundJI,
        ds1_ds2_matchingTADs = ds1_ds2_matchingTADs,
        ds2_ds1_boundJI =ds2_ds1_boundJI,
        ds2_ds1_matchingTADs = ds2_ds1_matchingTADs,
        # ds2_ds1_VI = ds2_ds1_VI,
        stringsAsFactors = FALSE
      )
      # ,data.frame(
      #   ds1 = hicds2,
      #   ds2 = hicds1, 
      #   moc =ds1_ds2_moc,  # symmetric
      #   binJI =ds1_ds2_binJI, # symmetric
      #   ds1_ds2_boundJI =ds2_ds1_boundJI,
      #   ds1_ds2_matchingTADs = ds2_ds1_matchingTADs,
      #   ds1_ds2_VI = ds2_ds1_VI,
      #   stringsAsFactors = FALSE
      # ))
    } # end iterating over chromo
    chr_dt
  } # end iterating over pair
  
  outFile <- file.path(outFolder, "all_pairs_dt.Rdata")
  save(all_pairs_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
}else {
  outFile <- file.path(outFolder, "all_pairs_dt.Rdata")
  all_pairs_dt <- get(load(outFile))
  cat(paste0("nrow(all_pairds_dt)=", nrow(all_pairs_dt), "\n"))
}


stopifnot(nrow(all_pairs_dt) == 22 * ncol(all_ds_pairs))

all_plot_cols <- colnames(all_pairs_dt)[!colnames(all_pairs_dt) %in% c("ds1", "ds2", "chromo")]
toplot=all_plot_cols[1]


# load("REVISION_CELLLINES_SIM/all_pairs_dt.Rdata")
stopifnot(!duplicated(file.path(all_pairs_dt$ds1, all_pairs_dt$ds2, all_pairs_dt$chromo)))
all_pairs_dt$matchPair <- file.path(all_pairs_dt$ds1, all_pairs_dt$ds2)

all_pairs_init <- all_pairs
all_pairs <- file.path( dirname(dirname(all_pairs_init)),   basename(dirname(all_pairs_init)))
rev_pairs <- file.path(  basename(dirname(all_pairs_init)), dirname(dirname(all_pairs_init)))
stopifnot(rev_pairs %in% all_pairs_dt$matchPair |
            all_pairs %in% all_pairs_dt$matchPair )



all_pairs_dt$tissue1 <- all_tissues[all_pairs_dt$ds1]
stopifnot(!is.na(all_pairs_dt$tissue1))
all_pairs_dt$tissue2 <- all_tissues[all_pairs_dt$ds2]
stopifnot(!is.na(all_pairs_dt$tissue2))

all_pairs_dt$pairType <- ifelse(all_pairs_dt$tissue1 == all_pairs_dt$tissue2,  "same_tissue", "diff_tissue")
stopifnot(which(all_pairs_dt$matchPair %in% c(rev_pairs, all_pairs)) %in% which(all_pairs_dt$pairType=="same_tissue"))
all_pairs_dt$pairType[all_pairs_dt$matchPair %in% c(rev_pairs, all_pairs)] <-  "norm_tumor_pair"

ncmps <- length(unique(file.path(all_pairs_dt$ds1,all_pairs_dt$ds2 )))

myset <- paste0("(TADcoverMatch=",coverTADmatchRatio, "; bdTolRad=", boundariesJI_tolRad, ")") 


all_pairs_dt$pairType <- factor(all_pairs_dt$pairType, levels=c("diff_tissue", "same_tissue", "norm_tumor_pair"))
stopifnot(!is.na(all_pairs_dt$pairType))

for(toplot in all_plot_cols) {
  
  mysub2 <- paste0(names(table(all_pairs_dt$pairType)), "=", as.numeric(table(all_pairs_dt$pairType)), 
                   collapse=";")
  
  plotTit <- paste0("ds1 vs. ds2 - ", toplot, " (# cmps = ", ncmps,")")
  # mySub <- paste0( "1 dot/chromo - onlyPipTADs=", as.character(onlyPipTADs),myset)
  mySub <- paste0( "1 dot/chromo ",myset)
  legTitle <- ""
  
  p3 <-  ggboxplot(all_pairs_dt,
                   y = paste0(toplot),
                   outlier.shape=NA,
                   add="jitter",
                   rug = FALSE,  
                   ylab = paste0(toplot),
                   xlab ="",
                   x="pairType",                   # fill = paste0("tad_signif"),
                   color = "pairType",
                   # fill = "signif_lab2",
                   palette = "aaas") +
    ggtitle(plotTit, subtitle = mySub) +
    labs(color=paste0(legTitle),fill=paste0(legTitle)) +
    guides(color=FALSE, fill=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  
  outFile <- file.path(outFolder, paste0("ds1_vs_d2_", toplot, "_boxplot.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
}


### Look at agreement among metrics

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

# for(var1 in all_plot_cols[1:(length(all_plot_cols)-1)]) {
#   for(var2 in all_plot_cols[2:(length(all_plot_cols))]) {
for(i_var1 in c(1:(length(all_plot_cols)-1))) {
  var1 <- all_plot_cols[i_var1]
  for(i_var2 in c((1+i_var1):length(all_plot_cols))) {
    var2 <- all_plot_cols[i_var2]
    outFile <- file.path(outFolder, paste0("check_", var2, "_vs_", var1, ".", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x=all_pairs_dt[, var1],
      y=all_pairs_dt[, var2],
      xlab=var1,
      ylab=var2,
      main=paste0(var2, " vs. ", var1),
      cex.lab=1.2,
      cex.main=1.2,
      cex.axis=1.2
    )
    addCorr(x=all_pairs_dt[, var1],
            y=all_pairs_dt[, var2],
            legPos="topleft",
            bty="n")
    mtext(side=3, text=paste0("# cmps=", ncmps, "; # points=", nrow(all_pairs_dt)))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }
}

##############################################3
cat("***Done\n")
cat(paste0(startTime , " - ", Sys.time(), "\n"))