startTime <- Sys.time()
require(ggplot2)
require(ggsci)
require(ggpubr)
source("revision_settings.R")

source("../MANUSCRIPT_FIGURES/full_dataset_names.R")
all_hicds <- names(hicds_names)

require(GenomicRanges)

binSize <- 40* 10^3
match_tolRad <- 2*binSize
match_coverRatio <- 0.8

source("revision_sim_metrics.R")

plotType <- "png"
myHeightGG <- 5
myWidthGG <- 6


buildTable <- FALSE

# Rscript revision_cellLines_sim_signif_allPairs.R

outFolder <- file.path("REVISION_CELLLINES_SIM_SIGNIF_ALLPAIRS")
dir.create(outFolder, recursive=TRUE )

# all_hicds=all_hicds[1:2]

pvalthresh <- 0.01
runFolder <- "."
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
resultData$region_id <- file.path(resultData$dataset, resultData$region)
resultData$signif <- resultData$adjPvalComb <= pvalthresh
signif_tads <- resultData$region_id[resultData$signif ]

final_normTum_dt <- resultData[grepl("_norm_", resultData$exprds),]
stopifnot(nrow(final_normTum_dt) > 0)
final_normTum_dt$region_id_hicds <- file.path(final_normTum_dt$hicds, final_normTum_dt$region)
normtum_signif_tads <- final_normTum_dt$region_id_hicds[final_normTum_dt$signif ]

require(foreach)
require(doMC)
registerDoMC(50)

onlyPipTADs <- FALSE
if(onlyPipTADs) stop("--not correctly implemented")

runFold <- "."

myset <- paste0("(TADcoverMatch=",match_coverRatio, "; bdTolRad=", match_tolRad, ")") 


if(buildTable){
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
  
  all_ds_pairs <- combn(x=all_hicds, m=2)
  stopifnot(nrow(all_ds_pairs) == 2)
  stopifnot(ncol(all_ds_pairs) == length(all_hicds) * (length(all_hicds) -1) * 0.5)
  
  all_pairs_dt <- foreach(i = 1:ncol(all_ds_pairs), .combine='rbind') %dopar% {
    
    
    hicds1 <- all_ds_pairs[1,i]
    hicds2 <- all_ds_pairs[2,i]
    
    stopifnot(hicds1 %in% names(all_tads))
    stopifnot(hicds2 %in% names(all_tads))
    
    
    hicds1_tads <- all_tads[[paste0(hicds1)]]
    hicds2_tads <- all_tads[[paste0(hicds2)]]
    
    
    stopifnot(grepl("_TAD", hicds1_tads$region))
    # hicds1_tads$region <- NULL
    stopifnot(grepl("_TAD", hicds2_tads$region))
    # hicds2_tads$region <- NULL
    
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
      
      
      dt_set1_boundaryMatch <- get_set1_boundaryMatch(set1DT=hicds1_chr_tads, 
                                                      set2DT=hicds2_chr_tads, 
                                                      tolRad=match_tolRad)
      
      head(dt_set1_boundaryMatch)
      
      # because they are converted
      dt_set1_boundaryMatch$start <- dt_set1_boundaryMatch$start + 1
      
      cat(paste0("nrow dt_set1_boundaryMatch = ", nrow(dt_set1_boundaryMatch), "\n"))
      
      dt_set1_tadCoverMatch <- get_set1_tadCoverMatch(set1DT=hicds1_chr_tads, 
                                                      set2DT=hicds2_chr_tads,
                                                      coverMatchRatioThresh= match_coverRatio) 
      
      head(dt_set1_tadCoverMatch)
      
      cat(paste0("nrow dt_set1_tadCoverMatch = ", nrow(dt_set1_tadCoverMatch), "\n"))
      
      # save(dt_set1_boundaryMatch, file="dt_set1_boundaryMatch", version=2)
      # save(dt_set1_tadCoverMatch, file="dt_set1_tadCoverMatch", version=2)
      match_dt1 <- merge(dt_set1_boundaryMatch, dt_set1_tadCoverMatch, by=c("chromo", "start", "end", "region"))
      # match_dt1 <- full_join(dt_set1_boundaryMatch, dt_set1_tadCoverMatch, by=c("chromo", "start", "end", "region"))
      
      cat(paste0("nrow match_dt1 = ", nrow(match_dt1), "\n"))
      
      
      stopifnot(nrow(match_dt1) == nrow(dt_set1_tadCoverMatch))
      stopifnot(nrow(match_dt1) == nrow(dt_set1_boundaryMatch))
      stopifnot(nrow(match_dt1) == nrow(hicds1_chr_tads))
      
      match_dt1$ref_ds <- hicds1
      match_dt1$match_ds <- hicds2
      
      ################## do the same but ds2 vs ds1
      
      
      
      dt_set2_boundaryMatch <- get_set1_boundaryMatch(set1DT=hicds2_chr_tads, 
                                                      set2DT=hicds1_chr_tads, 
                                                      tolRad=match_tolRad)
      
      head(dt_set2_boundaryMatch)
      
      # because they are converted
      dt_set2_boundaryMatch$start <- dt_set2_boundaryMatch$start + 1
      
      cat(paste0("nrow dt_set2_boundaryMatch = ", nrow(dt_set2_boundaryMatch), "\n"))
      
      dt_set2_tadCoverMatch <- get_set1_tadCoverMatch(set1DT=hicds2_chr_tads, 
                                                      set2DT=hicds1_chr_tads,
                                                      coverMatchRatioThresh= match_coverRatio) 
      
      head(dt_set2_tadCoverMatch)
      
      cat(paste0("nrow dt_set2_tadCoverMatch = ", nrow(dt_set2_tadCoverMatch), "\n"))
      
      # save(dt_set2_boundaryMatch, file="dt_set2_boundaryMatch", version=2)
      # save(dt_set2_tadCoverMatch, file="dt_set2_tadCoverMatch", version=2)
      match_dt2 <- merge(dt_set2_boundaryMatch, dt_set2_tadCoverMatch, by=c("chromo", "start", "end", "region"))
      # match_dt2 <- full_join(dt_set2_boundaryMatch, dt_set2_tadCoverMatch, by=c("chromo", "start", "end", "region"))
      
      cat(paste0("nrow match_dt2 = ", nrow(match_dt2), "\n"))
      
      
      stopifnot(nrow(match_dt2) == nrow(dt_set2_tadCoverMatch))
      stopifnot(nrow(match_dt2) == nrow(dt_set2_boundaryMatch))
      stopifnot(nrow(match_dt2) == nrow(hicds2_chr_tads))
      
      match_dt2$ref_ds <- hicds2
      match_dt2$match_ds <- hicds1
      
      
      match_dt <- rbind(match_dt1, match_dt2)
      match_dt
      
      
    } # end iterating over chromo
    chr_dt
  } # end iterating over pair
  
  outFile <- file.path(outFolder, "all_pairs_dt.Rdata")
  save(all_pairs_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  cat(paste0("... loading data\n"))
  outFile <- file.path(outFolder, "all_pairs_dt.Rdata")
  all_pairs_dt <- get(load(outFile))
  outFile <- file.path(outFolder, "all_tads.Rdata")
  all_tads <- get(load(outFile))
}

all_plot_vars <- colnames(all_pairs_dt)[!colnames(all_pairs_dt)%in% c("ref_ds", "match_ds", "chromo", "start", "end", "region")]

stopifnot(!duplicated(file.path(all_pairs_dt$chromo, all_pairs_dt$start, all_pairs_dt$end, 
                                all_pairs_dt$region, all_pairs_dt$ref_ds, all_pairs_dt$match_ds, all_pairs_dt$exprds)))
stopifnot(which(all_pairs_dt$end1_vs_end2_match) %in% which(all_pairs_dt$end1_vs_all2_match))   
stopifnot(all(which(all_pairs_dt$start1_vs_start2_match) %in% which(all_pairs_dt$start1_vs_all2_match)))

all_pairs_dt$matchPair <- file.path(all_pairs_dt$ref_ds, all_pairs_dt$match_ds)

all_pairs_init <- all_pairs
all_pairs <- file.path( dirname(dirname(all_pairs_init)),   basename(dirname(all_pairs_init)))
rev_pairs <- file.path(  basename(dirname(all_pairs_init)), dirname(dirname(all_pairs_init)))
stopifnot(rev_pairs %in% all_pairs_dt$matchPair |
            all_pairs %in% all_pairs_dt$matchPair )

all_pairs_dt$ref_tissue <- all_tissues[all_pairs_dt$ref_ds]
stopifnot(!is.na(all_pairs_dt$ref_tissue))
all_pairs_dt$match_tissue <- all_tissues[all_pairs_dt$match_ds]
stopifnot(!is.na(all_pairs_dt$match_tissue))

all_pairs_dt$pairType <- ifelse(all_pairs_dt$ref_tissue == all_pairs_dt$match_tissue,  "same_tissue", "diff_tissue")
stopifnot(which(all_pairs_dt$matchPair %in% c(rev_pairs, all_pairs)) %in% which(all_pairs_dt$pairType=="same_tissue"))
all_pairs_dt$pairType[all_pairs_dt$matchPair %in% c(rev_pairs, all_pairs)] <-  "norm_tumor_pair"

ncmps <- length(unique(file.path(all_pairs_dt$ref_ds,all_pairs_dt$match_ds )))

toplot <- "start1_vs_start2_match"


for(toplot in all_plot_vars) {

  cat(paste0("... ", toplot, " - agg data \n"))

  agg_dt <- aggregate(as.formula(paste0(toplot, " ~ ref_ds + match_ds + pairType")),
                      FUN = function(x) sum(x)/length(x), data=all_pairs_dt)

  mysub2 <- paste0(names(table(agg_dt$pairType)), "=", as.numeric(table(agg_dt$pairType)),
                   collapse=";")

  # plotTit <- paste0("ds1 vs. ds2 - ", toplot)
  plotTit <- paste0("ds1 vs. ds2 - ", toplot, " (# cmps = ", ncmps,")")
  # mySub <- paste0("# cmps = ", ncmps, "1 dot/comparison - onlyPipTADs=", as.character(onlyPipTADs) )
  mySub <- paste0( "1 dot/comparison ",myset)
  
  legTitle <- ""

  p3 <-  ggboxplot(agg_dt,
                   y = paste0(toplot),
                   outlier.shape=NA,
                   add="jitter",
                   rug = FALSE,
                   ylab = paste0("ratio ", toplot, "TADs"),
                   xlab ="",
                   x="pairType",                   # fill = paste0("tad_signif"),
                   color = "pairType",
                   # fill = "signif_lab2",
                   palette = "aaas") +
    ggtitle(plotTit, subtitle = mySub) +
    labs(color=paste0(legTitle),fill=paste0(legTitle)) +
    guides(color=FALSE, fill=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

  outFile <- file.path(outFolder, paste0("ref_vs_match_", toplot, "_", "ratio_boxplot.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))

}

sub_dt <- all_pairs_dt[all_pairs_dt$pairType == "norm_tumor_pair",]
sub_dt$region_id_hicds <- file.path(sub_dt$ref_ds, sub_dt$region)
sub_dt$signif_region <- sub_dt$region_id_hicds %in% normtum_signif_tads

stopifnot(nrow(sub_dt) > 0)

save(sub_dt, file="sub_dt.Rdata", version=2)

for(toplot in all_plot_vars) {
  
  cat(paste0("...SUB ", toplot, " - agg data sub\n"))
  
  agg_dt <- aggregate(as.formula(paste0(toplot, " ~ ref_ds + match_ds + pairType + signif_region")), 
                      FUN = function(x) sum(x)/length(x), data=sub_dt)
  
  mysub2 <- paste0(names(table(file.path(agg_dt$pairType,agg_dt$signif_region ))), "=", 
                   as.numeric(table(file.path(agg_dt$pairType, agg_dt$signif_region))), 
                   collapse=";")
  
  plotTit <- paste0("ds1 vs. ds2 - ", toplot, " (# cmps = ", ncmps,")")
  # mySub <- paste0( "1 dot/comparison - onlyPipTADs=", as.character(onlyPipTADs) )
  mySub <- paste0( "1 dot/comparison ",myset)
  
  legTitle <- "signif. "
  
  p3 <-  ggboxplot(agg_dt,
                   y = paste0(toplot),
                   outlier.shape=NA,
                   add="jitter",
                   rug = FALSE,  
                   ylab = paste0("ratio ", toplot, "TADs"),
                   xlab ="",
                   x="pairType",                   # fill = paste0("tad_signif"),
                   color = "signif_region",
                   # fill = "signif_lab2",
                   palette = "aaas") +
    ggtitle(plotTit, subtitle = mySub) +
    labs(color=paste0(legTitle),fill=paste0(legTitle)) +
    # guides(color=TRUE)+FALSE, fill=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  
  outFile <- file.path(outFolder, paste0("ref_vs_match_", toplot, "_", "ratio_bySignif_boxplot.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
}
# get the ratio

cat("***Done\n")
cat(paste0(startTime , " - ", Sys.time(), "\n"))