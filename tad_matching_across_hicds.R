options(scipen=100)

# Rscript tad_matching_across_hicds.R

script_name <- "tad_matching_across_hicds.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
registerDoMC(40)

require(data.table)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

buildTable <- TRUE

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 25
myHeightGG <- 25


hm.palette <- colorRampPalette(brewer.pal(9, 'YlOrRd'), space='Lab')


script0_name <- "0_prepGeneData"

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)


mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]

file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- "TAD_MATCHING_ACROSS_HICDS_NOPERMUT"
dir.create(outFolder, recursive = TRUE)

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

# => best TAD matching
# in # of genes
# in bp

####################################################################################################################################### >>> prepare the data


if(buildTable) {
  cat("... start preparing data before matching \n")

    hicds = all_hicds[1]
  all_data_list <- foreach(hicds = all_hicds) %dopar% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(hicds_file))
    tadpos_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
    tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]
    stopifnot(nrow(tadpos_dt) > 0 )
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    
    stopifnot(g2t_dt$entrezID %in% gff_dt$entrezID)
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      gene_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(gene_file))  
      
      region_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
      stopifnot(file.exists(region_file))  
      
      geneList <- get(load(gene_file))
      regionList <- get(load(region_file))
      
      stopifnot(geneList %in% g2t_dt$entrezID)
      ref_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      stopifnot(setequal(ref_g2t_dt$region, regionList))
      
      stopifnot(regionList %in% tadpos_dt$region)
      ref_tadpos_dt <- tadpos_dt[tadpos_dt$region %in% regionList,]
      stopifnot(setequal(ref_g2t_dt$region, ref_tadpos_dt$region))
      
      list(
        dataset_geneList = geneList,
        dataset_g2t_dt =ref_g2t_dt ,
        dataset_tadpos_dt = ref_tadpos_dt
      )
    }
    names(exprds_list) <- file.path(hicds, all_exprds[[paste0(hicds)]])
    exprds_list
  }
  all_data_list <- unlist(all_data_list, recursive=FALSE)
  stopifnot(length(all_data_list) == length(all_datasets))
  
  ####################################################################################################################################### >>> collect # of TADs and TAD size
  
  
  cat("... start matching\n")
  
  ref_dataset = all_datasets[1]
  # all_datasets=all_datasets[1:2]
  all_matching_dt <- foreach(ref_dataset = all_datasets, .combine='rbind') %do% {
    
    stopifnot(ref_dataset %in% names(all_data_list))
    ref_geneList <- all_data_list[[paste0(ref_dataset)]][["dataset_geneList"]]
    ref_g2t_dt <- all_data_list[[paste0(ref_dataset)]][["dataset_g2t_dt"]]
    ref_tadpos_dt <- all_data_list[[paste0(ref_dataset)]][["dataset_tadpos_dt"]]
    
    stopifnot(setequal(ref_g2t_dt$region, ref_tadpos_dt$region))
    ref_tads <- as.character(ref_tadpos_dt$region)
    
    size_dt <- ref_tadpos_dt
    size_dt$ref_totBp <- size_dt$end - size_dt$start + 1
    
    ngenes_dt <- aggregate(entrezID~region, FUN=length, data=ref_g2t_dt)
    colnames(ngenes_dt)[colnames(ngenes_dt)=="entrezID"] <- "ref_nGenes"
    
    ref_tadpos_dt_info <- merge(size_dt[,c("region", "chromo", "ref_totBp", "start", "end")], ngenes_dt, by="region")
    stopifnot(is.numeric(ref_tadpos_dt_info$start))
    stopifnot(is.numeric(ref_tadpos_dt_info$end))
    ref_tadpos_dt_info <- ref_tadpos_dt_info[order(ref_tadpos_dt_info$chromo, ref_tadpos_dt_info$start, ref_tadpos_dt_info$end),]
    colnames(ref_tadpos_dt_info)[colnames(ref_tadpos_dt_info) == "region"] <- "refID"
    
    
    query_datasets <- all_datasets[all_datasets != ref_dataset]
    stopifnot(length(query_datasets)+1 == length(all_datasets))
    query_dataset=query_datasets[2]
    
    ref_dt <- foreach(query_dataset = query_datasets, .combine='rbind') %do% {
      cat("> Start ", ref_dataset, " vs. ", query_dataset, "\n")  
      
      query_geneList <- all_data_list[[paste0(query_dataset)]][["dataset_geneList"]]
      query_g2t_dt <- all_data_list[[paste0(query_dataset)]][["dataset_g2t_dt"]]
      query_tadpos_dt <- all_data_list[[paste0(query_dataset)]][["dataset_tadpos_dt"]]
      
      ### PREPARE THE BP OVERLAP
      cat("... preparing bp matching overlap\n")
      ref_GR <- GRanges(seqnames=ref_tadpos_dt$chromo, ranges=IRanges(start=ref_tadpos_dt$start, end=ref_tadpos_dt$end, names=ref_tadpos_dt$region))
      query_GR <- GRanges(seqnames=query_tadpos_dt$chromo, ranges=IRanges(start=query_tadpos_dt$start, end=query_tadpos_dt$end, names=query_tadpos_dt$region))
      
      # determine which features from the query overlap which features in the subject
      overlap_GR <- findOverlaps(query=query_GR, subject=ref_GR)
      
      
      IDoverlap_hits_all <- findOverlaps(query=query_GR,
                                         subject=query_GR)
      IDoverlap_hits <- IDoverlap_hits_all[queryHits(IDoverlap_hits_all) != subjectHits(IDoverlap_hits_all)] 
      
      refID <- names(ref_GR[subjectHits(overlap_GR)])
      queryID <- names(query_GR[queryHits(overlap_GR)])
      
      overlapDT <- data.frame(
        refID = refID,
        queryID = queryID,
        overlapBp = width(pintersect(ref_GR[refID], query_GR[queryID])),
        stringsAsFactors = FALSE)
      stopifnot(overlapDT$refID %in% ref_tadpos_dt$region)
      stopifnot(overlapDT$queryID %in% query_tadpos_dt$region)
      
      # take the max overlap
      maxBp_overlapDT <- do.call(rbind, by(overlapDT, overlapDT$refID, function(subDT){  # WARNING: if exactly same bp overlap -> will take the 1st one...
        subDT[which.max(subDT$overlapBp),]
      }))
      stopifnot(rownames(maxBp_overlapDT) == maxBp_overlapDT$refID)
      rownames(maxBp_overlapDT) <- NULL
      stopifnot(setequal(overlapDT$refID, maxBp_overlapDT$refID))
      
      
      ### PREPARE THE GENE OVERLAP
      cat("... preparing gene matching overlap\n")
      ref_tad_genes <- split(ref_g2t_dt$entrezID, ref_g2t_dt$region)
      query_tad_genes <- split(query_g2t_dt$entrezID, query_g2t_dt$region)
      
      # all_gene_overlap_list <- lapply(ref_tad_genes, function(tad_genes) { # => FASTER TO PARALLELIZE THIS PART ???
      #   query_matches <- sapply(query_tad_genes, function(query_genes) {
      #     nmatch <- sum(unlist(tad_genes) %in% query_genes)
      #     ifelse(nmatch == 0, NA, nmatch)
      #   })
      #   # if there is no match -> vector of length 1: [1] "-Inf"
      #   # if there is match -> vector of length 2: [1] "chr10_TAD108" "5"
      #   nMatch <- max(query_matches, na.rm=TRUE)
      #   if(nMatch == "-Inf") {
      #     intersect_symbols <- NULL
      #   }else {
      #     matching_query_genes <- query_tad_genes[which.max(query_matches)]
      #     intersect_genes <- intersect(unlist(tad_genes), unlist(matching_query_genes))
      #     stopifnot(length(intersect_genes) == nMatch)
      #     stopifnot(intersect_genes %in% gff_dt$entrezID)
      #     intersect_symbols <- sort(gff_dt$symbol[gff_dt$entrezID %in% intersect_genes])
      #   }
      #   c(names(query_tad_genes)[which.max(query_matches)], nMatch, paste0(intersect_symbols, collapse=","))
      # })
      
      all_gene_overlap_list <- foreach(tad_genes=ref_tad_genes) %dopar% {
        query_matches <- sapply(query_tad_genes, function(query_genes) {
          nmatch <- sum(unlist(tad_genes) %in% query_genes)
          ifelse(nmatch == 0, NA, nmatch)
        })
        # if there is no match -> vector of length 1: [1] "-Inf"
        # if there is match -> vector of length 2: [1] "chr10_TAD108" "5"
        nMatch <- max(query_matches, na.rm=TRUE)
        if(nMatch == "-Inf") {
          intersect_symbols <- NULL
        }else {
          matching_query_genes <- query_tad_genes[which.max(query_matches)]
          intersect_genes <- intersect(unlist(tad_genes), unlist(matching_query_genes))
          stopifnot(length(intersect_genes) == nMatch)
          stopifnot(intersect_genes %in% gff_dt$entrezID)
          intersect_symbols <- sort(gff_dt$symbol[gff_dt$entrezID %in% intersect_genes])
        }
        c(names(query_tad_genes)[which.max(query_matches)], nMatch, paste0(intersect_symbols, collapse=","))
      }
      names(all_gene_overlap_list) <- names(ref_tad_genes)
      
      gene_overlap_list <- Filter(function(x) x[1] != "-Inf", all_gene_overlap_list)
      stopifnot(lengths(gene_overlap_list) == 3)
      
      maxGenes_overlapDT <- data.frame(do.call(rbind, gene_overlap_list))
      maxGenes_overlapDT$refID <- rownames(maxGenes_overlapDT)
      colnames(maxGenes_overlapDT)[1] <- "queryID"
      colnames(maxGenes_overlapDT)[2] <- "nOverlapGenes"
      colnames(maxGenes_overlapDT)[3] <- "overlapGenes"
      maxGenes_overlapDT$nOverlapGenes <- as.numeric(as.character(maxGenes_overlapDT$nOverlapGenes))
      stopifnot(maxGenes_overlapDT$refID %in% ref_g2t_dt$region)
      stopifnot(maxGenes_overlapDT$queryID %in% query_g2t_dt$region)
      stopifnot(!is.na(maxGenes_overlapDT$nOverlapGenes))  
      
      colnames(maxBp_overlapDT)[colnames(maxBp_overlapDT) == "queryID"] <- "matchingID_maxOverlapBp"
      colnames(maxGenes_overlapDT)[colnames(maxGenes_overlapDT) == "queryID"] <-  "matchingID_maxOverlapGenes"
      
      
      stopifnot(maxBp_overlapDT$refID %in% ref_tadpos_dt_info$refID)
      stopifnot(maxGenes_overlapDT$refID %in% ref_tadpos_dt_info$refID)
      
      
      tmp_dt <- merge(maxBp_overlapDT, maxGenes_overlapDT, by ="refID", all = TRUE)
      
      all_match_dt <- merge(tmp_dt, ref_tadpos_dt_info[,c("refID", "ref_totBp", "ref_nGenes"), drop=FALSE], by="refID", all=TRUE)
      stopifnot(nrow(all_match_dt) == nrow(ref_tadpos_dt_info))
      
      first_cols <- colnames(all_match_dt)[grepl("^ref",colnames(all_match_dt))]
      last_cols <- colnames(all_match_dt)[!grepl("^ref",colnames(all_match_dt))]
      all_match_dt <- all_match_dt[,c(first_cols, last_cols)]
      all_cols <- colnames(all_match_dt)
      
      all_match_dt$ref_dataset <- ref_dataset
      all_match_dt$matching_dataset <- query_dataset
      
      rownames(all_match_dt) <- NULL
      
      all_match_dt[,c("ref_dataset", "matching_dataset", all_cols)]
    } # end iterating over query_datasets
    ref_dt
  } # end iterating over ref_datasets
  
  outFile <- file.path(outFolder, "all_matching_dt.Rdata")
  save(all_matching_dt, file=outFile, version=2)  
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_matching_dt.Rdata")
  cat("... load data\n")
  all_matching_dt <- get(load(outFile))
}

cat(paste0(length(unique(all_matching_dt$ref_dataset))), "\n")

colnames(all_matching_dt)[colnames(all_matching_dt) == "overlapBP"] <- "overlapBp"

# need to select the best for each pair !
all_matching_dt$overlapBp_ratio <- all_matching_dt$overlapBp/all_matching_dt$ref_totBp
all_matching_dt$nOverlapGenes_ratio <- all_matching_dt$nOverlapGenes/all_matching_dt$ref_nGenes

nTotDS <- length(unique(all_matching_dt$ref_dataset))
nTotDS2 <- length(unique(all_matching_dt$matching_dataset))

stopifnot(nTotDS == nTotDS2)

cat(paste0(length(unique(all_matching_dt$ref_dataset))), "\n")

 # when there was no match -> NA was set; set 0
all_matching_dt$overlapBp_ratio[is.na(all_matching_dt$overlapBp_ratio)] <- 0
all_matching_dt$nOverlapGenes_ratio[is.na(all_matching_dt$nOverlapGenes_ratio)] <- 0

stopifnot(all_matching_dt$overlapBp_ratio >= 0 & all_matching_dt$overlapBp_ratio <= 1)
stopifnot(all_matching_dt$nOverlapGenes_ratio >= 0 & all_matching_dt$nOverlapGenes_ratio <= 1)

cat(paste0(length(unique(all_matching_dt$ref_dataset))), "\n")
# 
outFile <- file.path(outFolder, paste0("ratio_matching_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(list(
  overlapBp_ratio = all_matching_dt$overlapBp_ratio,
  nOverlapGenes_ratio = all_matching_dt$nOverlapGenes_ratio
  ),
  plotTit = paste0("ratio overlap across datasets")
)
mtext(text = paste0("nDS = ", nTotDS), font=3, side=3)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

overlapType <- "overlapBp"
overlapType <- "nOverlapGenes"
all_overlapTypes <- c("overlapBp", "nOverlapGenes")


all_levels <- sort(as.character(unique(all_matching_dt$ref_dataset)))

# all_matching_dt <- all_matching_dt[all_matching_dt$ref_dataset == all_matching_dt$ref_dataset[1],]

myxlab <- "ref_dataset_label"
myylab <- "matching_dataset"


for(overlapType in all_overlapTypes){
  
  # select the best match for each TAD
  cat("... select best for each TAD\n")
  
  cat(paste0(Sys.time(), " - "))
  
  
  all_matching_dt_nbr <- all_matching_dt[,c("ref_dataset", "refID")]
  all_matching_dt_nbr <- unique(all_matching_dt_nbr)
  nds_dt <- aggregate(refID ~ ref_dataset, data = all_matching_dt_nbr, FUN=length)
  colnames(nds_dt)[colnames(nds_dt) == "refID"] <- "nTADs"
  
  
  all_matching_dt_pairAgg <- aggregate(as.formula(paste0(overlapType, "_ratio ~ ref_dataset + matching_dataset")),
                                       FUN = mean, data = all_matching_dt)
  
  newcolname <- paste0(overlapType, "_ratio_avg")
  colnames(all_matching_dt_pairAgg)[colnames(all_matching_dt_pairAgg) == paste0(overlapType, "_ratio")] <- newcolname
  
  plot_dt <- merge(all_matching_dt_pairAgg, nds_dt, all.x=TRUE, by ="ref_dataset")
  stopifnot(nrow(plot_dt) == nrow(all_matching_dt_pairAgg))
  
  plot_dt$ref_dataset <- factor(plot_dt$ref_dataset, levels = all_levels)
  plot_dt$matching_dataset <- factor(plot_dt$matching_dataset, levels = all_levels)
  
  plot_dt <- plot_dt[order(as.numeric(plot_dt$ref_dataset)),]
  
  plot_dt$ref_dataset_label <- paste0(as.character(plot_dt$ref_dataset), "\n", plot_dt$nTADs)
  plot_dt$ref_dataset_label <- factor(plot_dt$ref_dataset_label, levels=unique(as.character(plot_dt$ref_dataset_label)))
  
  plot_var <- newcolname
  
  plot_dt[,paste0(plot_var)] <- round(plot_dt[,paste0(plot_var)], 2)
  
  agg_matching_plot <- ggplot(plot_dt, 
                             aes_string(x="ref_dataset_label", y="matching_dataset", fill = paste0(plot_var)))+
    geom_tile(color = "white")+
    geom_text(aes_string(label = paste0(plot_var)), color = "black", size = 3, fontface="bold")+
    scale_fill_gradientn(colours = hm.palette(100))+
    guides(fill = FALSE)+ 
    labs(title = paste0(overlapType, " - ", plot_var),
         subtitle = paste0("nDS = ", nTotDS),
         x = paste0(myxlab),
         y = paste0(myylab)
    )+
    theme(
      axis.text.x = element_text(size=12, angle=90, hjust = 1, vjust=0.5),
      axis.text.y = element_text(size=12),
      axis.title = element_text(size=14, face="bold"),
      plot.title = element_text(size=18, face="bold", hjust=0.5),
      plot.subtitle = element_text(size=16, face="italic", hjust=0.5),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      # panel.background = element_blank()
      panel.background=element_rect(fill="grey50", colour="grey50")
    )
  # coord_equal() # make the x and y at the same scale => same as coord_fixed(ratio = 1)
  outFile <- file.path(outFolder, paste0( "allMatchingMean", "_", overlapType, "_ratio_heatmap.", plotType))
  ggsave(agg_matching_plot, filename = outFile, height = myHeightGG, width = myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))  

  


  
  
  # too slow... 
  # try with data.table; see other version in _v2 and _v3 version
  # the v3 version will yield similar results, but will discard the NA rows directly
  
  all_matching_dt <- all_matching_dt[order(all_matching_dt[,paste0(overlapType)], decreasing=TRUE),]
  
  cat(paste0(length(unique(all_matching_dt$ref_dataset))), "\n")
  
  all_matching_dt_DT <- data.table(all_matching_dt)

  
  cat(paste0(Sys.time(), " - "))
  

  all_matching_onlyBest_dt <- all_matching_dt_DT[, lapply(.SD, head, n = 1L), by = list(ref_dataset, refID)] # NB: much faster than aggregate!
  
  cat(paste0(length(unique(all_matching_onlyBest_dt$ref_dataset))), "\n")
  
  
  # as they are sorted -> take the 1st
  # NB: ties are not handled
  all_matching_onlyBest_dt <- data.frame(all_matching_onlyBest_dt)
  
  cat(paste0(Sys.time(), "\n"))
  cat(paste0(length(unique(all_matching_onlyBest_dt$ref_dataset))), "\n")
  
  rownames(all_matching_onlyBest_dt) <- NULL
  all_matching_onlyBest_dt <- all_matching_onlyBest_dt[order(all_matching_onlyBest_dt$ref_dataset, all_matching_onlyBest_dt$refID),]
  all_matching_onlyBest_dt <- data.frame(all_matching_onlyBest_dt)
  cat(paste0(Sys.time(), "\n"))
  
  cat(paste0(length(unique(all_matching_onlyBest_dt$ref_dataset))), "\n")
  
  rownames(all_matching_onlyBest_dt) <- NULL
  all_matching_onlyBest_dt <- all_matching_onlyBest_dt[order(all_matching_onlyBest_dt$ref_dataset, all_matching_onlyBest_dt$refID),]
  
  cat(paste0(length(unique(all_matching_onlyBest_dt$ref_dataset))), "\n")
  
  
  outFile <- file.path(outFolder, paste0(overlapType, "_all_matching_onlyBest_dt.Rdata"))
  save(all_matching_onlyBest_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  

  stopifnot(!(duplicated(paste0(all_matching_onlyBest_dt[, c("ref_dataset")], all_matching_onlyBest_dt[, c("refID")] ))))
  
  all_matching_onlyBest_dt_onlyRef <- all_matching_onlyBest_dt[,c("ref_dataset", "refID")]
  all_matching_onlyBest_dt_onlyRef <- unique(all_matching_onlyBest_dt_onlyRef)
  nds_dt <- aggregate(refID ~ ref_dataset, data = all_matching_onlyBest_dt_onlyRef, FUN=length)
  colnames(nds_dt)[colnames(nds_dt) == "refID"] <- "nTADs"
  
  cat(paste0(length(unique(all_matching_onlyBest_dt_onlyRef$ref_dataset))), "\n")
  
  cat("... aggregation nbr\n")
    
  # agg_matching_nbr_dt <- aggregate( as.formula(paste0(overlapType, "_ratio ~ ref_dataset + matching_dataset")) , 
  #                                  data=all_matching_onlyBest_dt, FUN=length)
  # do not take the row where this is 0 -> because the 1st line was taken
  # if a TAD had no match at all -> will penalize the first dataset that pop up
  
  cat(paste0(length(unique(all_matching_onlyBest_dt$ref_dataset))), "\n")
  
  
  noMatch_dt <- all_matching_onlyBest_dt[all_matching_onlyBest_dt[paste0(overlapType, "_ratio")] == 0,]
  outFile <- file.path(outFolder, "noMatch_dt.txt")
  write.table(noMatch_dt, file = outFile, col.names=T, row.names=F, quote=F, append=F, sep="\t")
  cat(paste0("... written: ", outFile, "\n"))
  
  cat(paste0("... TADs with no match at all: ", sum(all_matching_onlyBest_dt[paste0(overlapType, "_ratio")] == 0), "\n"))
  
  agg_matching_nbr_dt <- aggregate( as.formula(paste0(overlapType, "_ratio ~ ref_dataset + matching_dataset")) , 
                                    data=all_matching_onlyBest_dt[all_matching_onlyBest_dt[paste0(overlapType, "_ratio")] > 0,], 
                                    FUN=length)
  
  colnames(agg_matching_nbr_dt)[colnames(agg_matching_nbr_dt) == paste0(overlapType, "_ratio")] <- "nBestMatchs"
  
  cat(paste0(length(unique(agg_matching_nbr_dt$ref_dataset))), "\n")
  
  
  cat("... aggregation ratio\n")
  # agg_matching_avgRatio_dt <- aggregate(as.formula(paste0(overlapType, "_ratio ~ ref_dataset + matching_dataset")) , 
  #                                       data=all_matching_onlyBest_dt, 
  #                                       FUN=mean)
  agg_matching_avgRatio_dt <- aggregate(as.formula(paste0(overlapType, "_ratio ~ ref_dataset + matching_dataset")) , 
                                        data=all_matching_onlyBest_dt[all_matching_onlyBest_dt[paste0(overlapType, "_ratio")] > 0,], 
                                        FUN=mean)
  colnames(agg_matching_avgRatio_dt)[colnames(agg_matching_avgRatio_dt) == paste0(overlapType, "_ratio")] <- "meanBestRatio"
  
  
  cat(paste0(length(unique(agg_matching_avgRatio_dt$ref_dataset))), "\n")
  
  agg_matching_nbr_ratio_dt <- merge(agg_matching_nbr_dt, nds_dt, by="ref_dataset", all =TRUE)
  stopifnot(nrow(agg_matching_nbr_ratio_dt) == nrow(agg_matching_nbr_dt))
  agg_matching_nbr_ratio_dt$nBestMatchs_ratio <- agg_matching_nbr_ratio_dt$nBestMatchs/agg_matching_nbr_ratio_dt$nTADs
  
  stopifnot( agg_matching_nbr_ratio_dt$nBestMatchs_ratio >= 0)
  stopifnot( agg_matching_nbr_ratio_dt$nBestMatchs_ratio <= 1)
  
  cat(paste0(length(unique(agg_matching_nbr_ratio_dt$ref_dataset))), "\n")
  
  # agg_matching_dt <- aggregate(matchingID_maxOverlapBp ~ ref_dataset + matching_dataset, data=all_matching_onlyBest_dt, FUN=length)
  # 
  # tot_matching_dt <- data.frame(ref_dataset = as.character(names(table(all_matching_onlyBest_dt$ref_dataset))), 
  #                               totTADs = as.numeric(table(all_matching_onlyBest_dt$ref_dataset)), stringsAsFactors = FALSE)
  # 
  # agg_dt <- merge(agg_matching_dt, tot_matching_dt, by="ref_dataset")
  # 
  # agg_dt <- agg_dt[order(agg_dt$ref_dataset),]
  
  agg_matching_nbr_ratio_dt$ref_dataset <- factor(agg_matching_nbr_ratio_dt$ref_dataset, levels = all_levels)
  agg_matching_nbr_ratio_dt$matching_dataset <- factor(agg_matching_nbr_ratio_dt$matching_dataset, levels = all_levels)
  
  stopifnot(!is.na(agg_matching_nbr_ratio_dt$ref_dataset))
  stopifnot(!is.na(agg_matching_nbr_ratio_dt$matching_dataset))
  
  
  agg_matching_nbr_ratio_dt <- agg_matching_nbr_ratio_dt[order(as.numeric(agg_matching_nbr_ratio_dt$ref_dataset)),]
  
  agg_matching_nbr_ratio_dt$ref_dataset_label <- paste0(as.character(agg_matching_nbr_ratio_dt$ref_dataset), "\n", 
                                                        agg_matching_nbr_ratio_dt$nTADs)
  
  agg_matching_nbr_ratio_dt$ref_dataset_label <- factor(agg_matching_nbr_ratio_dt$ref_dataset_label, 
                                                  levels = unique(as.character(agg_matching_nbr_ratio_dt$ref_dataset_label)))
  
  plot_var <- "nBestMatchs"
  plot_vars <- c("nBestMatchs", "nBestMatchs_ratio")
  for(plot_var in plot_vars) {
    
    agg_matching_nbr_ratio_dt[,paste0(plot_var)] <- round(agg_matching_nbr_ratio_dt[,paste0(plot_var)], 2)
    
    
    ds_matching_plot <- ggplot(agg_matching_nbr_ratio_dt, 
                      aes_string(x="ref_dataset_label", y="matching_dataset", fill = paste0(plot_var)))+
      geom_tile(color = "white")+
      geom_text(aes_string(label = paste0(plot_var)), color = "black", size = 3, fontface="bold")+
      scale_fill_gradientn(colours = hm.palette(100))+
      guides(fill = FALSE)+ 
      labs(title = paste0(overlapType, " - ", plot_var),
           subtitle = paste0("nDS = ", nTotDS),
           x = paste0(myxlab),
           y = paste0(myylab)
      )+
      theme(
        axis.text.x = element_text(size=12, angle=90, hjust = 1, vjust=0.5),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(size=18, face="bold", hjust=0.5),
        plot.subtitle = element_text(size=16, face="italic", hjust=0.5),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        # panel.background = element_blank()
        panel.background=element_rect(fill="grey50", colour="grey50")
      )
    # coord_equal() # make the x and y at the same scale => same as coord_fixed(ratio = 1)
    outFile <- file.path(outFolder, paste0( "bestMatching", "_", plot_var, "_", overlapType, "_heatmap.", plotType))
    ggsave(ds_matching_plot, filename = outFile, height = myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))  
  } # end-for iterating over plot_vars
  
  
  
  

} # end-for iterating over overlap types



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
