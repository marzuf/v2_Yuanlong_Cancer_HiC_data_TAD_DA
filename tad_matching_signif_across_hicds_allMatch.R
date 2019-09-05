options(scipen=100)

# Rscript tad_matching_across_hicds_allMatch.R

script_name <- "tad_matching_signif_across_hicds_allMatch.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

buildTable <- TRUE

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 12


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
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH")
dir.create(outFolder, recursive = TRUE)

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

# => best TAD matching
# in # of genes
# in bp

final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))

signif_column <- "adjPvalComb"

signifThresh <- 0.01

signifcol <- paste0(signif_column, "_", signifThresh)

final_dt[, paste0(signifcol)] <- final_dt[, paste0(signif_column)] <= signifThresh


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
      
      
      signif_tads <- final_dt$region[final_dt$hicds == hicds & final_dt$exprds == exprds & final_dt[, paste0(signifcol)] ]
      
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
      
      # take only the data from signif TADs
      stopifnot(setequal(ref_g2t_dt$entrezID, geneList))
      signif_ref_g2t_dt <- ref_g2t_dt[ref_g2t_dt$region %in% signif_tads,]
      signif_geneList <- geneList[geneList %in% as.character(signif_ref_g2t_dt$entrezID)]
      
      signif_ref_tadpos_dt <- ref_tadpos_dt[ref_tadpos_dt$region %in% signif_tads,]
      
      stopifnot(setequal(signif_ref_g2t_dt$entrezID, signif_geneList))
      stopifnot(setequal(signif_ref_g2t_dt$region, signif_ref_tadpos_dt$region))
      
      list(
        dataset_geneList = signif_geneList,
        dataset_g2t_dt = signif_ref_g2t_dt ,
        dataset_tadpos_dt = signif_ref_tadpos_dt
      )
    }
    names(exprds_list) <- file.path(hicds, all_exprds[[paste0(hicds)]])
    exprds_list
  }
  all_data_list <- unlist(all_data_list, recursive=FALSE)
  stopifnot(length(all_data_list) == length(all_datasets))
  save(all_data_list, file = file.path(outFolder, "all_data_list.Rdata"), version=2)
  ####################################################################################################################################### >>> collect # of TADs and TAD size
  
  
  cat("... start matching\n")
  
  ref_dataset = all_datasets[1]
  # all_datasets=all_datasets[1:2]
  all_signif_matching_dt <- foreach(ref_dataset = all_datasets, .combine='rbind') %dopar% {
    
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
    query_dataset=query_datasets[15]
    
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
      
      if(length(overlap_GR) == 0) return(NULL)
      
      IDoverlap_hits_all <- findOverlaps(query=query_GR,
                                         subject=query_GR)
      IDoverlap_hits <- IDoverlap_hits_all[queryHits(IDoverlap_hits_all) != subjectHits(IDoverlap_hits_all)] 
      
      refID <- names(ref_GR[subjectHits(overlap_GR)])
      queryID <- names(query_GR[queryHits(overlap_GR)])
      
      overlapDT_bp <- data.frame(
        refID = refID,
        queryID = queryID,
        overlapBp = width(pintersect(ref_GR[refID], query_GR[queryID])),
        stringsAsFactors = FALSE)
      stopifnot(overlapDT_bp$refID %in% ref_tadpos_dt$region)
      stopifnot(overlapDT_bp$queryID %in% query_tadpos_dt$region)
      
      nmat <- length(refID)
      stopifnot(nmat == length(queryID))
      
      i=1
      i=2
      overlapDT_genes <- foreach(i=1:nmat, .combine='rbind') %do% {
        ref_genes <- final_dt$region_genes[final_dt$hicds == dirname(ref_dataset) & 
                                             final_dt$exprds == basename(ref_dataset) & 
                                             final_dt$region == refID[i]]
        ref_genes_ul <- unlist(strsplit(ref_genes, split=","))
        query_genes <- final_dt$region_genes[final_dt$hicds == dirname(query_dataset) & 
                                               final_dt$exprds == basename(query_dataset) & 
                                               final_dt$region == queryID[i]]
        query_genes_ul <- unlist(strsplit(query_genes, split=","))
        common_genes <- intersect(ref_genes_ul, query_genes_ul)
        data.frame(
          refID = refID[i],
          queryID = queryID[i],
          overlapGenes = paste0(common_genes, collapse=","),
          nOverlapGenes = length(common_genes),
        stringsAsFactors=FALSE
        )
      }
      
      
      overlapDT <- merge(overlapDT_bp, overlapDT_genes, by = c("refID", "queryID"), all=TRUE)
      
      overlapDT <- overlapDT[order(overlapDT$refID, overlapDT$queryID, overlapDT$overlapBp, overlapDT$nOverlapGenes),]
      
      all_cols <- colnames(overlapDT)
      
      overlapDT$ref_dataset <- ref_dataset
      overlapDT$query_dataset <- query_dataset
      
      overlapDT[, c("ref_dataset", "query_dataset", all_cols)]
      
      # all_match_dt$ref_dataset <- ref_dataset
      # all_match_dt$matching_dataset <- query_dataset
      
      
      # DO NOT TAKE THE MAX OVERLAP, RETAIN ALL
      # take the max overlap
      # maxBp_overlapDT <- do.call(rbind, by(overlapDT, overlapDT$refID, function(subDT){  # WARNING: if exactly same bp overlap -> will take the 1st one...
      #   subDT[which.max(subDT$overlapBp),]
      # }))
      # stopifnot(rownames(maxBp_overlapDT) == maxBp_overlapDT$refID)
      # rownames(maxBp_overlapDT) <- NULL
      # stopifnot(setequal(overlapDT$refID, maxBp_overlapDT$refID))
      
      
      ### PREPARE THE GENE OVERLAP
      # cat("... preparing gene matching overlap\n")
      # ref_tad_genes <- split(ref_g2t_dt$entrezID, ref_g2t_dt$region)
      # query_tad_genes <- split(query_g2t_dt$entrezID, query_g2t_dt$region)
      
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
      # tad_genes=ref_tad_genes[1]
      # all_gene_overlap_list <- foreach(tad_genes=ref_tad_genes) %dopar% {
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
      # }
      # names(all_gene_overlap_list) <- names(ref_tad_genes)
      # 
      # gene_overlap_list <- Filter(function(x) x[1] != "-Inf", all_gene_overlap_list)
      # stopifnot(lengths(gene_overlap_list) == 3)
      # 
      # maxGenes_overlapDT <- data.frame(do.call(rbind, gene_overlap_list))
      # maxGenes_overlapDT$refID <- rownames(maxGenes_overlapDT)
      # colnames(maxGenes_overlapDT)[1] <- "queryID"
      # colnames(maxGenes_overlapDT)[2] <- "nOverlapGenes"
      # colnames(maxGenes_overlapDT)[3] <- "overlapGenes"
      # maxGenes_overlapDT$nOverlapGenes <- as.numeric(as.character(maxGenes_overlapDT$nOverlapGenes))
      # stopifnot(maxGenes_overlapDT$refID %in% ref_g2t_dt$region)
      # stopifnot(maxGenes_overlapDT$queryID %in% query_g2t_dt$region)
      # stopifnot(!is.na(maxGenes_overlapDT$nOverlapGenes))  
      # 
      # colnames(maxBp_overlapDT)[colnames(maxBp_overlapDT) == "queryID"] <- "matchingID_maxOverlapBp"
      # colnames(maxGenes_overlapDT)[colnames(maxGenes_overlapDT) == "queryID"] <-  "matchingID_maxOverlapGenes"
      # 
      # 
      # stopifnot(maxBp_overlapDT$refID %in% ref_tadpos_dt_info$refID)
      # stopifnot(maxGenes_overlapDT$refID %in% ref_tadpos_dt_info$refID)
      # 
      # 
      # tmp_dt <- merge(maxBp_overlapDT, maxGenes_overlapDT, by ="refID", all = TRUE)
      # 
      # all_match_dt <- merge(tmp_dt, ref_tadpos_dt_info[,c("refID", "ref_totBp", "ref_nGenes"), drop=FALSE], by="refID", all=TRUE)
      # stopifnot(nrow(all_match_dt) == nrow(ref_tadpos_dt_info))
      # 
      # first_cols <- colnames(all_match_dt)[grepl("^ref",colnames(all_match_dt))]
      # last_cols <- colnames(all_match_dt)[!grepl("^ref",colnames(all_match_dt))]
      # all_match_dt <- all_match_dt[,c(first_cols, last_cols)]
      # all_cols <- colnames(all_match_dt)
      
      # all_match_dt$ref_dataset <- ref_dataset
      # all_match_dt$matching_dataset <- query_dataset
      # rownames(all_match_dt) <- NULL
      # all_match_dt[,c("ref_dataset", "matching_dataset", all_cols)]
      
    } # end iterating over query_datasets
    ref_dt
  } # end iterating over ref_datasets
  
  outFile <- file.path(outFolder, paste0("all_signif_matching_dt_", signifcol, ".Rdata"))
  save(all_signif_matching_dt, file=outFile, version=2)  
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, paste0("all_signif_matching_dt_", signifcol, ".Rdata"))
  cat("... load data\n")
  all_signif_matching_dt <- get(load(outFile))
}

# nrow(all_signif_matching_dt) # 12132
# length(unique(paste0(all_signif_matching_dt$ref_dataset, all_signif_matching_dt$refID))) # 1303
# length(unique(paste0(all_signif_matching_dt$query_dataset, all_signif_matching_dt$queryID))) # 1303


all_signif_matching_dt$refID_full <- paste0(all_signif_matching_dt$ref_dataset, "-", all_signif_matching_dt$refID)
all_signif_matching_dt$queryID_full <- paste0(all_signif_matching_dt$query_dataset, "-", all_signif_matching_dt$queryID)


# get the number of matches, to have an idea of multiple matching
nMatchByTad <- aggregate(overlapBp ~ refID + ref_dataset, FUN=length, data = all_signif_matching_dt)
range(nMatchByTad$overlapBp) # 1 - 55

all_signif_tads_with_match <- unique(all_signif_matching_dt$refID_full)
length(all_signif_tads_with_match) # 1303 # cmp with 
sum(final_dt[,paste0(signifcol)]) # 1453 3 -> 150 match with no other

tad=all_signif_tads_with_match[1]

# set_of_tads <- by(all_signif_matching_dt, all_signif_matching_dt$refID_full, function(x) {
#   myref <- unique(x$refID_full)
#   stopifnot(length(myref) == 1)
#   c(myref, as.character(x$queryID_full))
# })


set_of_tads2 <- foreach(tad = all_signif_tads_with_match) %dopar% {
  matching_tads <- all_signif_matching_dt$queryID_full[all_signif_matching_dt$refID_full == tad]
  c(tad, matching_tads)
}
names(set_of_tads2) <- all_signif_tads_with_match

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
