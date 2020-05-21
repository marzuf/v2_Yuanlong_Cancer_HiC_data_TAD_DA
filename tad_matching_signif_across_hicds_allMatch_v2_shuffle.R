options(scipen=100)

set.seed(20052020)

# v2 version
# -> one single script for all and for cmpType
# -> if share more than 80% of the genes -> merge conserved regions
# -> count as conserved in one dataset at the exprds level
# -> min conserved region

SSHFS=F

# Rscript tad_matching_signif_across_hicds_allMatch_v2_shuffle.R
# Rscript tad_matching_signif_across_hicds_allMatch_v2_shuffle.R norm_vs_tumor
# Rscript tad_matching_signif_across_hicds_allMatch_v2_shuffle.R subtypes
# Rscript tad_matching_signif_across_hicds_allMatch_v2_shuffle.R wt_vs_mut

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  data_cmpType <- ""
  file_prefix <- ""
} else if(length(args) == 1) {
  data_cmpType <- args[1]  
  stopifnot(data_cmpType %in% c("norm_vs_tumor", "subtypes", "wt_vs_mut"))
  file_prefix <- paste0(data_cmpType, "_")
} else {
 stop("error\n") 
}

nPermut <- 1000

script_name <- "tad_matching_signif_across_hicds_allMatch_v2_shuffle.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("plot_lolliTAD_funct.R")
source("my_heatmap.2.R")

buildTable <- TRUE

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4


source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

script0_name <- "0_prepGeneData"

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)

all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]

stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_SHUFFLE", data_cmpType)
dir.create(outFolder, recursive = TRUE)

logFile <- file.path(outFolder, "tad_matching_signif_across_hicds_logFile.txt")
file.remove(logFile)

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

# to retrieve the colors of cmpType for plotting heatmap
dataset_dt <- data.frame(dataset = all_datasets, hicds = dirname(all_datasets),exprds = basename(all_datasets),stringsAsFactors = FALSE)
dataset_dt$cmpType <- all_cmps[paste0(dataset_dt$exprds)]
dataset_dt$subtype_col <- all_cols[paste0(dataset_dt$cmpType)]
stopifnot(!is.na(dataset_dt$subtype_col))
dataset_dt <- dataset_dt[order(dataset_dt$cmpType, dataset_dt$hicds, dataset_dt$exprds),]

if(length(args)>0) {
  
  stopifnot(data_cmpType %in% dataset_dt$cmpType)
  dataset_dt <- dataset_dt[dataset_dt$cmpType == data_cmpType,]
  
  all_datasets <- all_datasets[basename(all_datasets) %in% dataset_dt$exprds & dirname(all_datasets) %in% dataset_dt$hicds]
  stopifnot(length(all_datasets) > 0)
  
  all_hicds <- all_hicds[all_hicds %in% dirname(all_datasets)]
  all_exprds <- lapply(all_exprds, function(x) x[x %in% basename(all_datasets)])
  
  stopifnot(length(unlist(all_exprds)) == length(all_datasets))
}


dataset_dt$dataset_label <- gsub("/", "\n", dataset_dt$dataset)
ds_label_levels <- dataset_dt$dataset_label
ds_levels <- dataset_dt$dataset

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

final_dt <- final_dt[final_dt$hicds %in% dirname(all_datasets) & final_dt$exprds %in% basename(all_datasets),]
stopifnot(nrow(final_dt) > 0)

stopifnot(length(unique(paste0(final_dt$hicds, final_dt$exprds))) == length(all_datasets))

minOverlapBpRatio <- 0.8
minIntersectGenes <- 3
gene_matching_fuse_threshold <- 0.8

txt <- paste0("> signif_column\t=\t", signif_column, "\n")
cat(txt)
cat(txt, file=logFile, append=T)
txt <- paste0("> signifThresh\t=\t", signifThresh, "\n")
cat(txt)
cat(txt, file=logFile, append=T)
txt <- paste0("> minOverlapBpRatio\t=\t", minOverlapBpRatio, "\n")
cat(txt)
cat(txt, file=logFile, append=T)
txt <- paste0("> minIntersectGenes\t=\t", minIntersectGenes, "\n")
cat(txt)
cat(txt, file=logFile, append=T)

true_signif_tads <- file.path(final_dt$hicds[final_dt[, paste0(signifcol)] ],
                         final_dt$exprds[final_dt[, paste0(signifcol)] ],
                         final_dt$region[final_dt[, paste0(signifcol)] ])
# 
# outFile <- file.path(outFolder, paste0(file_prefix, "signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
# save(signif_tads, file=outFile, version=2)
# cat(paste0("... written: ", outFile, "\n"))

####################################################################################################################################### >>> prepare the data
  
  foo_permut <- foreach(i_perm=1:nPermut) %dopar% {
    
    dir.create(file.path(outFolder,i_perm), recursive=TRUE)
    
    if(buildTable) {
      cat("... start preparing data before matching \n")
      
    
    hicds = all_hicds[1]
    all_data_list <- foreach(hicds = all_hicds) %do% {
      
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
        
        signif_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds & final_dt[, paste0(signifcol)], ]
        not_signif_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds & !final_dt[, paste0(signifcol)], ]
        
        # save(signif_dt, file="signif_dt.Rdata",version=2)
        # save(not_signif_dt, file="not_signif_dt.Rdata",version=2)
        # save(not_signif_dt, file="final_dt.Rdata",version=2)
        
        stopifnot(nrow(signif_dt) + nrow(not_signif_dt) == nrow(final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds,]))
        stopifnot(length(intersect(signif_dt$region, not_signif_dt$region)) == 0)
        
        # sample in the not signif TADs as many as in the signif TADs
        not_signif_tads <- not_signif_dt$region
        signif_tads <- signif_dt$region
        
        pseudosignif_tads <- sample(not_signif_tads, size=length(signif_tads))
        stopifnot(length(pseudosignif_tads) == length(signif_tads))
        stopifnot(!pseudosignif_tads %in% signif_tads)
        
        # save(pseudosignif_tads, file="pseudosignif_tads.Rdata",version=2)
        
        # x=paste0(outFolder,
        #           i_perm,
        #           paste0(file_prefix, hicds, "_", exprds, "_pseudosignif_tads", 
        #                  signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
        # save(x, file="x.Rdata", version=2)

                
        # outFile <- file.path(outFolder,
        #                      i_perm,
        #                      paste0(file_prefix, hicds, "_", exprds, "_pseudosignif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
        # save(pseudosignif_tads, file=outFile, version=2)
        # cat(paste0("... written: ", outFile, "\n"))
        signif_tads <- pseudosignif_tads
        
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
    save(all_data_list, file = file.path(outFolder, i_perm, "all_data_list.Rdata"), version=2)
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
      
      ref_tad_size <- setNames(c(ref_tadpos_dt$end-ref_tadpos_dt$start+1), ref_tadpos_dt$region)
      
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
      query_dataset=query_datasets[1]
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
        
        stopifnot(refID %in% names(ref_tad_size))
        
        overlapDT_bp <- data.frame(
          refID = refID,
          queryID = queryID,
          overlapBp = width(pintersect(ref_GR[refID], query_GR[queryID])),
          overlapBpRatio = width(pintersect(ref_GR[refID], query_GR[queryID]))/ref_tad_size[refID],
          stringsAsFactors = FALSE)
        stopifnot(overlapDT_bp$refID %in% ref_tadpos_dt$region)
        stopifnot(overlapDT_bp$queryID %in% query_tadpos_dt$region)
        stopifnot(!is.na(overlapDT_bp))
        stopifnot(overlapDT_bp$overlapBpRatio >= 0)
        stopifnot(overlapDT_bp$overlapBpRatio <= 1)
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
      } # end iterating over query_datasets
      ref_dt
    } # end iterating over ref_datasets
    
    outFile <- file.path(outFolder, 
                         i_perm,
                         paste0(file_prefix, "all_signif_matching_dt_", signifcol, ".Rdata"))
    save(all_signif_matching_dt, file=outFile, version=2)  
    cat(paste0("... written: ", outFile, "\n"))
    
  } else {
    outFile <- file.path(outFolder, 
                         i_perm,
                         paste0(file_prefix, "all_signif_matching_dt_", signifcol, ".Rdata"))
    cat("... load data\n")
    all_signif_matching_dt <- get(load(outFile))
    outFile <- file.path(outFolder, 
                         i_perm,
                         "all_data_list.Rdata")
    all_data_list <- get(load(outFile))
  }
  all_signif_matching_dt$refID_full <- file.path(all_signif_matching_dt$ref_dataset, all_signif_matching_dt$refID)
  all_signif_matching_dt$queryID_full <- file.path(all_signif_matching_dt$query_dataset, all_signif_matching_dt$queryID)
  
  
  #######################################################################################################
  ####################################################################################################### CONSERVED REGIONS REFINEMENT
  #######################################################################################################
  
  ### > Filter1: retain only matches with >= *minOverlapBpRatio* (80%) bp overlap
  
  txt <- paste0("> FILTER 1 - # of match\t=\t", nrow(all_signif_matching_dt), "\n")
  cat(txt)
  cat(txt, file=logFile, append=T)
  all_signif_matching_dt  <- all_signif_matching_dt[all_signif_matching_dt$overlapBpRatio >= minOverlapBpRatio,]
  txt <- paste0("> FILTER 1 - # of match with overlap bp >= ", minOverlapBpRatio, "\t=\t", nrow(all_signif_matching_dt), "\n")
  cat(txt)
  cat(txt, file=logFile, append=T)
  
  all_signif_tads_with_match <- unique(all_signif_matching_dt$refID_full)
  length(all_signif_tads_with_match) # 
  sum(final_dt[,paste0(signifcol)]) # 
  
  tad=all_signif_tads_with_match[1]
  set_of_tads <- foreach(tad = all_signif_tads_with_match) %dopar% {
    # all_signif_matching_dt[all_signif_matching_dt$refID_full == tad,]
    matching_tads <- all_signif_matching_dt$queryID_full[all_signif_matching_dt$refID_full == tad]
    sort(c(tad, matching_tads))
  }
  names(set_of_tads) <- all_signif_tads_with_match
  
  ### Remove duplicated sets
  
  txt <- paste0("... remove duplicated sets:\t", length(set_of_tads) , " -> ")
  cat(txt)
  cat(txt, file=logFile, append=T)
  set_of_tads <- unique(set_of_tads) # loose the name !
  length(set_of_tads) # 368 ! # 514
  txt <- paste0(length(set_of_tads) , "\n")
  cat(txt)
  cat(txt, file=logFile, append=T)
  
  ### Remove nested sets
  not_nested_sets <- unlist(lapply(1:length(set_of_tads), function(i_set_tads) {
    set_tads <- set_of_tads[[i_set_tads]]
    !any(sapply(set_of_tads[-i_set_tads], function(x) all(set_tads %in% x) ))
  }))
  length(set_of_tads)
  not_nested_set_of_tads <- set_of_tads[not_nested_sets]
  length(not_nested_set_of_tads)
  stopifnot(!is.na(not_nested_set_of_tads))
  stopifnot(setequal(unlist(set_of_tads), unlist(not_nested_set_of_tads)))
  stopifnot(length(set_of_tads) >= length(not_nested_set_of_tads))
  stopifnot(names(not_nested_set_of_tads) %in% names(set_of_tads))
  
  txt <- paste0("... remove nested sets:\t", length(set_of_tads) , " -> ")
  cat(txt)
  cat(txt, file=logFile, append=T)
  
  set_of_tads <- not_nested_set_of_tads
  txt <- paste0(length(set_of_tads) , "\n")
  cat(txt)
  cat(txt, file=logFile, append=T)
  
  
  ### > Filter2: retain only "conserved regions" with >= *minIntersectGenes*(3) genes at the intersect
  txt <- paste0("> FILTER 2 - # of conserved regions\t=\t", length(set_of_tads), "\n")
  cat(txt)
  cat(txt, file=logFile, append=T)
  
  names(set_of_tads) <- paste0("conserved_region_", seq_along(set_of_tads))
  cr=names(set_of_tads)[1]
  set_of_tads_intersect_genes <- foreach(cr = names(set_of_tads)) %dopar% {
    all_regions <- set_of_tads[[paste0(cr)]]
    ds_genes <- sapply(all_regions, function(tad) {
      dt <- all_data_list[[paste0(dirname(tad))]][["dataset_g2t_dt"]]
      stopifnot(basename(tad) %in% dt$region)
      as.character(dt$entrezID[dt$region == basename(tad)])
    })
    is_genes <- Reduce(intersect, ds_genes)
    stopifnot(is_genes %in% gff_dt$entrezID)
    stopifnot(is_genes %in% names(entrez2symb))
    setNames(entrez2symb[paste0(is_genes)], is_genes)
  }
  names(set_of_tads_intersect_genes) <- names(set_of_tads)
  nIntersectGenesByRegions <- lengths(set_of_tads_intersect_genes)
  regionsWithMinGenes <- names(nIntersectGenesByRegions) [nIntersectGenesByRegions >= minIntersectGenes]
  stopifnot(regionsWithMinGenes %in% names(set_of_tads))
  all_set_of_tads  <- set_of_tads
  set_of_tads <- set_of_tads[paste0(regionsWithMinGenes)]
  stopifnot(setequal(regionsWithMinGenes, names(set_of_tads)))
  
  ### !!! >>> v2 update here 
  ### Merge conserved regions that have >= *gene_matching_fuse_threshold* % (80%) gene overlap
  set_of_tads_intersect_genes <- foreach(cr = names(set_of_tads)) %dopar% {
    all_regions <- set_of_tads[[paste0(cr)]]
    ds_genes <- sapply(all_regions, function(tad) {
      dt <- all_data_list[[paste0(dirname(tad))]][["dataset_g2t_dt"]]
      stopifnot(basename(tad) %in% dt$region)
      as.character(dt$entrezID[dt$region == basename(tad)])
    })
    is_genes <- Reduce(intersect, ds_genes)
    stopifnot(is_genes %in% gff_dt$entrezID)
    stopifnot(is_genes %in% names(entrez2symb))
    setNames(entrez2symb[paste0(is_genes)], is_genes)
  }
  names(set_of_tads_intersect_genes) <- names(set_of_tads)
  
  cr = names(set_of_tads_intersect_genes)[1]
  set_of_tads_matching_ratio <- foreach(cr = names(set_of_tads_intersect_genes)) %dopar% {
    region_genes <- set_of_tads_intersect_genes[[paste0(cr)]]
    matching_ratio_other_regions <- unlist(lapply(set_of_tads_intersect_genes, function(x) sum(region_genes %in% x)/length(region_genes)))
    tofuse <- names(matching_ratio_other_regions)[matching_ratio_other_regions >= gene_matching_fuse_threshold]
    tofuse
  }
  names(set_of_tads_matching_ratio) <- names(set_of_tads_intersect_genes)
  list_to_fuse <- Filter(function(x) length(x) > 1, set_of_tads_matching_ratio)  
  names(list_to_fuse) <- NULL
  list_to_fuse <- unique(list_to_fuse)
  
  # list_to_fuse <- list(
  #   cr1 = c("cr1", "cr2", "cr3"),
  #   cr2 = c("cr2", "cr4"),
  #   cr5=c("cr5", "cr6"),
  #   cr6=c("cr6", "cr5", "cr7")
  # )
  # unique(sapply(list_to_fuse, function(x) 
  #   unique(unlist(list_to_fuse[sapply(list_to_fuse, function(y) 
  #     any(x %in% y))])), simplify=FALSE))
  regions_to_fuse <- unique(sapply(list_to_fuse, function(x) 
    unique(unlist(list_to_fuse[sapply(list_to_fuse, function(y) 
      any(x %in% y))])), simplify=FALSE))
  
  stopifnot(!duplicated(unlist(regions_to_fuse)))
  
  length_before_merging <- length(set_of_tads)
  stopifnot(length_before_merging == length(set_of_tads_intersect_genes))
  set_of_tads_before_merging <- set_of_tads
  set_of_tads_intersect_genes_before_merging <- set_of_tads_intersect_genes
  
  txt <- paste0("... # of merging to perform:\t", length(regions_to_fuse) , "\n")
  cat(txt)
  cat(txt, file=logFile, append=T)
  
  i=1
  for(i in seq_along(regions_to_fuse)) {
    tomerge <- regions_to_fuse[[i]]
    txt <- paste0("... merging ", i, "- ", paste0(tomerge, collapse=","), "\n")
    cat(txt)
    cat(txt, file=logFile, append=T)
    foo <- sapply(unname(set_of_tads_intersect_genes[names(set_of_tads_intersect_genes) %in% tomerge]), function(x) { 
      txt <- paste0(as.character(x), "\n")
      cat(txt)
      cat(txt, file=logFile, append=T)
    })
    
    new_region <- paste0(tomerge, collapse="_")
    new_entry_genes <- unique(unlist(set_of_tads_intersect_genes[names(set_of_tads_intersect_genes) %in% tomerge]))
    new_entry <- sort(unique(unlist(set_of_tads[names(set_of_tads) %in% tomerge])))
    
    
    set_of_tads <- set_of_tads[! names(set_of_tads) %in% tomerge]
    set_of_tads[[paste0(new_region)]] <- new_entry
    
    set_of_tads_intersect_genes <- set_of_tads_intersect_genes[! names(set_of_tads_intersect_genes) %in% tomerge]
    set_of_tads_intersect_genes[[paste0(new_region)]] <- new_entry_genes
    
  }
  stopifnot(length(set_of_tads) == (length_before_merging - length(unlist(regions_to_fuse)) + length(regions_to_fuse)))
  stopifnot(length(set_of_tads) == length(set_of_tads_intersect_genes))
  stopifnot(setequal(names(set_of_tads), names(set_of_tads_intersect_genes)))
  stopifnot(setequal(unlist(set_of_tads_before_merging), unlist(set_of_tads)))
  stopifnot(setequal(unlist(set_of_tads_intersect_genes_before_merging), unlist(set_of_tads_intersect_genes)))
  
  txt <- paste0("... merge overlapping sets:\t", length(set_of_tads_before_merging) , " -> ",  length(set_of_tads), "\n")
  cat(txt)
  cat(txt, file=logFile, append=T)
  
  
  names(set_of_tads) <- paste0("conserved_region_", seq_along(set_of_tads))
  set_of_tads_intersect_genes <- foreach(cr = names(set_of_tads)) %dopar% {
    all_regions <- set_of_tads[[paste0(cr)]]
    ds_genes <- sapply(all_regions, function(tad) {
      dt <- all_data_list[[paste0(dirname(tad))]][["dataset_g2t_dt"]]
      stopifnot(basename(tad) %in% dt$region)
      as.character(dt$entrezID[dt$region == basename(tad)])
    })
    is_genes <- Reduce(intersect, ds_genes)
    stopifnot(is_genes %in% gff_dt$entrezID)
    stopifnot(is_genes %in% names(entrez2symb))
    setNames(entrez2symb[paste0(is_genes)], is_genes)
  }
  names(set_of_tads_intersect_genes) <- names(set_of_tads)
  
  # just ensure: "conserved region" has at least 1 match
  stopifnot(lengths(set_of_tads) >= 2)
  
  # => for the GO analysis, save conserved regions
  conserved_signif_tads <- set_of_tads
  outFile <- file.path(outFolder,
                       i_perm,
                       paste0(file_prefix, "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
  save(conserved_signif_tads, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  # => final set of conserved regions here !
  stopifnot(!unlist(conserved_signif_tads) %in% true_signif_tads)
  
  #######################################################################################################################################
  ######################################################################################################################### PLOT - TAD occurence
  #######################################################################################################################################                                       
  
  all_tad_vect <- as.character(unlist(set_of_tads))
  tad_occurences <- setNames(as.numeric(table(all_tad_vect)), names(table(all_tad_vect)))
  # outfile <- file.path(outFolder, paste0(file_prefix, "tad_occurences_in_conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".",plotType))
  # do.call(plotType, list(outfile, height=myHeight, width=myHeight*1.2))
  # plot(density(tad_occurences), main=paste0(data_cmpType))
  # mtext(side=3, text="TAD occurence")
  # foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  tail(sort(tad_occurences))
  txt <- paste0("> summary(tad_occurences):\n")
  cat(txt)
  cat(txt, file=logFile, append=T)
  txt <- paste0(names(summary(tad_occurences)))
  cat(txt)
  cat(txt, file=logFile, append=T)
  txt <- paste0(summary(tad_occurences), "\n")
  cat(txt)
  cat(txt, file=logFile, append=T)
  
  #######################################################################################################################################
  ######################################################################################################################### PLOT - overview conservation (heatmap)
  #######################################################################################################################################                                       
  all_datasets <- unique(all_signif_matching_dt$ref_dataset)
  
  
  # BINARY: 1 if the dataset has a TAD belonging to the conserved region
  cons_region = set_of_tads[[15]]
  matching_dt <- foreach(cons_region = set_of_tads, .combine='cbind') %dopar% {
    all_matching_ds <- dirname(cons_region)
    matching_ds <- as.numeric(ds_levels %in% all_matching_ds )
    # stopifnot(sum(matching_ds) == length(all_matching_ds)) # not TRUE -> 1 TAD can match to multiple TAD in the query
    stopifnot(sum(matching_ds) <= length(all_matching_ds)) # not TRUE -> 1 TAD can match to multiple TAD in the query
    matching_ds
  }
  stopifnot(!duplicated(ds_levels))
  rownames(matching_dt) <- ds_levels
  colnames(matching_dt) <- names(set_of_tads)
  
  ### !!! >>> v2 update here -> # of conservations based on exprds, not hicds+exprds
  # levels to sort the columns by decreasing # of matches
  # tmp <- colSums(matching_dt)
  # region_levels <- names(sort(-tmp))
  tmp_dt <- data.frame(matching_dt)
  tmp_dt$hicds <- dirname(rownames(matching_dt))
  tmp_dt$exprds <- basename(rownames(matching_dt))
  tmp_dt_m <- melt(tmp_dt, id=c("hicds", "exprds"))
  # tmp_dt_m2b <- aggregate(value ~ exprds+variable, data=tmp_dt_m, FUN=sum)
  # tmp_dt_m3b <- aggregate(value ~ variable, data=tmp_dt_m2b, FUN=function(x) sum(x>0))
  # tmp_dt_m3b <- tmp_dt_m3b[order(tmp_dt_m3b$value, decreasing=TRUE),]
  # stopifnot(tmp_dt_m3$exprds == tmp_dt_m3b$value)
  # stopifnot(tmp_dt_m3$variable == tmp_dt_m3b$variable)
  tmp_dt_m2 <- tmp_dt_m[tmp_dt_m$value > 0,]
  tmp_dt_m3 <- aggregate(exprds~variable, data=tmp_dt_m2, FUN=function(x) length(unique(x)))
  tmp_dt_m3 <- tmp_dt_m3[order(tmp_dt_m3$exprds, decreasing=TRUE),]
  region_levels <- tmp_dt_m3$variable
  exprds_match_dt <- tmp_dt_m3
  
  out_df <- data.frame(
    conserved_region = names(set_of_tads),
    corresp_tads = as.character(unlist(lapply(set_of_tads, function(x) paste0(x, collapse=",")))),
    intersect_genes_symbol = as.character(unlist(lapply(set_of_tads_intersect_genes, function(x) paste0(x, collapse=",")))),
    intersect_genes_entrez = as.character(unlist(lapply(set_of_tads_intersect_genes, function(x) paste0(names(x), collapse=",")))),
    stringsAsFactors=FALSE
  )
  out_df$conserved_region <- factor(out_df$conserved_region, levels = region_levels)
  stopifnot(!is.na(out_df$conserved_region))
  out_df <- out_df[order(as.numeric(out_df$conserved_region)),]
  outFile <- file.path(outFolder,
                       i_perm,
                       paste0(file_prefix, "conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_df.txt"))
  write.table(out_df, append=F, quote=F, sep="\t", file = outFile, col.names=TRUE, row.names=FALSE)
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, 
                       i_perm,
                       paste0(file_prefix, "conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
  save(out_df, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  #### heatmap version
  plot_matching_dt <- matching_dt
  # rownames(plot_matching_dt) <- gsub("/", "\n", rownames(plot_matching_dt))

    save(plot_matching_dt, file = file.path(outFolder, i_perm,
                                            paste0("plot_matching_dt_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata")), 
       version=2)
  
  
  
  #######################################################################################################################################
  ######################################################################################################################### lolli plot nTop regions by cmp types
  #######################################################################################################################################                   
  
  cmpType = unique(as.character(dataset_dt$cmpType))[1]
  for(cmpType in unique(as.character(dataset_dt$cmpType))) {
    
    ds_to_keep <- dataset_dt$dataset[as.character(dataset_dt$cmpType) == cmpType]
    stopifnot(ds_to_keep %in% rownames(matching_dt))
    
    cmp_matching_dt <- matching_dt[paste0(ds_to_keep),]
    stopifnot(nrow(cmp_matching_dt) == length(ds_to_keep))
    
    #>>> v2 update: region levels defined according to exprds
    # cmp_region_levels <- names(sort(-colSums(cmp_matching_dt)))
    cmp_tmp_dt <- data.frame(cmp_matching_dt)
    cmp_tmp_dt$hicds <- dirname(rownames(cmp_matching_dt))
    cmp_tmp_dt$exprds <- basename(rownames(cmp_matching_dt))
    cmp_tmp_dt_m <- melt(cmp_tmp_dt, id=c("hicds", "exprds"))
    cmp_tmp_dt_m2 <- cmp_tmp_dt_m[cmp_tmp_dt_m$value > 0,]
    cmp_tmp_dt_m3 <- aggregate(exprds~variable, data=cmp_tmp_dt_m2, FUN=function(x) length(unique(x)))
    cmp_tmp_dt_m3 <- cmp_tmp_dt_m3[order(cmp_tmp_dt_m3$exprds, decreasing=TRUE),]
    cmp_region_levels <- cmp_tmp_dt_m3$variable
    cmp_exprds_match_dt <- cmp_tmp_dt_m3
    
    nDS_cmp <- nrow(cmp_matching_dt)
    nExprds_cmp <- length(unique(basename(rownames(cmp_matching_dt))))
    
    
    cmp_out_df <- out_df
    cmp_out_df$conserved_region <- factor(cmp_out_df$conserved_region, levels = cmp_region_levels)
    cmp_out_df <- cmp_out_df[order(as.numeric(cmp_out_df$conserved_region)),]
    outFile <- file.path(outFolder, 
                         i_perm,
                         paste0(file_prefix, "_", cmpType, "_conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_df.txt"))
    write.table(cmp_out_df, append=F, quote=F, sep="\t", file = outFile, col.names=TRUE, row.names=FALSE)
    cat(paste0("... written: ", outFile, "\n"))
  
    
    
  }
  
    
  }
  



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

