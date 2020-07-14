

# x1=get(load("TMP_FUNC_CONS_ONESHOT/conserved_signif_tad_minBpRatio0.8_minInterGenes3.Rdata"))
# x2=get(load("TMP_FUNC_CONS_v0//conserved_signif_tad_minBpRatio0.8_minInterGenes3.Rdata"))
# signif_conserv_data= get(load("TMP_FUNC_CONS_vFUNC/signif_conserv_data.Rdata"))
# x3 =  signif_conserv_data[["conserved_signif_tads"]]
# all.equal(x1,x2)
# all.equal(x1,x3)
# 
# x0 = get(load("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2/conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
# names(x0) <- NULL
# x0t <- unlist(lapply(x0, function(x) paste0(sort(x), collapse=",")))
# names(x1) <- NULL
# x1t <- unlist(lapply(x1, function(x) paste0(sort(x), collapse=",")))
# all.equal(sort(x0t),sort(x1t))
# 
# # 
# x1=get(load("TMP_FUNC_CONS_ONESHOT/all_signif_matching_dt2.Rdata"))
# x2=get(load("TMP_FUNC_CONS_v0//all_signif_matching_dt2.Rdata"))
# all.equal(x1,x2)
# signif_conserv_data=get(load("TMP_FUNC_CONS_vFUNC/signif_conserv_data.Rdata"))
# x3 =  signif_conserv_data[["matching_table"]]
# all.equal(x2,x3)
# all.equal(x1,x3)

# 
# x1=get(load("TMP_FUNC_CONS_ONESHOT/all_signif_matching_dt2.Rdata"))
# x1 <- x1[order(as.character(x1$refID_full), as.character(x1$queryID_full), as.character(x1$overlapBp) , as.character(x1$nOverlapGenes)),]
# x2=get(load("TMP_FUNC_CONS_v0//all_signif_matching_dt2.Rdata"))
# x2 <- x2[order(as.character(x2$refID_full), as.character(x2$queryID_full), x2$overlapBp , x2$nOverlapGenes),]
# all.equal(x1,x2)
# x3 =  signif_conserv_data[["matching_table"]]
# x3 <- x3[order(as.character(x3$refID_full), as.character(x3$queryID_full), x3$overlapBp , x3$nOverlapGenes),]
# # x1 = x1[order(x1$ref_dataset, x1$query_dataset, x1$refID, x1$queryID),]
# # x2 = x2[order(x2$ref_dataset, x2$query_dataset, x2$refID, x2$queryID),]
# # x3 = x3[order(x3$ref_dataset, x3$query_dataset, x3$refID, x3$queryID),]
# rownames(x1) <- NULL
# rownames(x2) <- NULL
# rownames(x3) <- NULL
# all.equal(x2,x3)
# all.equal(x1,x3)
# 
# 
# c1=get(load("TMP_FUNC_CONS_v0//overlapDT.Rdata"))
# c1 = c1[order(c1$ref_dataset, c1$refID, c1$query_dataset, c1$queryID),]
# head(c1)
# 
# c2=get(load("overlapDT.Rdata"))
# head(c2)






# INPUT:
#   - table with signif. TADs:
#   dataset / region / 
#   
#   - coordinates of the TADs:
#   dataset / region / chromo / start / end 
# 
#  - gene-to-TAD dt [should be processed to contain only the pip genes ! ]
#     dataset / entrezID / chromo / start / end / region / symbol

# column order does not matter

# Rscript FUNC_tad_matching_signif_across_hicds_allMatch_v2.R  



### A VOIR; MAIS JE CROIS QUE TOUTE LA PARTIE AVEC ngenes_dt NE SERT A RIEN

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Retrieve consensus regions for a list of TADs.
#'
#' Retrieve region corresponding to overlapping TADs.
#'
#' @param signif_dt Dataframe with the list of TAds that have to be matched (columns: dataset/region).
#' @param all_tad_pos_dt Coordinates of the TADs hold in signif_dt (colums:  dataset/chromo/region/start/end).
#' @param all_g2t_dt Gene-to-TAD assignment for the TADs hold in signif_dt (columns: dataset/entrezID/chromo/start/end/region/symbol).
#' @param minOverlapBpRatio Retain only matches with >= minOverlapBpRatio bp overlap (filter 1).
#' @param minIntersectGenes Retain only conserved regions with >= minIntersectGenes genes at the intersect (filter 2).
#' @param gene_matching_fuse_threshold Merge conserved regions that have >= gene_matching_fuse_threshold % gene overlap (after filter 2).
#' @param nCpu Number of CPU available.
#' @param logFile If provided, write logs in this file.
#' @param verbose Verbose or not function execution.
#' @return A list of four elements: the all-vs-all matching table ("matching_table"), the conserved regions with the corresponding TADs ("conserved_signif_tads"), the conserved regions with the corresponding intersect genes ("conserved_signif_intersect_genes"), the list of parameters used for the filters  ("parameters").
#' @export
#' 

get_conservedRegion <- function(
  signif_dt, all_tad_pos_dt, all_g2t_dt,
  ### > Filter1: retain only matches with >= *minOverlapBpRatio* (80%) bp overlap
  minOverlapBpRatio = 0.8,
  ### > Filter2: retain only "conserved regions" with >= *minIntersectGenes*(3) genes at the intersect
  minIntersectGenes = 3,
  ### !!! >>> v2 update here 
  ### Merge conserved regions that have >= *gene_matching_fuse_threshold* % (80%) gene overlap
  gene_matching_fuse_threshold = 0.8,
  nCpu = 2,
  logFile=NULL,
  verbose=TRUE
){
  
  if(!suppressPackageStartupMessages(require("foreach"))) stop("-- foreach package required\n")  
  if(!suppressPackageStartupMessages(require("doMC"))) stop("-- doMC package required\n")  
  if(!suppressPackageStartupMessages(require("GenomicRanges"))) stop("-- GenomicRanges package required\n")  
  registerDoMC(nCpu)
  
  all_datasets <- unique(signif_dt$dataset)
  
  # ensure there are all charaters and not vector
  signif_dt$region <- as.character(signif_dt$region)
  signif_dt$dataset <- as.character(signif_dt$dataset)
  all_g2t_dt$region <- as.character(all_g2t_dt$region)
  all_g2t_dt$dataset <- as.character(all_g2t_dt$dataset)
  all_tad_pos_dt$region <- as.character(all_tad_pos_dt$region)
  all_tad_pos_dt$dataset <- as.character(all_tad_pos_dt$dataset)
  
  # for all,  use as ID dataset + TAD label -> this should be unique
  signif_dt$regID <- file.path(signif_dt$dataset, signif_dt$region)
  stopifnot(!duplicated(signif_dt$regID))
  all_g2t_dt$regID <- file.path(all_g2t_dt$dataset, all_g2t_dt$region)
  stopifnot(!duplicated(file.path(all_g2t_dt$entrezID, all_g2t_dt$regID)))
  stopifnot(signif_dt$regID %in% all_g2t_dt$regID)
  all_tad_pos_dt$regID <- file.path(all_tad_pos_dt$dataset, all_tad_pos_dt$region)
  stopifnot(signif_dt$regID %in% all_tad_pos_dt$regID)
  stopifnot(!duplicated(all_tad_pos_dt$regID))
  
  # ensure I have all information for the signif regions
  stopifnot(signif_dt$regID %in% all_g2t_dt$regID)
  stopifnot(signif_dt$regID %in% all_tad_pos_dt$regID)
  
  # => I am only interested in the TADs that are passed in the signif tables, 
  # so filter the other tables in case they have extra information
  all_g2t_dt <- all_g2t_dt[all_g2t_dt$regID %in% signif_dt$regID,]
  stopifnot(nrow(all_g2t_dt) > 0)
  
  all_tad_pos_dt <- all_tad_pos_dt[all_tad_pos_dt$regID %in% signif_dt$regID,]
  stopifnot(nrow(all_tad_pos_dt) > 0)
  
  all_tad_pos_dt$ref_totBp <- all_tad_pos_dt$end - all_tad_pos_dt$start + 1
  
  tad_size <- setNames(all_tad_pos_dt$ref_totBp, all_tad_pos_dt$regID)
  
  
  # ngenes_dt <- aggregate(entrezID~regID, FUN=length, data=all_g2t_dt)
  # tad_ngenes <- setNames(ngenes_dt$entrezID,ngenes_dt$regID)
  # stopifnot(length(tad_size) == length(tad_ngenes))
  # stopifnot(setequal(names(tad_size), names(tad_ngenes)))
  
  ### all vs all matching -> match with itself
  ### PREPARE THE BP OVERLAP
  
  txt <- paste0("... start bp matching overlap\n")
  outTxt(txt, verbose,logFile)
  
  ref_GR <- GRanges(seqnames=all_tad_pos_dt$chromo, ranges=IRanges(start=all_tad_pos_dt$start, end=all_tad_pos_dt$end, names=all_tad_pos_dt$regID)) # change here the names to be regID instead of region and all_tad_pos
  query_GR <- GRanges(seqnames=all_tad_pos_dt$chromo, ranges=IRanges(start=all_tad_pos_dt$start, end=all_tad_pos_dt$end, names=all_tad_pos_dt$regID)) # change here the names to be regID instead of region
  
  # determine which features from the query overlap which features in the subject
  overlap_GR <- findOverlaps(query=query_GR, subject=ref_GR)
  
  if(length(overlap_GR) == 0) return(NULL)
  
  IDoverlap_hits_all <- findOverlaps(query=query_GR,
                                     subject=query_GR)
  IDoverlap_hits <- IDoverlap_hits_all[queryHits(IDoverlap_hits_all) != subjectHits(IDoverlap_hits_all)]
  
  refID <- names(ref_GR[subjectHits(overlap_GR)])
  queryID <- names(query_GR[queryHits(overlap_GR)])
  
  stopifnot(refID %in% names(tad_size)) # change here only tad_size not ref_tad_size
  
  # do not take matching with itself
  toTake <- refID != queryID
  refID <- refID[toTake]
  queryID <- queryID[toTake]
  
  overlapDT_bp <- data.frame(
    refID = refID,
    queryID = queryID,
    overlapBp = width(pintersect(ref_GR[refID], query_GR[queryID])),
    overlapBpRatio = width(pintersect(ref_GR[refID], query_GR[queryID]))/tad_size[refID],
    stringsAsFactors = FALSE)
  
  # ensure only same chromo are compared      
  overlapDT_bp$chromo_ref <- gsub("_TAD.+", "", basename(as.character(overlapDT_bp$refID))) # here take basename only
  overlapDT_bp$chromo_query <- gsub("_TAD.+", "", basename(as.character(overlapDT_bp$queryID)))
  stopifnot(overlapDT_bp$chromo_ref == overlapDT_bp$chromo_query)
  overlapDT_bp$chromo_ref <- NULL
  overlapDT_bp$chromo_query <- NULL
  
  stopifnot(overlapDT_bp$refID %in% all_tad_pos_dt$regID) # no ref_ and match regID not region
  stopifnot(overlapDT_bp$queryID %in% all_tad_pos_dt$regID) # no query_
  stopifnot(!is.na(overlapDT_bp))
  stopifnot(overlapDT_bp$overlapBpRatio >= 0)
  stopifnot(overlapDT_bp$overlapBpRatio <= 1)
  nmat <- length(refID)
  stopifnot(nmat == length(queryID))
  i=1
  i=2
  overlapDT_genes <- foreach(i=1:nmat, .combine='rbind') %do% {
    # ref_genes <- final_dt$region_genes[final_dt$hicds == dirname(ref_dataset) &
    #                                      final_dt$exprds == basename(ref_dataset) &
    #                                      final_dt$region == refID[i]]
    # ref_genes_ul <- unlist(strsplit(ref_genes, split=","))
    # query_genes <- final_dt$region_genes[final_dt$hicds == dirname(query_dataset) &
    #                                        final_dt$exprds == basename(query_dataset) &
    #                                        final_dt$region == queryID[i]]
    # query_genes_ul <- unlist(strsplit(query_genes, split=","))
    
    ref_genes_ul <- all_g2t_dt$symbol[all_g2t_dt$regID == refID[i] ]
    query_genes_ul <- all_g2t_dt$symbol[all_g2t_dt$regID == queryID[i] ]
    #      ref_genes_ul <- all_g2t_dt$symbol[all_g2t_dt$dataset == ref_dataset & all_g2t_dt$region == refID[i] ]
    #      query_genes_ul <- all_g2t_dt$symbol[all_g2t_dt$dataset == query_dataset & all_g2t_dt$region == queryID[i] ]
    
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
  #    overlapDT$ref_dataset <- ref_dataset
  #    overlapDT$query_dataset <- query_dataset
  overlapDT$ref_dataset <- dirname(overlapDT$refID)
  overlapDT$query_dataset <- dirname(overlapDT$queryID)
  overlapDT$refID <- basename(overlapDT$refID)
  overlapDT$queryID <- basename(overlapDT$queryID)
  overlapDT=overlapDT[, c("ref_dataset", "query_dataset", all_cols)]
  #  } # end iterating over query_datasets
  #  ref_dt
  #} # end iterating over ref_datasets
  
  overlapDT <- overlapDT[order(overlapDT$ref_dataset, overlapDT$refID, overlapDT$query_dataset, overlapDT$queryID),]
  rownames(overlapDT) <- NULL
  all_signif_matching_dt <- overlapDT
  # outFile <- file.path(outFolder, paste0("all_signif_matching_dt.Rdata"))
  # save(all_signif_matching_dt, file=outFile, version=2)
  # cat(paste0("... written: ", outFile, "\n"))
  # stop("---ok")
  
  # si je fais ça au début -> pas besoin de faire le match query_dataset by query_dataset ??
  all_signif_matching_dt$refID_full <- file.path(all_signif_matching_dt$ref_dataset, all_signif_matching_dt$refID)
  all_signif_matching_dt$queryID_full <- file.path(all_signif_matching_dt$query_dataset, all_signif_matching_dt$queryID)
  
  
  all_signif_matching_dt <- all_signif_matching_dt[order(as.character(all_signif_matching_dt$refID_full), 
                                                         as.character(all_signif_matching_dt$queryID_full), 
                                                         as.numeric(as.character(all_signif_matching_dt$overlapBp)), 
                                                         as.numeric(as.character(all_signif_matching_dt$nOverlapGenes))),]
  rownames(all_signif_matching_dt) <- NULL
  all_vs_all_matching_dt <- all_signif_matching_dt

save(all_vs_all_matching_dt, file="all_vs_all_matching_dt.Rdata", version=2)

  
  #######################################################################################################
  ####################################################################################################### CONSERVED REGIONS REFINEMENT
  #######################################################################################################
  
  ### > Filter1: retain only matches with >= *minOverlapBpRatio* (80%) bp overlap
  
  txt <- paste0("> FILTER 1 - # of match\t=\t", nrow(all_signif_matching_dt), "\n")
  outTxt(txt, verbose,logFile)
  
  all_signif_matching_dt  <- all_signif_matching_dt[all_signif_matching_dt$overlapBpRatio >= minOverlapBpRatio,]
  txt <- paste0("> FILTER 1 - # of match with overlap bp >= ", minOverlapBpRatio, "\t=\t", nrow(all_signif_matching_dt), "\n")
  outTxt(txt, verbose,logFile)
  
  all_signif_tads_with_match <- unique(all_signif_matching_dt$refID_full)
  
  tad=all_signif_tads_with_match[1]
  set_of_tads <- foreach(tad = all_signif_tads_with_match) %dopar% {
    # all_signif_matching_dt[all_signif_matching_dt$refID_full == tad,]
    matching_tads <- all_signif_matching_dt$queryID_full[all_signif_matching_dt$refID_full == tad]
    sort(c(tad, matching_tads))
  }
  names(set_of_tads) <- all_signif_tads_with_match
  
  ### Remove duplicated sets
  
  txt <- paste0("... remove duplicated sets:\t", length(set_of_tads) , " -> ")
  outTxt(txt, verbose,logFile)
  
  set_of_tads <- unique(set_of_tads) # loose the name !
  txt <- paste0(length(set_of_tads) , "\n")
  outTxt(txt, verbose,logFile)
  
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
  outTxt(txt, verbose,logFile)
  
  set_of_tads <- not_nested_set_of_tads
  txt <- paste0(length(set_of_tads) , "\n")
  outTxt(txt, verbose,logFile)
  
  
  ### > Filter2: retain only "conserved regions" with >= *minIntersectGenes*(3) genes at the intersect
  txt <- paste0("> FILTER 2 - # of conserved regions\t=\t", length(set_of_tads), "\n")
  outTxt(txt, verbose,logFile)
  
save(set_of_tads, file="set_of_tads.Rdata", version=2)

  names(set_of_tads) <- paste0("conserved_region_", seq_along(set_of_tads))
  cr=names(set_of_tads)[1]
  set_of_tads_intersect_genes <- foreach(cr = names(set_of_tads)) %dopar% { # look at genes at the intersect to filter those with < minIntersectGenes at the intersect  
    all_regions <- set_of_tads[[paste0(cr)]]
    ds_genes <- sapply(all_regions, function(tad) {
      dt <- all_g2t_dt[all_g2t_dt$dataset == paste0(dirname(tad)),]
      stopifnot(basename(tad) %in% dt$region)
      setNames(as.character(dt$symbol[dt$region == basename(tad)]), as.character(dt$entrezID[dt$region == basename(tad)]))
    })
    Reduce(intersect, ds_genes)
  }
save(set_of_tads_intersect_genes, file="set_of_tads_intersect_genes.Rdata", version=2)

  names(set_of_tads_intersect_genes) <- names(set_of_tads)
  nIntersectGenesByRegions <- lengths(set_of_tads_intersect_genes)
  regionsWithMinGenes <- names(nIntersectGenesByRegions) [nIntersectGenesByRegions >= minIntersectGenes]
  stopifnot(regionsWithMinGenes %in% names(set_of_tads))
  # all_set_of_tads  <- set_of_tads
  set_of_tads <- set_of_tads[paste0(regionsWithMinGenes)]  ## set_of_tads updated here !
  stopifnot(setequal(regionsWithMinGenes, names(set_of_tads)))
  
  ### !!! >>> v2 update here 
  ### Merge conserved regions that have >= *gene_matching_fuse_threshold* % (80%) gene overlap
  
  txt <- paste0("> Merging if gene overlap - # of conserved regions\t=\t", length(set_of_tads), "\n")
  outTxt(txt, verbose,logFile)
  
  set_of_tads_intersect_genes <- foreach(cr = names(set_of_tads)) %dopar% {  # same as before to retrieve the genes at the intersect (NB: set_of_tads has been updated in the meantime !)
    all_regions <- set_of_tads[[paste0(cr)]]
    ds_genes <- sapply(all_regions, function(tad) {
      dt <- all_g2t_dt[all_g2t_dt$dataset == paste0(dirname(tad)),]
      stopifnot(basename(tad) %in% dt$region)
      setNames(as.character(dt$symbol[dt$region == basename(tad)]), as.character(dt$entrezID[dt$region == basename(tad)]))
    })
    Reduce(intersect, ds_genes)
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
  outTxt(txt, verbose,logFile)
  
  i=1
  for(i in seq_along(regions_to_fuse)) {
    tomerge <- regions_to_fuse[[i]]
    txt <- paste0("... merging ", i, "- ", paste0(tomerge, collapse=","), "\n")
    outTxt(txt, verbose,logFile)
    
    foo <- sapply(unname(set_of_tads_intersect_genes[names(set_of_tads_intersect_genes) %in% tomerge]), function(x) { 
       txt <- paste0(as.character(x), "\n")
       outTxt(txt, verbose,logFile)
       
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
  outTxt(txt, verbose,logFile)
  
  
  names(set_of_tads) <- paste0("conserved_region_", seq_along(set_of_tads))
  set_of_tads_intersect_genes <- foreach(cr = names(set_of_tads)) %dopar% {
    all_regions <- set_of_tads[[paste0(cr)]]
    ds_genes <- sapply(all_regions, function(tad) {
      # dt <- all_data_list[[paste0(dirname(tad))]][["dataset_g2t_dt"]]
      dt <- all_g2t_dt[all_g2t_dt$dataset == paste0(dirname(tad)),]
      stopifnot(basename(tad) %in% dt$region)
      setNames(as.character(dt$symbol[dt$region == basename(tad)]), as.character(dt$entrezID[dt$region == basename(tad)]))
    })
    # is_genes <- Reduce(intersect, ds_genes)
    # stopifnot(is_genes %in% gff_dt$entrezID)
    # stopifnot(is_genes %in% names(entrez2symb))
    # setNames(entrez2symb[paste0(is_genes)], is_genes)
    Reduce(intersect, ds_genes)
  }
  names(set_of_tads_intersect_genes) <- names(set_of_tads)
  
  # just ensure: "conserved region" has at least 1 match
  stopifnot(lengths(set_of_tads) >= 2)
  
  # => for the GO analysis, save conserved regions
  conserved_signif_tads <- set_of_tads
  # outFile <- file.path(outFolder, paste0("conserved_signif_tad_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
  # save(conserved_signif_tads, file=outFile, version=2)
  # cat(paste0("... written: ", outFile, "\n"))
  
  # => final set of conserved regions here !
  
  # stop("STOP HERE\n")
  
  stopifnot(setequal(names(conserved_signif_tads), names(set_of_tads_intersect_genes)))
  
  return(list(
    matching_table=all_vs_all_matching_dt,
    conserved_signif_intersect_genes=set_of_tads_intersect_genes,
    conserved_signif_tads=conserved_signif_tads,
    parameters=setNames(c(minIntersectGenes, gene_matching_fuse_threshold, minOverlapBpRatio), c("minIntersectGenes","gene_matching_fuse_threshold","minOverlapBpRatio"))
    ))
  
}



