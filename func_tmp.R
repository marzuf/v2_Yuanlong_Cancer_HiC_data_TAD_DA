
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

require(foreach)
require(doMC)
require(GenomicRanges)

nCpu = 40

logFile=""

### > Filter1: retain only matches with >= *minOverlapBpRatio* (80%) bp overlap
minOverlapBpRatio = 0.8


### > Filter2: retain only "conserved regions" with >= *minIntersectGenes*(3) genes at the intersect
minIntersectGenes = 3


### !!! >>> v2 update here 
### Merge conserved regions that have >= *gene_matching_fuse_threshold* % (80%) gene overlap
gene_matching_fuse_threshold = 0.8

registerDoMC(nCpu)

script0_name <- "0_prepGeneData"

outFolder <- "TMP_FUNC_CONS_oneshot"
dir.create(outFolder)

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
entrez2symb <- setNames(gff_dt$symbol,gff_dt$entrezID)

# for testing
mainFolder <- file.path(".")
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("PERMUT", all_hicds) & !grepl("RANDOM", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

# for testing - prep the signif. tables
final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))
final_dt$dataset <- file.path(final_dt$hicds,final_dt$exprds) 
signif_dt <- final_dt[final_dt$adjPvalComb <= 0.01,c("dataset", "region")]
save(signif_dt, file=file.path(outFolder, "signif_dt.Rdata"))
load(file.path(outFolder, "signif_dt.Rdata"))

# for testing, prep the coordinates of TADs
all_tad_pos_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
  hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
  stopifnot(file.exists(hicds_file))
  tadpos_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
  tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]
  stopifnot(nrow(tadpos_dt) > 0 )
  tadpos_dt$dataset <- file.path(hicds, exprds)
  tadpos_dt <- tadpos_dt[, c("dataset", "chromo", "region", "start", "end")]
  # take only the signif ones
  signif_tads <- signif_dt$region[signif_dt$dataset == file.path(hicds, exprds)]
  tadpos_dt[tadpos_dt$region %in% signif_tads,]
  }
  exprds_dt
}
save(all_tad_pos_dt, file=file.path(outFolder, "all_tad_pos_dt.Rdata"))
load(file.path(outFolder, "all_tad_pos_dt.Rdata"))

# for testing, prep the g-2-t assignment TADs
# for testing, prep the coordinates of TADs
all_g2t_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
    stopifnot(geneList %in% g2t_dt$entrezID)
    g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
    # take only the signif ones
    signif_tads <- signif_dt$region[signif_dt$dataset == file.path(hicds, exprds)]
    g2t_dt$dataset <- file.path(hicds, exprds)
    g2t_dt <- g2t_dt[, c("dataset", "entrezID", "chromo", "start", "end", "region")]
    g2t_dt[g2t_dt$region %in% signif_tads,]
  }
  exprds_dt
}
all_g2t_dt$symbol <- entrez2symb[all_g2t_dt$entrezID]
save(all_g2t_dt, file=file.path(outFolder, "all_g2t_dt.Rdata"))
load(file.path(outFolder, "all_g2t_dt.Rdata"))

# stop("--ok\n")
    
# > head(signif_dt)
# dataset       region
# 1514 Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr6_TAD633
# 121  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr1_TAD551
# 302  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas chr11_TAD238
# 344  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas chr11_TAD384
# 692  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas   chr16_TAD6
# 932  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas chr19_TAD244
# > head(all_tad_pos_dt)
# dataset chromo       region     start       end
# 68   Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr10  chr10_TAD67  14880001  15040000
# 821  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr11 chr11_TAD238  61920001  62160000
# 971  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr11 chr11_TAD384  95720001  96160000
# 996  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr11 chr11_TAD409 102160001 102560000
# 1006 Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr11 chr11_TAD419 104680001 105080000
# 2364 Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas  chr14 chr14_TAD227  74360001  74560000
# > head(all_g2t_dt)
# dataset entrezID chromo    start      end       region
# 147  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas    51182  chr10 14880159 14913740  chr10_TAD67
# 149  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas    79723  chr10 14920782 14946314  chr10_TAD67
# 150  Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas    64421  chr10 14948870 14996106  chr10_TAD67
# 2671 Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas    10647  chr11 62009102 62012280 chr11_TAD238
# 2672 Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas     4250  chr11 62037630 62040628 chr11_TAD238
# 2675 Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas    80150  chr11 62104774 62160887 chr11_TAD238

all_datasets <- unique(signif_dt$dataset)
ref_dataset <- all_datasets[1]
# do the all vs. all matching

signif_dt$region <- as.character(signif_dt$region)
signif_dt$dataset <- as.character(signif_dt$dataset)
all_g2t_dt$region <- as.character(all_g2t_dt$region)
all_g2t_dt$dataset <- as.character(all_g2t_dt$dataset)
all_tad_pos_dt$region <- as.character(all_tad_pos_dt$region)
all_tad_pos_dt$dataset <- as.character(all_tad_pos_dt$dataset)

# => I am only interested in the TADs that are passed in the signif tables, 
# so filter the other tables in case they have extra information
signif_dt$regID <- file.path(signif_dt$dataset, signif_dt$region)
stopifnot(!duplicated(signif_dt$regID))

all_g2t_dt$regID <- file.path(all_g2t_dt$dataset, all_g2t_dt$region)
stopifnot(!duplicated(file.path(all_g2t_dt$entrezID, all_g2t_dt$regID)))
stopifnot(signif_dt$regID %in% all_g2t_dt$regID)
all_g2t_dt <- all_g2t_dt[all_g2t_dt$regID %in% signif_dt$regID,]
stopifnot(nrow(all_g2t_dt) > 0)

all_tad_pos_dt$regID <- file.path(all_tad_pos_dt$dataset, all_tad_pos_dt$region)
stopifnot(signif_dt$regID %in% all_tad_pos_dt$regID)
stopifnot(!duplicated(all_tad_pos_dt$regID))
all_tad_pos_dt <- all_tad_pos_dt[all_tad_pos_dt$regID %in% signif_dt$regID,]
stopifnot(nrow(all_tad_pos_dt) > 0)

all_tad_pos_dt$ref_totBp <- all_tad_pos_dt$end - all_tad_pos_dt$start + 1

tad_size <- setNames(all_tad_pos_dt$ref_totBp, all_tad_pos_dt$regID)


ngenes_dt <- aggregate(entrezID~regID, FUN=length, data=all_g2t_dt)

tad_ngenes <- setNames(ngenes_dt$entrezID,ngenes_dt$regID)

stopifnot(length(tad_size) == length(tad_ngenes))
stopifnot(setequal(names(tad_size), names(tad_ngenes)))

### all vs all matching -> match with itself
    ### PREPARE THE BP OVERLAP
    cat("... preparing bp matching overlap\n")
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

overlapDT = overlapDT[order(overlapDT$ref_dataset, overlapDT$refID, overlapDT$query_dataset, overlapDT$queryID),]
rownames(overlapDT) =NULL
all_signif_matching_dt=overlapDT



outFile <- file.path(outFolder, paste0("all_signif_matching_dt.Rdata"))
save(all_signif_matching_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
  stop("---ok")

 # si je fais ça au début -> pas besoin de faire le match query_dataset by query_dataset ??
all_signif_matching_dt$refID_full <- file.path(all_signif_matching_dt$ref_dataset, all_signif_matching_dt$refID)
all_signif_matching_dt$queryID_full <- file.path(all_signif_matching_dt$query_dataset, all_signif_matching_dt$queryID)




























































all_signif_matching_dt <- foreach(ref_dataset = all_datasets, .combine='rbind') %dopar% {
  

  # retrieve list of genes for this dataset
  ref_geneList <- all_g2t_dt$entrezID[all_g2t_dt$dataset == ref_dataset]
  ref_g2t_dt <- all_g2t_dt[all_g2t_dt$dataset == ref_dataset,]
  stopifnot(!duplicated(ref_g2t_dt$entrezID))
  ref_tadpos_dt <- all_tad_pos_dt[all_tad_pos_dt$dataset == ref_dataset,]
  
  stopifnot(setequal(ref_g2t_dt$region, ref_tadpos_dt$region))
  ref_tads <- signif_dt$region[signif_dt$dataset == ref_dataset]
  stopifnot(setequal(as.character(ref_tadpos_dt$region), ref_tads))
  
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
  
  
  ### !!!! > DO I REALLY NEED TO DO THIS QUERY_DATASET BY QUERY_DATASET ??? ###
  
  
  ref_dt <- foreach(query_dataset = query_datasets, .combine='rbind') %do% {
    cat("> Start ", ref_dataset, " vs. ", query_dataset, "\n")  
    
    query_geneList <- all_g2t_dt$entrezID[all_g2t_dt$dataset == query_dataset]
    query_g2t_dt <- all_g2t_dt[all_g2t_dt$dataset == query_dataset,]
    stopifnot(!duplicated(query_g2t_dt$entrezID))
    query_tadpos_dt <- all_tad_pos_dt[all_tad_pos_dt$dataset == query_dataset,]
    
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

    
    # ensure only same chromo are compared      
    overlapDT_bp$chromo_ref <- gsub("_TAD.+", "", as.character(overlapDT_bp$refID))
    overlapDT_bp$chromo_query <- gsub("_TAD.+", "", as.character(overlapDT_bp$queryID))
    stopifnot(overlapDT_bp$chromo_ref == overlapDT_bp$chromo_query)
    overlapDT_bp$chromo_ref <- NULL
    overlapDT_bp$chromo_query <- NULL
    
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
      # ref_genes <- final_dt$region_genes[final_dt$hicds == dirname(ref_dataset) &
      #                                      final_dt$exprds == basename(ref_dataset) &
      #                                      final_dt$region == refID[i]]
      # ref_genes_ul <- unlist(strsplit(ref_genes, split=","))
      # query_genes <- final_dt$region_genes[final_dt$hicds == dirname(query_dataset) &
      #                                        final_dt$exprds == basename(query_dataset) &
      #                                        final_dt$region == queryID[i]]
      # query_genes_ul <- unlist(strsplit(query_genes, split=","))
      
      ref_genes_ul <- all_g2t_dt$symbol[all_g2t_dt$dataset == ref_dataset & all_g2t_dt$region == refID[i] ]
      query_genes_ul <- all_g2t_dt$symbol[all_g2t_dt$dataset == query_dataset & all_g2t_dt$region == queryID[i] ]
      
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

outFile <- file.path(outFolder, paste0("all_signif_matching_dt.Rdata"))
save(all_signif_matching_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
  

 # si je fais ça au début -> pas besoin de faire le match query_dataset by query_dataset ??
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
nIntersectGenesByRegions <- lengths(set_of_tads_intersect_genes)
regionsWithMinGenes <- names(nIntersectGenesByRegions) [nIntersectGenesByRegions >= minIntersectGenes]
stopifnot(regionsWithMinGenes %in% names(set_of_tads))
all_set_of_tads  <- set_of_tads
set_of_tads <- set_of_tads[paste0(regionsWithMinGenes)]
stopifnot(setequal(regionsWithMinGenes, names(set_of_tads)))

### !!! >>> v2 update here 
### Merge conserved regions that have >= *gene_matching_fuse_threshold* % (80%) gene overlap

txt <- paste0("> Merging if gene overlap - # of conserved regions\t=\t", length(set_of_tads), "\n")
cat(txt)
cat(txt, file=logFile, append=T)

set_of_tads_intersect_genes <- foreach(cr = names(set_of_tads)) %dopar% {
  all_regions <- set_of_tads[[paste0(cr)]]
  ds_genes <- sapply(all_regions, function(tad) {
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
outFile <- file.path(outFolder, paste0(file_prefix, "conserved_signif_tad_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
save(conserved_signif_tads, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

# => final set of conserved regions here !

stop("STOP HERE\n")

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
outFile <- file.path(outFolder, paste0(file_prefix, "conserved_regions_with_genes_signif_tads_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_df.txt"))
write.table(out_df, append=F, quote=F, sep="\t", file = outFile, col.names=TRUE, row.names=FALSE)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0(file_prefix, "conserved_regions_with_genes_signif_tads_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
save(out_df, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


#### heatmap version
plot_matching_dt <- matching_dt
# rownames(plot_matching_dt) <- gsub("/", "\n", rownames(plot_matching_dt))
initmar <- par()$mar

# modifMar <- c(0,-4.1,-4.1,15)
modifMar <- c(0,0,0,0)

outFile <- file.path(outFolder, paste0(file_prefix, "heatmap_conservedRegions_rowReorder_signif_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".", plotType))
do.call(plotType, list(outFile,  height=myHeightHeat, width=myWidthHeat))
par(mar=(initmar+ modifMar))
# dev.off()
my_heatmap.2(plot_matching_dt[,paste0(region_levels)], col=c("white", "black"),
          dendrogram = "none",
          trace = "none",
          key=FALSE,
          Rowv=TRUE,
          Colv=FALSE,
          labCol = FALSE,
          lmat = rbind(c(5,0,4), c(3,1,2)),
          # lwid=c(1.5,0.2,4),
          lwid=c(0.2,0.2,5),
          lhei=c(0.01,5),
          margins=c(2,18),
          adjRow=c(0,0.5),
          #lwid= c(0.1, 0.1, 5),
          colRow = dataset_dt$subtype_col,
          RowSideColors = dataset_dt$subtype_col,
          rowAxis=4
          )


# default
# lmat: 
#   5 3 0 1 4 2 
# lhei: 
#   1.5 4 
# lwid: 
#   1.5 0.2 4 

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

save(plot_matching_dt, file = file.path(outFolder, paste0("plot_matching_dt_signif_tads_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata")), 
     version=2)

par(mar = initmar)

outFile <- file.path(outFolder, paste0(file_prefix,"heatmap_conservedRegions_noReorder_signif_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".", plotType))
do.call(plotType, list(outFile,  height=myHeightHeat, width=myWidthHeat))
par(mar=initmar+modifMar)

my_heatmap.2(plot_matching_dt[,paste0(region_levels)], col=c("white", "black"),
          dendrogram = "none",
          trace = "none",
          key=F,
          Rowv=FALSE,
          Colv=FALSE,
          labCol = FALSE,
          # cexRow = 0.4,
          rowAxis=2,
          lmat = rbind(c(4,3), c(2,1)),
           lwid=c(1.5,4),
           lhei=c(0.01,1),
           margins=c(5,5),
          adjRow=c(1,0.5),
          # RowSideColors = dataset_dt$subtype_col,
          colRow = dataset_dt$subtype_col
          )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




