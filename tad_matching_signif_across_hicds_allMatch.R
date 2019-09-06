options(scipen=100)

SSHFS=F

# Rscript tad_matching_across_hicds_allMatch.R

script_name <- "tad_matching_signif_across_hicds_allMatch.R"

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

buildTable <- FALSE

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightHeat <- myHeight * 1.8
myWidthHeat <- myWidth * 1.8

myWidthGG <- 12
myHeightGG <- 12


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

minOverlapBpRatio <- 0.8
minIntersectGenes <- 3

nRegionLolli <- 10

cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> signifThresh\t=\t", signifThresh, "\n"))
cat(paste0("> minOverlapBpRatio\t=\t", minOverlapBpRatio, "\n"))
cat(paste0("> minIntersectGenes\t=\t", minIntersectGenes, "\n"))
cat(paste0("> nRegionLolli\t=\t", nRegionLolli, "\n"))

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
  
  outFile <- file.path(outFolder, paste0("all_signif_matching_dt_", signifcol, ".Rdata"))
  save(all_signif_matching_dt, file=outFile, version=2)  
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, paste0("all_signif_matching_dt_", signifcol, ".Rdata"))
  cat("... load data\n")
  all_signif_matching_dt <- get(load(outFile))
  
  all_data_list <- get(load(file.path(outFolder, "all_data_list.Rdata")))
  
}

# nrow(all_signif_matching_dt) # 12132
# length(unique(paste0(all_signif_matching_dt$ref_dataset, all_signif_matching_dt$refID))) # 1303
# length(unique(paste0(all_signif_matching_dt$query_dataset, all_signif_matching_dt$queryID))) # 1303

all_signif_matching_dt$refID_full <- file.path(all_signif_matching_dt$ref_dataset, all_signif_matching_dt$refID)
all_signif_matching_dt$queryID_full <- file.path(all_signif_matching_dt$query_dataset, all_signif_matching_dt$queryID)

cat(paste0("> FILTER 1 - # of match\t=\t", nrow(all_signif_matching_dt), "\n"))
all_signif_matching_dt  <- all_signif_matching_dt[all_signif_matching_dt$overlapBpRatio >= minOverlapBpRatio,]
cat(paste0("> FILTER 1 - # of match with overlap bp >= ", minOverlapBpRatio, "\t=\t", nrow(all_signif_matching_dt), "\n"))

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


set_of_tads <- foreach(tad = all_signif_tads_with_match) %dopar% {
  # all_signif_matching_dt[all_signif_matching_dt$refID_full == tad,]
  matching_tads <- all_signif_matching_dt$queryID_full[all_signif_matching_dt$refID_full == tad]
  sort(c(tad, matching_tads))
}
names(set_of_tads) <- all_signif_tads_with_match
length(set_of_tads) # 1303
set_of_tads <- unique(set_of_tads)
length(set_of_tads) # 368 !

# check how many times each TAD is present -> to see if a TAD is involved in multiple conserved region
tad_occurence <- table(unlist(set_of_tads))
range(tad_occurence)
# GSE109229_SKBR3_40kb/TCGAbrca_lum_bas-chr6_TAD125
tmpx = "GSE109229_SKBR3_40kb/TCGAbrca_lum_bas/chr6_TAD125"
tmpx = names(head(sort(-tad_occurence),1))
Filter(function(x) any(grepl(tmpx, x)), set_of_tads)

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

cat(paste0("> FILTER 2 - # of conserved regions\t=\t", length(set_of_tads), "\n"))
stopifnot(regionsWithMinGenes %in% names(set_of_tads))
all_set_of_tads  <- set_of_tads
set_of_tads <- set_of_tads[paste0(regionsWithMinGenes)]
cat(paste0("> FILTER 2 - # of regions with # intersect genes >= ", minIntersectGenes, "\t=\t", length(set_of_tads), "\n"))

#######################################################################################################################################
######################################################################################################################### PLOTTING - overview conservation
#######################################################################################################################################                                       
all_datasets <- unique(all_signif_matching_dt$ref_dataset)

# to retrieve the colors
dataset_dt <- data.frame(dataset = all_datasets, hicds = dirname(all_datasets),exprds = basename(all_datasets),stringsAsFactors = FALSE)
dataset_dt$cmpType <- all_cmps[paste0(dataset_dt$exprds)]
dataset_dt$subtype_col <- all_cols[paste0(dataset_dt$cmpType)]
stopifnot(!is.na(dataset_dt$subtype_col))
dataset_dt <- dataset_dt[order(dataset_dt$cmpType, dataset_dt$hicds, dataset_dt$exprds),]
dataset_dt$dataset_label <- gsub("/", "\n", dataset_dt$dataset)
ds_label_levels <- dataset_dt$dataset_label
ds_levels <- dataset_dt$dataset

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
# levels to sort the columns by decreasing # of matches
tmp <- colSums(matching_dt)
region_levels <- names(sort(-tmp))

#### heatmap version
plot_matching_dt <- matching_dt
# rownames(plot_matching_dt) <- gsub("/", "\n", rownames(plot_matching_dt))
initmar <- par()$mar

# modifMar <- c(0,-4.1,-4.1,15)
modifMar <- c(0,0,0,0)

outFile <- file.path(outFolder, paste0("heatmap_conservedRegions_rowReorder_signif", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".", plotType))
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

save(plot_matching_dt, file = file.path(outFolder, "plot_matching_dt.Rdata"), version=2)

par(mar = initmar)

outFile <- file.path(outFolder, paste0("heatmap_conservedRegions_noReorder_signif", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".", plotType))
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

stop("ok\n")

#### ggplot version
# plot_dt <- melt(matching_dt)
# save(plot_dt, file ="plot_dt.Rdata", version=2)
# 
# colnames(plot_dt) <- c("dataset", "region", "signif") 
# 
# matching_dt <- matching_dt[,region_levels]
# image(matching_dt)
# 
# plot_dt$dataset_label <- gsub("/", "\n", plot_dt$dataset)
# plot_dt$dataset_label <- factor(plot_dt$dataset_label, levels=ds_label_levels)
# stopifnot(!is.na(plot_dt$dataset_label))
# 
# plot_dt$region <- factor(plot_dt$region, levels=region_levels)
# stopifnot(!is.na(plot_dt$region))
# plot_dt$signif <- factor(plot_dt$signif, levels = c(1,0))
# 
# g_match <- ggplot(plot_dt, aes(x = region, y = dataset_label, fill = signif, color=signif)) +
#   geom_point() + 
#   scale_fill_manual(values = c("black", "white"))+
#   scale_colour_manual(values = c("black", "white"))

#######################################################################################################################################
######################################################################################################################### lolli plot nTop regions all DS
#######################################################################################################################################                   

nDS <- nrow(matching_dt)

cr= region_levels[1]
foo <- foreach(i_cr = 1:nRegionLolli) %dopar% {
  
  cr <- region_levels[i_cr]
  
  stopifnot(cr %in% names(set_of_tads))
  tads_to_plot <- set_of_tads[[paste0(cr)]]
  
  stopifnot(cr %in% colnames(matching_dt))
  nMatch <- sum(matching_dt[,paste0(cr)])
  
  nPlotted <- length(tads_to_plot)
  
  outHeightGG <- min(c(7 * nPlotted/2, 49))
  outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
  outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
  
  plotList <- list()
  i_tad=1
  for(i_tad in 1:nPlotted) {
  
    hicds <- dirname(dirname(tads_to_plot[i_tad]))
    exprds <- basename(dirname(tads_to_plot[i_tad]))
    tad <- basename(tads_to_plot[i_tad])
      
    mytit <- paste0( hicds, " - ", exprds, " - ", tad)
    
    
    plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                          hicds = hicds,
                                          all_TADs = tad,
                                          orderByLolli = "startPos", mytitle=mytit)
  } # end-for iterating over TADs to plot
  save(plotList, file="plotList.Rdata")
  
  outFile <- file.path(outFolder, paste0("allCmps_conservedRegions", i_cr, "_signif", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_lolli.", plotType))
  
  mytit <- paste0("Conserved region ", i_cr, " (all cmps) - ", nMatch, "/", nDS)
  all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
  outHeightGG <- min(c(7 * nPlotted/2, 49))
  outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
  outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
  
  ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
  cat("... written: ", outFile, "\n")
  
}


#######################################################################################################################################
######################################################################################################################### lolli plot nTop regions by cmp types
#######################################################################################################################################                   

cmpType = unique(as.character(dataset_dt$cmpType))[1]
for(cmpType in unique(as.character(dataset_dt$cmpType))) {
  
  ds_to_keep <- dataset_dt$dataset[as.character(dataset_dt$cmpType) == cmpType]
  stopifnot(ds_to_keep %in% rownames(matching_dt))
  
  cmp_matching_dt <- matching_dt[paste0(ds_to_keep),]
  stopifnot(nrow(cmp_matching_dt) == length(ds_to_keep))
  
  cmp_region_levels <- names(sort(-colSums(cmp_matching_dt)))
  
  nDS_cmp <- nrow(cmp_matching_dt)
  
  
  
  foo <- foreach(i_cr = 1:nRegionLolli) %dopar% {
    
    cr <- region_levels[i_cr]
    
    stopifnot(cr %in% names(set_of_tads))
    tads_to_plot <- set_of_tads[[paste0(cr)]]
    
    stopifnot(cr %in% colnames(matching_dt))
    nMatch <- sum(matching_dt[,paste0(cr)])
    
    nMatch_cmp <- sum(cmp_matching_dt[,paste0(cr)])
    
    nPlotted <- length(tads_to_plot)
    
    outHeightGG <- min(c(7 * nPlotted/2, 49))
    outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
    outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
    
    plotList <- list()
    i_tad=1
    for(i_tad in 1:nPlotted) {
      
      hicds <- dirname(dirname(tads_to_plot[i_tad]))
      exprds <- basename(dirname(tads_to_plot[i_tad]))
      tad <- basename(tads_to_plot[i_tad])
      
      mytit <- paste0( hicds, " - ", exprds, " - ", tad)
      
      
      plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                            hicds = hicds,
                                            all_TADs = tad,
                                            orderByLolli = "startPos", mytitle=mytit)
    } # end-for iterating over TADs to plot
    save(plotList, file="plotList.Rdata")
    
    outFile <- file.path(outFolder, paste0(cmpType, "_conservedRegions", i_cr, "_signif", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_lolli.", plotType))
    
    mytit <- paste0("Conserved region ", i_cr, " (", cmpType, ") - ", nMatch_cmp, "/", nDS_cmp, " (", nMatch, "/", nDS, ")")
    all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
    outHeightGG <- min(c(7 * nPlotted/2, 49))
    outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
    outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
    
    ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
    cat("... written: ", outFile, "\n")
    
  }
  
  
}



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

