
### v2 version (04.03.2018)
# => add threshold intersect/union >= 0.8
# => ! to run: better to set a signif. thresh than a # of top
# => sort lolli by start position !!!

# Rscript intersect_topTADs_acrossDS_v2.R 3
# Rscript intersect_topTADs_acrossDS_v2.R 0.05
# Rscript intersect_topTADs_acrossDS_v2.R 1

startTime <- Sys.time()

cat("> START intersect_topTADs_acrossDS.R \n")

# saved:
# save(signifTADs_allDS_data, file = outFile)
# save(all_matchDT, file = outFile)
# save(all_bestMatchDT, file = outFile)
# save(ratio_matchingSignifTAD_DT, file = outFile)
# save(hicds_exprds_asMatch_DT, file = outFile)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "~/media/electron", "")

require(foreach)
require(doMC)
require(IRanges)
require(GenomicRanges)
require(ggplot2)

source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg"), "set_dataset_colors.R"))
head(score_DT)
source( file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_11_18"), "analysis_utils.R"))
source( file.path("colors_utils.R"))
dataset_proc_colors <- setNames(score_DT$proc_col, score_DT$dataset)
length(dataset_proc_colors)
SSHFS <- FALSE
setDir <- ifelse(SSHFS, "~/media/electron", "")

source("utils_fct.R")
source("plot_lolliTAD_funct.R")
SSHFS <- FALSE
setDir <- ifelse(SSHFS, "~/media/electron", "")

source("utils_plot_fcts.R")

# to be considered as matching TAD: intersect_genes/union_genes should be >= ratioGenMatchingThreshold
ratioGeneMatchingThreshold <- 0.8

registerDoMC(ifelse(SSHFS, 2, 40))

build_signifTADs_allDS_data <- TRUE

topThresh <- 3

# plot only the 0.95 most numerous matching
lolliPlotThreshQuantile <- 0.95
cat(paste0("!!! hard-coded: lolliPlotThreshQuantile = ", lolliPlotThreshQuantile, "\n"))

lolliPlotNtop <- 3
cat(paste0("!!! hard-coded: lolliPlotNtop = ", lolliPlotNtop, "\n"))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0)
topThresh <- as.numeric(args[1])

# if(topThresh == 1) {
#   warning("topThresh == 1 is ambiguous; will be considered as a nTop not pval thresh !\n")
# }
# changed to investigate conserved TADs
if(topThresh == 1) {
  warning("topThresh == 1 is ambiguous; will be considered as pval thresh !\n")
}


outFolder <- file.path("INTERSECT_topTADs_ACROSSDS_v2", paste0("top", topThresh))
dir.create(outFolder, recursive = TRUE)

plotCex <- 1.2
plotType <- "svg"
myHeightDens <- ifelse(plotType=="png", 400, 7)
myWidthDens <- ifelse(plotType=="png", 600, 10)
myHeight <- myWidth <- myHeightDens

myHeightGG <- 7
myWidthGG <- 10

logFile=""

stopifnot(!is.na(topThresh))
stopifnot(topThresh > 0)

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=TRUE, stringsAsFactors = FALSE)
entrezDT$entrezID <- as.character(entrezDT$entrezID)
head(entrezDT)

#PIPELINE/OUTPUT_FOLDER/GSE105318_DLD1_40kb/TCGAcoad_msi_mss//emp_pval_combined.Rdata
script0_name <- "0_prepGeneData"
script11_name <- "11_runEmpPvalCombined"

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

all_hicexpr_ds <- unname(unlist(sapply(list.files(pipOutFolder, full.names = TRUE), function(x) file.path(basename(x),  list.files(x)))))
stopifnot(dir.exists(file.path(pipOutFolder, all_hicexpr_ds)))


ds=all_hicexpr_ds[1]

if(build_signifTADs_allDS_data){
  
cat("... start building signifTADs_allDS_data \n")
signifTADs_allDS_data <- foreach(ds = all_hicexpr_ds) %dopar% {
  
  cat("... start: ", ds, "\n")
  
  hicds <- dirname(ds)
  exprds <- basename(ds)
  stopifnot(dir.exists(hicds))
  dsPipOutDir <- file.path(pipOutFolder, ds)
  stopifnot(dir.exists(dsPipOutDir))
  gene2tadDT_file <- file.path(hicds, "genes2tad/all_genes_positions.txt")
  TADposDT_file <- file.path(hicds, "genes2tad/all_assigned_regions.txt")
  stopifnot(file.exists(gene2tadDT_file))
  stopifnot(file.exists(TADposDT_file))
  
  # RETRIEVE SIGNIF. TADs 
  stopifnot(dir.exists(file.path(dsPipOutDir, script11_name)))
  tad_pvalFile <- file.path(dsPipOutDir, script11_name, "emp_pval_combined.Rdata")
  stopifnot(file.exists(tad_pvalFile))
  tad_pvals <- eval(parse(text = load(tad_pvalFile)))
  adj_tad_pvals <- sort(p.adjust(tad_pvals, method="BH"))
  
  # RETRIEVE PIPELINE GENES
  # subset the genes used in the pipeline
  stopifnot(dir.exists(file.path(dsPipOutDir, script0_name)))
  # geneList: values are entrez, names can be ENSEMBL, entrez, etc.
  geneList_file <- file.path(dsPipOutDir, script0_name, "pipeline_geneList.Rdata")
  stopifnot(file.exists(geneList_file))
  geneList <- eval(parse(text = load(geneList_file)))
  
  # if(topThresh >= 1) {
  if(topThresh > 1) {
    topThresh <- min(c(topThresh, length(adj_tad_pvals)))
    pvalThresh <- as.numeric(adj_tad_pvals[topThresh])
    stopifnot(!is.na(pvalThresh))
    stopifnot(pvalThresh <= 1)
  } else {
    pvalThresh <- topThresh
  }
  top_pvals <- adj_tad_pvals[adj_tad_pvals <= pvalThresh]
  stopifnot(!is.na(top_pvals))
  
  # if(topThresh >= 1) stopifnot(length(unique(top_pvals)) <= topThresh)
  # !!! CHANGED
  if(topThresh > 1) stopifnot(length(unique(top_pvals)) <= topThresh)
  
  
  top_tads <- names(top_pvals)
  
  # RETRIEVE RANKS OF SIGNIF TADs (-> signifDT)
  tad_ranks <- rank(top_pvals, ties="min")
  stopifnot(top_tads %in% names(tad_ranks))
  tad_ranks <- tad_ranks[top_tads]
  
  # RETRIEVE POSITION OF SIGNIF TADs (-> posDT)
  TADposDT <- read.delim(TADposDT_file, stringsAsFactors = FALSE, header=F, col.names=c("chromo","region", "start", "end"))
  stopifnot(is.numeric(TADposDT$start))
  TADposDT <- TADposDT[TADposDT$region %in% top_tads,]
  TADposDT <- TADposDT[match(top_tads, TADposDT$region),]
  stopifnot(nrow(TADposDT) > 0)
  stopifnot(!is.na(TADposDT))
  stopifnot(TADposDT$region %in% top_tads)
  stopifnot(top_tads %in% TADposDT$region)
  
  stopifnot(TADposDT$region == top_tads)
  head(TADposDT)
  
  # RETRIEVE GENES ID AND GENE SYMBOLS FROM SIGNIF TADs (-> geneDT)
  g2tDT <- read.delim(gene2tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
  stopifnot(geneList %in% g2tDT$entrezID)
  stopifnot(top_tads %in% g2tDT$region)
  g2tDT <- g2tDT[g2tDT$entrezID %in% geneList &
                   g2tDT$region %in% top_tads,]
  stopifnot(g2tDT$entrezID %in% entrezDT$entrezID)
  g2t2sDT <- merge(g2tDT[, c("entrezID", "region")], entrezDT[,c("entrezID", "symbol")], by="entrezID", all.x = TRUE, all.y = FALSE )
  stopifnot(!is.na(g2t2sDT))
  g2t2sDT$entrezID <- as.character(g2t2sDT$entrezID)
  g2t2sDT$region <- as.character(g2t2sDT$region)
  g2t2sDT$symbol <- as.character(g2t2sDT$symbol)
  head(g2t2sDT)
  stopifnot(!duplicated(g2t2sDT$entrezID))
  stopifnot(top_tads %in% g2t2sDT$region)
  head(g2t2sDT)
  #geneDT_tmp <- g2t2sDT[order(match(top_tads, g2t2sDT$region), g2t2sDT$entrezID),] # will not work because of the  duplicated values in g2t2sDT$region !!!
  
  geneDT_tmp <- g2t2sDT[order(unlist(sapply(g2t2sDT$region, function(x) which(top_tads == x) )), as.character(g2t2sDT$entrezID)),]
  
  stopifnot(nrow(geneDT_tmp) > 0)
  stopifnot(!is.na(geneDT_tmp))
  stopifnot(geneDT_tmp$region %in% top_tads)
  head(geneDT_tmp)
  stopifnot(top_tads  %in% geneDT_tmp$region)

  geneDT <- data.frame(
    ID = paste(hicds, exprds, as.character(geneDT_tmp$region), sep="_"),
    GENE = as.character(geneDT_tmp$entrezID), 
    SYMBOL = as.character(geneDT_tmp$symbol),
    stringsAsFactors = FALSE
    )
  
  # RETRIEVE MATCHING TADs
  
  # BUIlD THE TABLES:
  # (matchDT built later)
  id_col <- paste(hicds, exprds, top_tads, sep="_")
  
  idDT <- data.frame(
    ID = id_col,
    HICDS = hicds,
    EXPRDS = exprds,
    TAD = top_tads,
    stringsAsFactors = FALSE
  )  
  rownames(idDT) <- NULL
  head(idDT)
  
  signifDT <- data.frame(
    ID = id_col,
    RANK = tad_ranks,              # not sure necessary: could be retrieve rank(PVAL)
    PVAL = as.numeric(top_pvals),
    stringsAsFactors = FALSE
  )
  rownames(signifDT) <- NULL
  head(signifDT)
  
  stopifnot(TADposDT$region == top_tads)
  posDT <- data.frame(
    ID = id_col,
    CHROMO = TADposDT$chromo,
    START  = TADposDT$start,
    END = TADposDT$end,
    stringsAsFactors = FALSE
  )
  rownames(posDT) <- NULL
  head(posDT)
  
  stopifnot(id_col == unique(geneDT$ID))
  
  matchDT <- data.frame(
    ID = id_col,
    MATCHID = NA_character_,
    MATCHRATIO = NA_real_,
    stringsAsFactors = FALSE
  )
  
  list(
    ids=id_col,
    idDT = idDT,
    geneDT = geneDT,
    signifDT=signifDT,
    posDT=posDT,
    matchDT=matchDT
  )
}
names(signifTADs_allDS_data) <- all_hicexpr_ds
outFile <- file.path(outFolder, "signifTADs_allDS_data.Rdata")
save(signifTADs_allDS_data, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


} else { # end-if(build_signifTADs_allDS_data)
  outFile <- file.path(outFolder, "signifTADs_allDS_data.Rdata")
  signifTADs_allDS_data <- eval(parse(text = load(outFile)))
}

all_geneDT <- foreach(ds = names(signifTADs_allDS_data), .combine='rbind') %dopar% {
  signifTADs_allDS_data[[paste0(ds)]][["geneDT"]]
}
printVar("nrow(all_geneDT)")

all_signifDT <- foreach(ds = names(signifTADs_allDS_data), .combine='rbind') %dopar% {
  signifTADs_allDS_data[[paste0(ds)]][["signifDT"]]
}
printVar("nrow(all_signifDT)")

# now I have to build matchDT ...
ds=names(signifTADs_allDS_data)[1]
inputGR_DT <- foreach(ds = names(signifTADs_allDS_data), .combine='rbind') %dopar% {
  posDT <- signifTADs_allDS_data[[paste0(ds)]][["posDT"]]
  idDT <- signifTADs_allDS_data[[paste0(ds)]][["idDT"]]
  outDT <- merge(x=posDT, y=idDT, by="ID", all.x=TRUE, all.y=TRUE )
  head(outDT)
  stopifnot(!is.na(outDT))
  outDT
}

inputGR_DT$ID <- as.character(inputGR_DT$ID)
inputGR_DT$CHROMO <- as.character(inputGR_DT$CHROMO)
inputGR_DT$HICDS <- as.character(inputGR_DT$HICDS)
inputGR_DT$EXPRDS <- as.character(inputGR_DT$EXPRDS)
stopifnot(is.numeric(inputGR_DT$START))
stopifnot(is.numeric(inputGR_DT$END))
head(inputGR_DT)

nAllTADs <- nrow(inputGR_DT)
printVar("nrow(inputGR_DT)")

query_IR <- IRanges(start = inputGR_DT$START, 
                  width = (inputGR_DT$END - inputGR_DT$START + 1), 
                  names=inputGR_DT$ID)

query_GR <- GRanges(ranges = query_IR,
                    seqnames=inputGR_DT$CHROMO)
mcols(query_GR)$ID <- inputGR_DT$ID
stopifnot( names(query_GR) == mcols(query_GR)$ID )
head(mcols(query_GR))

cat("... compute overlap among all TADs")
IDoverlap_hits_all <- findOverlaps(query=query_GR,
                           subject=query_GR)
IDoverlap_hits <- IDoverlap_hits_all[queryHits(IDoverlap_hits_all) != subjectHits(IDoverlap_hits_all)] 

# remove matching with itself
printVar("length(IDoverlap_hits_all)")
printVar("length(IDoverlap_hits)")
stopifnot( (length(IDoverlap_hits_all)-length(IDoverlap_hits) ) == length(query_IR) )

queryID <- names(query_GR[queryHits(IDoverlap_hits)])
matchingID <- names(query_GR[subjectHits(IDoverlap_hits)])
IDsOverlapDT <- data.frame(
  queryID = queryID,
  matchingID = matchingID,
  overlapBP = width(pintersect(query_IR[queryID], query_IR[matchingID])),
    stringsAsFactors = FALSE)
head(IDsOverlapDT)

# check that the matching was done intrachromosomally only !
queryChromo <- gsub(".+_(.+?)_TAD.+?$","\\1", IDsOverlapDT$queryID)
matchingChromo <- gsub(".+_(.+?)_TAD.+?$","\\1", IDsOverlapDT$matchingID)
stopifnot(queryChromo == matchingChromo)

printVar("nAllTADs")
all_query_IDs <- unique(IDsOverlapDT$queryID)
printVar("length(all_query_IDs)")

query_id <- all_query_IDs[1]

head(all_geneDT)
stopifnot(is.character(all_geneDT$ID))
stopifnot(is.character(all_geneDT$GENE))
stopifnot(is.character(all_geneDT$SYMBOL))

cat("... start iterating over all queryIDs to build all_matchDT_noNA\n")

all_matchDT_noNA <- foreach(query_id = all_query_IDs, .combine='rbind') %dopar% {
  
  cat("...... start ", query_id, "\n")
  
  query_genes <- all_geneDT$GENE[all_geneDT$ID == query_id]
  stopifnot(length(query_genes) > 0)
  
  all_matching_IDs <- IDsOverlapDT$matchingID[IDsOverlapDT$queryID == query_id]
  nMatches <- length(all_matching_IDs)
  stopifnot(nMatches > 0)
  
  all_matching_ratio_unstd <- sapply(all_matching_IDs, function(matching_id) {
    matching_genes <- all_geneDT$GENE[all_geneDT$ID == matching_id]
    stopifnot(length(matching_genes) > 0)
    intersectGenes <- intersect(matching_genes, query_genes)
    unionGenes <- union(matching_genes, query_genes)
    length(intersectGenes)/length(unionGenes)
   })
  stopifnot(names(all_matching_ratio_unstd) == all_matching_IDs)
  all_matching_ratio <- sort(all_matching_ratio_unstd, decreasing = TRUE)
  stopifnot( length(all_matching_ratio) == nMatches)
  
  query_hicds <- inputGR_DT$HICDS[inputGR_DT$ID == query_id]
  query_exprds <- inputGR_DT$EXPRDS[inputGR_DT$ID == query_id]
  query_tad <- inputGR_DT$TAD[inputGR_DT$ID == query_id]
  stopifnot( paste0(query_hicds, "_", query_exprds, "_", query_tad) == query_id)
  
  sortedIDs <- names(all_matching_ratio)
  stopifnot(!duplicated(sortedIDs))
  
  matching_hicds <- as.character(sapply(sortedIDs, function(matching_id) {
    inputGR_DT$HICDS[inputGR_DT$ID == matching_id]
  }))
  stopifnot( length(matching_hicds) == length(sortedIDs) )
  
  matching_exprds <- as.character(sapply(sortedIDs, function(matching_id) {
    inputGR_DT$EXPRDS[inputGR_DT$ID == matching_id]
  }))
  stopifnot( length(matching_exprds) == length(sortedIDs) )

   matching_tads <- as.character(sapply(sortedIDs, function(matching_id) {
     inputGR_DT$TAD[inputGR_DT$ID == matching_id]
   }))
   stopifnot( length(matching_tads) == length(sortedIDs) )
   
  all_matching_ratio <- as.numeric(all_matching_ratio)
  stopifnot(!is.na(all_matching_ratio))
 
  matchDT <- data.frame(
    query_id = query_id,
   query_hicds = rep(query_hicds, nMatches),
   query_exprds = rep(query_exprds, nMatches),
   queryTAD = rep(query_tad, nMatches),
   
   matching_id = sortedIDs,
   matching_hicds = matching_hicds,
   matching_exprds = matching_exprds,
   matchingTAD = matching_tads,
   
   matchingRatio = all_matching_ratio,
   
   stringsAsFactors = FALSE
  )
}

# ADD NA FOR THE ONES WITHOUT ANY MATCH
all_ds_tads <- as.character(inputGR_DT$ID)
stopifnot( length(all_ds_tads) == nAllTADs)

noMatch_tads <- all_ds_tads[!all_ds_tads %in% all_matchDT_noNA$query_id]
printVar("nAllTADs")
printVar("length(noMatch_tads)")

cat("... start iterating over IDs without match queryIDs to build all_na_matchDT\n")

all_na_matchDT <- foreach(na_id = noMatch_tads, .combine='rbind') %dopar% {
  id_hicds <- inputGR_DT$HICDS[inputGR_DT$ID == na_id]
  id_exprds <- inputGR_DT$EXPRDS[inputGR_DT$ID == na_id]
  id_tad <- inputGR_DT$TAD[inputGR_DT$ID == na_id]
  stopifnot( paste0(id_hicds, "_", id_exprds, "_", id_tad) == na_id)
  na_matchDT <- data.frame(
    query_id = na_id,
    query_hicds = id_hicds,
    query_exprds = id_exprds,
    queryTAD = id_tad,
    matching_id = NA_character_,
    matching_hicds = NA_character_,
    matching_exprds = NA_character_,
    matchingTAD = NA_character_,
    matchingRatio = NA_real_,
    stringsAsFactors = FALSE
  )
}


all_matchDT <- rbind(all_matchDT_noNA, all_na_matchDT)
stopifnot( length(unique(all_matchDT$query_id)) == nAllTADs)
outFile <- file.path(outFolder, "all_matchDT.Rdata")
save(all_matchDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

cat(paste0("... start building all_bestMatchDT \n"))

# for a given hicds and exprds -> select the best match TAD
all_bestMatchDT <- do.call(rbind,
        lapply(split(all_matchDT,list(all_matchDT$query_id,all_matchDT$matching_hicds,all_matchDT$matching_exprds ),drop=T), 
          function(subDT) subDT[which.max(subDT$matchingRatio),]))
rownames(all_bestMatchDT) <- NULL

all_bestMatchDT[
  all_bestMatchDT$query_id == "ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_lowInf_highInf_chr7_TAD58" &
    all_bestMatchDT$matching_hicds == "K562_40kb" &
    all_bestMatchDT$matching_exprds == "TCGAlaml_wt_mutFLT3" ,
]
all_matchDT[
  all_matchDT$query_id == "ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_lowInf_highInf_chr7_TAD58" &
    all_matchDT$matching_hicds == "K562_40kb" &
    all_matchDT$matching_exprds == "TCGAlaml_wt_mutFLT3" ,
  ]

stopifnot(!duplicated(all_bestMatchDT[, c("query_id", "matching_hicds", "matching_exprds")]))

printVar("length(unique(all_bestMatchDT$query_id))")
printVar("length(unique(na.omit(all_matchDT)$query_id))")
printVar("nAllTADs")

stopifnot( length(unique(all_bestMatchDT$query_id)) == length(unique(na.omit(all_matchDT)$query_id)))
# stopifnot( length(unique(all_bestMatchDT$query_id)) == nAllTADs) # not TRUE because of the NA

outFile <- file.path(outFolder, "all_bestMatchDT.Rdata")
save(all_bestMatchDT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


#*************************************************************************************************
#************************************************************************************************* TAD-level analysis
#*************************************************************************************************

################################ 04.03.2019 ADDED HERE THE THRESHOLD FOR THE V2 VERSION ############################################################# 

all_bestMatchDT <- all_bestMatchDT[all_bestMatchDT$matchingRatio >= ratioGeneMatchingThreshold , ]

#############################################################

cat(paste0("... start computing statistics queryTAD level \n"))

# compute the number of datasets in which signif (with or without counting same exrpds)
all_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT)
head(all_nMatchDT)
colnames(all_nMatchDT) <- c("query_id", "all_nMatch")

diffExprds_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT[all_bestMatchDT$query_exprds != all_bestMatchDT$matching_exprds,])
head(diffExprds_nMatchDT)
colnames(diffExprds_nMatchDT) <- c("query_id", "diffExprds_nMatch")

queryID_matchDT <- merge(all_nMatchDT, diffExprds_nMatchDT, by="query_id", all.x = TRUE, all.y=TRUE)
stopifnot(!is.na(queryID_matchDT$all_nMatch))
sum(is.na(queryID_matchDT$diffExprds_nMatch))

stopifnot(!duplicated(queryID_matchDT$query_id))

####################
#***************************************** multiDens and cumsum
####################


### CUMSUM plots, also with lines for subtypes
subTypeDT <- data.frame(query_exprds=names(cancer_subAnnot), subtype=cancer_subAnnot, stringsAsFactors = FALSE, row.names = NULL)
colorDT <- data.frame(query_exprds=names(dataset_proc_colors), color=dataset_proc_colors, stringsAsFactors = FALSE, row.names = NULL)

tmpDT <- queryID_matchDT
tmpDT <- unique(merge(queryID_matchDT, all_matchDT[, c("query_id", "query_exprds")], by="query_id", all.x=T, all.y=F))
stopifnot(!is.na(tmpDT$query_exprds))
stopifnot(!duplicated(tmpDT$query_id))
tmpDT <- merge(tmpDT, subTypeDT[, c("query_exprds", "subtype")], by ="query_exprds", all.x=T, all.y=F)
tmpDT <- merge(tmpDT, colorDT[, c("query_exprds", "color")], by ="query_exprds", all.x=T, all.y=F)
head(tmpDT)
stopifnot(!is.na(tmpDT$query_exprds))
stopifnot(!duplicated(tmpDT$query_id))

outFile <- file.path(outFolder, paste0("cumNbrTADs_by_nbrDSmatch_allDS_and_onlyDiffExprDS_multidens", ".", plotType))
printVar("outFile")
printVar("myHeightDens")
printVar("myWidthDens")
do.call(plotType, list(outFile, height=myHeightDens, width=myWidthDens))
plot_multiDens(list(
  nMatch_all = queryID_matchDT$all_nMatch,
  nMatch_diffExprds = queryID_matchDT$diffExprds_nMatch),
  plotTit="# datasets with matching signif. TAD", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab=""
)
foo <- dev.off()
cat("... written: ", outFile, "\n")


outFile <- file.path(outFolder, paste0("cumNbrTADs_by_nbrDSmatch_allDS", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_cumMatch(queryID_matchDT, "all_nMatch")
foo <- dev.off()
cat("... written: ", outFile, "\n")

# WITH LINES FOR SUB
outFile <- file.path(outFolder, paste0("cumNbrTADs_by_nbrDSmatch_allDS_withSub", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_cumMatch(queryID_matchDT, "all_nMatch")
### ADD LINES BY SUBTYPE:
tomatch="all_nMatch"
curr_match <- na.omit(tmpDT[, tomatch])
xvect <- seq_len(max(curr_match))
lineVect <- lapply(split(tmpDT, tmpDT$subtype), function(x) {
  curr_match <- na.omit(x[, tomatch])
  yvect <- sapply(xvect, function(i){
    sum(curr_match >= i)
  })
  yvect
})
lineCols <- cancer_subColors[as.character(levels(as.factor(tmpDT$subtype)))]
foo <- sapply(seq_along(lineVect), function(i)  lines(x = rep(list(xvect), length(lineVect))[[i]], y = lineVect[[i]], col = lineCols[i]) )
legend("topright", bty="n", 
       legend = c("all", names(lineCols)),
       lty=1,
       col = c("black", lineCols)
)
foo <- dev.off()
cat("... written: ", outFile, "\n")


outFile <- file.path(outFolder, paste0("cumNbrTADs_by_nbrDSmatch_onlyDiffExprDS", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_cumMatch(queryID_matchDT, "diffExprds_nMatch")
foo <- dev.off()
cat("... written: ", outFile, "\n")

# WITH LINES FOR SUB
outFile <- file.path(outFolder, paste0("cumNbrTADs_by_nbrDSmatch_onlyDiffExprDS_withSub", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_cumMatch(queryID_matchDT, "diffExprds_nMatch")
tomatch="diffExprds_nMatch"
curr_match <- na.omit(tmpDT[, tomatch])
xvect <- seq_len(max(curr_match))
lineVect <- lapply(split(tmpDT, tmpDT$subtype), function(x) {
  curr_match <- na.omit(x[, tomatch])
  yvect <- sapply(xvect, function(i){
    sum(curr_match >= i)
  })
  yvect
})
lineCols <- cancer_subColors[as.character(levels(as.factor(tmpDT$subtype)))]
foo <- sapply(seq_along(lineVect), function(i)  lines(x = rep(list(xvect), length(lineVect))[[i]], y = lineVect[[i]], col = lineCols[i]) )
legend("topright", bty="n", 
       legend = c("all", names(lineCols)),
       lty=1,
       col = c("black", lineCols)
)
foo <- dev.off()
cat("... written: ", outFile, "\n")

outFile <- file.path(outFolder, paste0("matchingRatio_all_and_best_mutlidens", ".", plotType))
do.call(plotType, list(outFile, height=myHeightDens, width=myWidthDens))
plot_multiDens(list(
  all_matchingRatio = all_matchDT$matchingRatio,
  best_matchingRatio = all_bestMatchDT$matchingRatio),
  plotTit="matchingRatio", legTxt=NULL, legPos="topleft", my_ylab="density", my_xlab=""
)
foo <- dev.off()
cat("... written: ", outFile, "\n")


# ### Problem of dup ??? si un groupe de TADs sont tous best match entre eux -> compté à double dans la courbe ????????

# stop("--ok")

####################
#***************************************** boxplots
####################
outFile <- file.path(outFolder, "inputGR_DT.Rdata")
save(inputGR_DT, file = outFile)
cat("... written: ", outFile, "\n")

outFile <- file.path(outFolder, "queryID_matchDT.Rdata")
save(queryID_matchDT, file = outFile)
cat("... written: ", outFile, "\n")

#### boxplot for the dataset
tmpDT <- inputGR_DT
colnames(tmpDT)[colnames(tmpDT) == "ID"] <- "query_id"
boxplotDT <- merge(queryID_matchDT, tmpDT[, c("query_id", "HICDS", "EXPRDS")], by="query_id", all.x=T, all.y=F)
head(boxplotDT)
boxplotDT$dataset <- paste0(boxplotDT$HICDS, "\n", boxplotDT$EXPRDS)

plotylab <- paste0("# matching hicds_exprds_ datasets")
plotTit <- paste0("# of other datasets (all) with signifTADs matching by hicds")
plotSub <- paste0("(# hicds = ", length(unique(boxplotDT$HICDS)), ")")

outFile <- file.path(outFolder, paste0("hicds_all_nMatch_boxplot.", plotType))
p_all <- ggplot_boxplot_hicdsexprds(barDT=boxplotDT, xvar="HICDS", yvar="all_nMatch",
                           colvar=NULL,
                           myylab=plotylab, myTit=plotTit,
                           mySub=plotSub) 
ggsave(plot = p_all, filename = outFile, height=myHeightGG, width = myWidthGG)
cat("... written: ", outFile, "\n")

plotylab <- paste0("# matching hicds_exprds_ datasets")
plotTit <- paste0("# of other datasets (diffExprds) with signifTADs matching by hicds")
plotSub <- paste0("(# hicds = ", length(unique(boxplotDT$HICDS)), ")")

outFile <- file.path(outFolder, paste0("hicds_diffExprds_nMatch_boxplot.", plotType))
p_diff <- ggplot_boxplot_hicdsexprds(barDT=boxplotDT, xvar="HICDS", yvar="diffExprds_nMatch",
                           colvar=NULL,
                           myylab=plotylab, myTit=plotTit,
                           mySub=plotSub) 
ggsave(plot = p_diff, filename = outFile, height=myHeightGG, width = myWidthGG)
cat("... written: ", outFile, "\n")

# scatterplot
plotSub <- paste0("(# hicds = ", length(unique(na.omit(boxplotDT)$HICDS)), ")")
outFile <- file.path(outFolder, paste0("hicds_all_vs_diffExprds_nMatch_scatterplot.", plotType))
plotTit <- paste0("# of other datasets with signifTADs matching by hicds")
do.call(plotType, list(outFile, height=myHeightDens, width=myHeightDens))
plot(all_nMatch ~ diffExprds_nMatch, 
     pch=16, cex=0.7,
     cex.lab=plotCex,
     cex.axis=plotCex,
     data = na.omit(boxplotDT),
     main=plotTit
     )
mtext(side=3, text = plotSub)
foo <- dev.off()
cat("... written: ", outFile, "\n")


####################
#***************************************** lolliplot - quantile
####################

# load(file.path(outFolder, "all_bestMatchDT.Rdata"))
# 
# all_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT)
# head(all_nMatchDT)
# colnames(all_nMatchDT) <- c("query_id", "all_nMatch")
# 
# diffExprds_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT[all_bestMatchDT$query_exprds != all_bestMatchDT$matching_exprds,])
# head(diffExprds_nMatchDT)
# colnames(diffExprds_nMatchDT) <- c("query_id", "diffExprds_nMatch")
# 
# queryID_matchDT <- merge(all_nMatchDT, diffExprds_nMatchDT, by="query_id", all.x = TRUE, all.y=TRUE)
# stopifnot(!is.na(queryID_matchDT$all_nMatch))
# sum(is.na(queryID_matchDT$diffExprds_nMatch))
# stopifnot(!duplicated(queryID_matchDT$query_id))

srtd_DT <- queryID_matchDT[order(queryID_matchDT$all_nMatch, decreasing = T),]

nMatchToPlot <- as.numeric(quantile(srtd_DT$all_nMatch, probs = lolliPlotThreshQuantile))

top_to_plot <- max(which(srtd_DT$all_nMatch >= nMatchToPlot))

plotted_sets <- list()

cat(paste0("... start lolliplot plotting\n"))
i=1

for(i in seq_len(top_to_plot)) {

  plot_list <- list()
  
  tad_id <- srtd_DT$query_id[i]
  matching_ids <- all_bestMatchDT$matching_id[all_bestMatchDT$query_id == tad_id]
  stopifnot(!duplicated(matching_ids))
  
  to_plot_set <- c(tad_id, matching_ids)
  stopifnot(!duplicated(to_plot_set))
  
  if(any(unlist(lapply(plotted_sets, function(x) all(to_plot_set %in% x))))) {
    next
  } else {
    plotted_sets <- c(plotted_sets, list(to_plot_set))
    stopifnot(list(to_plot_set) %in% plotted_sets)
  }
  
  stopifnot(any(unlist(lapply(plotted_sets, function(x) length(setdiff(to_plot_set, x)) == 0))))
  
  plot_list[[1]] <- plot_lolliTAD_ds(exprds = unique(all_bestMatchDT$query_exprds[all_bestMatchDT$query_id == tad_id]), 
                                  hicds = unique(all_bestMatchDT$query_hicds[all_bestMatchDT$query_id == tad_id]), 
                                  all_TADs = unique(all_bestMatchDT$queryTAD[all_bestMatchDT$query_id == tad_id]))
  j <- 2
  for(matchTAD in matching_ids) {
    plot_list[[j]] <- plot_lolliTAD_ds(exprds = unique(all_bestMatchDT$query_exprds[all_bestMatchDT$query_id == matchTAD]), 
                                       hicds = unique(all_bestMatchDT$query_hicds[all_bestMatchDT$query_id == matchTAD]), 
                                       all_TADs = unique(all_bestMatchDT$queryTAD[all_bestMatchDT$query_id == matchTAD]))
    j <- j+1
  }
  
  all_plots <- do.call(grid.arrange, c(plot_list,  list(ncol=2, top=textGrob(paste(tad_id),
                                                                             gp=gpar(fontsize=20,font=2)))))
  # cat("myHeight =", outHeight, "\n")
  outWidth <- 20
  outHeight <- min(c(7 * length(plot_list)/2, 49))
  
  outFile <- file.path(outFolder, paste0(tad_id, "_quantile", lolliPlotThreshQuantile, "_", paste0("nMatch", length(matching_ids)), ".", plotType))
  ggsave(filename = outFile, all_plots, width=outWidth, height = outHeight)
  cat("... written: ", outFile, "\n")

}


####################
#***************************************** lolliplot - nTop
####################

# load(file.path(outFolder, "all_bestMatchDT.Rdata"))
#
# all_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT)
# head(all_nMatchDT)
# colnames(all_nMatchDT) <- c("query_id", "all_nMatch")
#
# diffExprds_nMatchDT <- aggregate(matching_id ~ query_id, FUN=length, data = all_bestMatchDT[all_bestMatchDT$query_exprds != all_bestMatchDT$matching_exprds,])
# head(diffExprds_nMatchDT)
# colnames(diffExprds_nMatchDT) <- c("query_id", "diffExprds_nMatch")
#
# queryID_matchDT <- merge(all_nMatchDT, diffExprds_nMatchDT, by="query_id", all.x = TRUE, all.y=TRUE)
# stopifnot(!is.na(queryID_matchDT$all_nMatch))
# sum(is.na(queryID_matchDT$diffExprds_nMatch))
# stopifnot(!duplicated(queryID_matchDT$query_id))

srtd_DT <- queryID_matchDT[order(queryID_matchDT$all_nMatch, decreasing = T),]



# lolliPlotNtop

plotted_sets <- list()

cat(paste0("... start lolliplot plotting\n"))
i=1

nPlotted <- i<-  0
while(nPlotted < lolliPlotNtop) {
  i <- i+1

  plot_list <- list()

  tad_id <- srtd_DT$query_id[i]
  matching_ids <- all_bestMatchDT$matching_id[all_bestMatchDT$query_id == tad_id]
  stopifnot(!duplicated(matching_ids))

  to_plot_set <- c(tad_id, matching_ids)
  stopifnot(!duplicated(to_plot_set))

  if(any(unlist(lapply(plotted_sets, function(x) all(to_plot_set %in% x))))) {
    next
  } else {
    plotted_sets <- c(plotted_sets, list(to_plot_set))
    stopifnot(list(to_plot_set) %in% plotted_sets)
    nPlotted <- nPlotted + 1
  }

  stopifnot(any(unlist(lapply(plotted_sets, function(x) length(setdiff(to_plot_set, x)) == 0))))

  plot_list[[1]] <- plot_lolliTAD_ds(exprds = unique(all_bestMatchDT$query_exprds[all_bestMatchDT$query_id == tad_id]),
                                     hicds = unique(all_bestMatchDT$query_hicds[all_bestMatchDT$query_id == tad_id]),
                                     all_TADs = unique(all_bestMatchDT$queryTAD[all_bestMatchDT$query_id == tad_id]))
  j <- 2
  for(matchTAD in matching_ids) {
    plot_list[[j]] <- plot_lolliTAD_ds(exprds = unique(all_bestMatchDT$query_exprds[all_bestMatchDT$query_id == matchTAD]),
                                       hicds = unique(all_bestMatchDT$query_hicds[all_bestMatchDT$query_id == matchTAD]),
                                       all_TADs = unique(all_bestMatchDT$queryTAD[all_bestMatchDT$query_id == matchTAD]))
    j <- j+1
  }

  all_plots <- do.call(grid.arrange, c(plot_list,  list(ncol=2, top=textGrob(paste(tad_id),
                                                                             gp=gpar(fontsize=20,font=2)))))
  # cat("myHeight =", outHeight, "\n")
  outWidth <- 20
  outHeight <- min(c(7 * length(plot_list)/2, 49))

  outFile <- file.path(outFolder, paste0(tad_id, "_nTop", lolliPlotNtop, "_", paste0("nMatch", length(matching_ids)), ".", plotType))
  ggsave(filename = outFile, all_plots, width=outWidth, height = outHeight)
  cat("... written: ", outFile, "\n")

}

#*************************************************************************************************
#************************************************************************************************* hicds_exprds-level analysis
#*************************************************************************************************

cat(paste0("... start computing statistics hicds_exprds level \n"))

hicds_exprds_DT <- unique(all_bestMatchDT[, c("query_hicds", "query_exprds")])

#***********
### Pour chaque hicds/exprds -> % de ces signif TADs qui matchent un nombre "x" de signif TADs d'autres datasets  ------ BUILD TABLE FOR BARPLOT ratio_matchingSignifTAD_DT
#***********
i=1
ratio_matchingSignifTAD_DT <- foreach(i = seq_len(nrow(hicds_exprds_DT)), .combine='rbind') %dopar% {
  i_hicds <- hicds_exprds_DT$query_hicds[i]
  i_exprds <- hicds_exprds_DT$query_exprds[i] 
  
  # retrieve the tad query_id for this hicds_exprds dataset
  signifTADs <- unique(all_bestMatchDT$query_id[all_bestMatchDT$query_hicds == i_hicds &
                                           all_bestMatchDT$query_exprds == i_exprds
                                           ])
  nSignif <- length(signifTADs)  
  # sapply(signifTADs, function(x) {
  #   sum(all_bestMatchDT$matching_id == x)
  # })
  
  # retrieve the tad query_id for this hicds_exprds dataset
  other_hicds_exprds_DT <- all_bestMatchDT[all_bestMatchDT$query_hicds != i_hicds &
                                             all_bestMatchDT$query_exprds != i_exprds,
                                           ]
  
  
  signifTADs_as_matching <- sapply(signifTADs, function(x) {
    as.numeric(any(all_bestMatchDT$matching_id == x))
  })
  
  signifTADs_as_matching_otherds <- sapply(signifTADs, function(x) {
    as.numeric(any(other_hicds_exprds_DT$matching_id == x))
  })
  
  nSignif_as_matching <- sum(signifTADs_as_matching)
  nSignif_as_matching_otherds <- sum(signifTADs_as_matching_otherds)
  
  data.frame(
    hicds = i_hicds,
    exprds = i_exprds,
    nSignifTADs = nSignif,
    ratioSignifTADs_as_bestMatch = nSignif_as_matching/nSignif,
    ratioSignifTADs_as_bestMatch_diffds = nSignif_as_matching_otherds/nSignif,
    stringsAsFactors = FALSE
  )
}

outFile <- file.path(outFolder, "ratio_matchingSignifTAD_DT.Rdata")
save(ratio_matchingSignifTAD_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

stopifnot(ratio_matchingSignifTAD_DT$ratioSignifTADs_as_bestMatch <= 1)
stopifnot(ratio_matchingSignifTAD_DT$ratioSignifTADs_as_bestMatch >= 0)
range(ratio_matchingSignifTAD_DT$ratioSignifTADs_as_bestMatch)

# > head(ratio_matchingSignifTAD_DT)
# hicds                   exprds nSignifTADs ratioSignifTADs_as_bestMatch
# 1     ENCSR079VIJ_G401_40kb       TCGAkich_norm_kich           2                       1.0000
# 2 ENCSR312KHQ_SK-MEL-5_40kb  TCGAskcm_lowInf_highInf          16                       0.9375
ratio_matchingSignifTAD_DT <- ratio_matchingSignifTAD_DT[order(ratio_matchingSignifTAD_DT$ratioSignifTADs_as_bestMatch, decreasing=TRUE),]
ratio_matchingSignifTAD_DT$dataset <- paste0(ratio_matchingSignifTAD_DT$hicds, "\n", ratio_matchingSignifTAD_DT$exprds)
ratio_matchingSignifTAD_DT$dataset <- factor(ratio_matchingSignifTAD_DT$dataset, levels = as.character(ratio_matchingSignifTAD_DT$dataset))

mynames <- ratio_matchingSignifTAD_DT$exprds
curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[mynames])])
stopifnot(!is.na(curr_colors))
ratio_matchingSignifTAD_DT$dscols <- curr_colors

plotTit <- paste0("ratio DS has bestMatch with other DS")
plotSub <- paste0("(# ds = ", length(unique(ratio_matchingSignifTAD_DT$dataset)), ")")

####################
#***************************************** BARPLOT - ratio_matchingSignifTAD_DT
####################

p_bestMatch <- ggplot_barplot_hicdsexprds(barDT=ratio_matchingSignifTAD_DT, 
                                               xvar="dataset", 
                                               yvar="ratioSignifTADs_as_bestMatch",
                                               xcolvar = "dscols",
                                               myxlab="", 
                                               myylab="ratioDS_is_bestMatch",
                                               myTit=plotTit,
                                               mySub=plotSub
)
outFile <- file.path(outFolder, paste0("ratioSignifTADs_as_bestMatch", "_", "barplot.", plotType))
ggsave(plot = p_bestMatch, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

####################
#***************************************** BARPLOT - ratio_matchingSignifTAD_DT - other hicds_exprds
####################
ratio_matchingSignifTAD_DT <- ratio_matchingSignifTAD_DT[order(ratio_matchingSignifTAD_DT$ratioSignifTADs_as_bestMatch_diffds, decreasing=TRUE),]
ratio_matchingSignifTAD_DT$dataset <- paste0(ratio_matchingSignifTAD_DT$hicds, "\n", ratio_matchingSignifTAD_DT$exprds)
ratio_matchingSignifTAD_DT$dataset <- factor(ratio_matchingSignifTAD_DT$dataset, levels = as.character(ratio_matchingSignifTAD_DT$dataset))

mynames <- ratio_matchingSignifTAD_DT$exprds
curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[mynames])])
stopifnot(!is.na(curr_colors))
ratio_matchingSignifTAD_DT$dscols <- curr_colors

plotTit <- paste0("ratio DS has bestMatch with other DS (diffDS)")

p_bestMatch <- ggplot_barplot_hicdsexprds(barDT=ratio_matchingSignifTAD_DT, 
                                          xvar="dataset", 
                                          yvar="ratioSignifTADs_as_bestMatch_diffds",
                                          xcolvar = "dscols",
                                          myxlab="", 
                                          myylab="ratioDS_is_bestMatch",
                                          myTit=plotTit,
                                          mySub=plotSub
)
outFile <- file.path(outFolder, paste0("ratioSignifTADs_as_bestMatch_diffds", "_", "barplot.", plotType))
ggsave(plot = p_bestMatch, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

#***********
### Pour chaque hicds/exprds, combien de fois ce dataset est dans le best matching TAD des autres datasets ------ BUILD TABLE FOR BARPLOT hicds_exprds_asMatch_DT
#***********
i=1
hicds_exprds_asMatch_DT <- foreach(i = seq_len(nrow(hicds_exprds_DT)), .combine='rbind') %dopar% {
  
  i_hicds <- hicds_exprds_DT$query_hicds[i]
  i_exprds <- hicds_exprds_DT$query_exprds[i] 
  
  # retrieve the tad query_id for this hicds_exprds dataset
  other_hicds_exprds_DT <- all_bestMatchDT[all_bestMatchDT$query_hicds != i_hicds &
                                                  all_bestMatchDT$query_exprds != i_exprds,
                                                ]
  
  nQuery <- length(unique(other_hicds_exprds_DT$query_id))
  
  nDS_asMatching <- sum(other_hicds_exprds_DT$matching_exprds == i_exprds & other_hicds_exprds_DT$matching_hicds == i_hicds)
  
  stopifnot(nDS_asMatching <= nQuery)
# sum(all_bestMatchDT$query_hicds == hicds_exprds_DT$query_hicds[i] &
#                                                   all_bestMatchDT$query_exprds == hicds_exprds_DT$query_exprds[i] 
#                                                 )
  data.frame(
    hicds = i_hicds,
    exprds = i_exprds,
    nQueryTADs = nQuery,
    ratioDS_as_queryBestMatch = nDS_asMatching/nQuery,
    stringsAsFactors = FALSE
  )
}

outFile <- file.path(outFolder, "hicds_exprds_asMatch_DT.Rdata")
save(hicds_exprds_asMatch_DT, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

stopifnot(hicds_exprds_asMatch_DT$ratioDS_as_queryBestMatch <= 1)
stopifnot(hicds_exprds_asMatch_DT$ratioDS_as_queryBestMatch >= 0)
range(hicds_exprds_asMatch_DT$ratioDS_as_queryBestMatch)
# > head(hicds_exprds_asMatch_DT)
# hicds                   exprds nQueryTADs ratioDS_as_queryBestMatch
# 1     ENCSR079VIJ_G401_40kb       TCGAkich_norm_kich        352                0.02272727

####################
#***************************************** BARPLOT - hicds_exprds_asMatch_DT
####################


hicds_exprds_asMatch_DT <- hicds_exprds_asMatch_DT[order(hicds_exprds_asMatch_DT$ratioDS_as_queryBestMatch, decreasing=TRUE),]
hicds_exprds_asMatch_DT$dataset <- paste0(hicds_exprds_asMatch_DT$hicds, "\n", hicds_exprds_asMatch_DT$exprds)
hicds_exprds_asMatch_DT$dataset <- factor(hicds_exprds_asMatch_DT$dataset, levels = as.character(hicds_exprds_asMatch_DT$dataset))

mynames <- hicds_exprds_asMatch_DT$exprds
curr_colors <- as.character(cancer_subColors[as.character(cancer_subAnnot[mynames])])
stopifnot(!is.na(curr_colors))
hicds_exprds_asMatch_DT$dscols <- curr_colors

plotTit <- paste0("ratio DS has queryBestMatch with other DS")
plotSub <- paste0("(# ds = ", length(unique(hicds_exprds_asMatch_DT$dataset)), ")")

p_queryBestMatch <- ggplot_barplot_hicdsexprds(barDT=hicds_exprds_asMatch_DT, 
                           xvar="dataset", 
                           yvar="ratioDS_as_queryBestMatch",
                           xcolvar = "dscols",
                           myxlab="", 
                           myylab="ratioDS_is_queryBestMatch",
                           myTit=plotTit,
                           mySub=plotSub
                           )
outFile <- file.path(outFolder, paste0("ratioDS_as_queryBestMatch", "_", "barplot.", plotType))
ggsave(plot = p_queryBestMatch, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))
  
  
# scatterplot
scatterDT <- merge(ratio_matchingSignifTAD_DT[, c("hicds", "exprds", "ratioSignifTADs_as_bestMatch","ratioSignifTADs_as_bestMatch_diffds" )], 
                   hicds_exprds_asMatch_DT[c("hicds", "exprds", "ratioDS_as_queryBestMatch")], by=c("hicds", "exprds"))
stopifnot(!is.na(scatterDT))

head(ratio_matchingSignifTAD_DT)
head(hicds_exprds_asMatch_DT)

# scatterplot
plotSub <- paste0("(# hicds_exprds = ", nrow(unique(na.omit(scatterDT))), ")")
outFile <- file.path(outFolder, paste0("ratioSignifTADs_as_bestMatch_vs_as_queryBestMatch_scatterplot.", plotType))
plotTit <- paste0("ratio query as_bestMatch vs. as_queryBestMatch")
do.call(plotType, list(outFile, height=myHeightDens, width=myHeightDens))
plot(ratioSignifTADs_as_bestMatch ~ ratioDS_as_queryBestMatch, 
     pch=16, cex=0.7,
     cex.lab=plotCex,
     cex.axis=plotCex,
     data = scatterDT,
     main=plotTit
)
mtext(side=3, text = plotSub)
foo <- dev.off()
cat("... written: ", outFile, "\n")

# scatterplot
plotSub <- paste0("(# hicds_exprds = ", nrow(unique(na.omit(scatterDT))), ")")
outFile <- file.path(outFolder, paste0("ratioSignifTADs_as_bestMatch_diffds_vs_as_queryBestMatch_scatterplot.", plotType))
plotTit <- paste0("ratio query as_bestMatch vs. as_queryBestMatch (diffDS)")
do.call(plotType, list(outFile, height=myHeightDens, width=myHeightDens))
plot(ratioSignifTADs_as_bestMatch_diffds ~ ratioDS_as_queryBestMatch, 
     pch=16, cex=0.7,
     cex.lab=plotCex,
     cex.axis=plotCex,
     data = scatterDT,
     main=plotTit
)
mtext(side=3, text = plotSub)
foo <- dev.off()
cat("... written: ", outFile, "\n")
                                                                                   
# ######################################################################################
# ######################################################################################
# ######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

#*************************************************************************************************
#************************************************************************************************* trash/debug
#*************************************************************************************************

# check why in lolliplot but no genes in common

#e.g. ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_lowInf_highInf_chr1_TAD300 with ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258

outFile <- file.path(outFolder, "signifTADs_allDS_data.Rdata")
load(outFile)

load(file.path(outFolder, "all_bestMatchDT.Rdata"))
load(outFile)

load(file.path(outFolder, "all_matchDT.Rdata"))
load(outFile)

queryID="ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258"
all_bestMatchDT[all_bestMatchDT$query_id == queryID,]

all_matchDT[all_matchDT$query_id == queryID,]

str(signifTADs_allDS_data[["ENCSR312KHQ_SK-MEL-5_40kb/TCGAskcm_lowInf_highInf"]])
names(signifTADs_allDS_data[["ENCSR312KHQ_SK-MEL-5_40kb/TCGAskcm_lowInf_highInf"]])
dt1 <- signifTADs_allDS_data[["ENCSR312KHQ_SK-MEL-5_40kb/TCGAskcm_lowInf_highInf"]][["posDT"]]
query1 <- "ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_lowInf_highInf_chr1_TAD300"
dt1[dt1$ID == query1,]
# ID CHROMO     START       END
# 2 ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_lowInf_highInf_chr1_TAD300   chr1 152840001 153280000

dt2 <- signifTADs_allDS_data[["ENCSR346DCU_LNCaP_40kb/TCGAprad_norm_prad"]][["posDT"]]
query2 <- "ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258"
dt2[dt2$ID == query2,]
# ID CHROMO     START       END
# 2 ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258   chr1 153200001 153640000


# outFile <- "INTERSECT_topTADs_ACROSSDS/top3/all_matchDT.Rdata"
# load(outFile)
# all_matchDT[1:5,1:5]
# 
# outFile <- "INTERSECT_topTADs_ACROSSDS/top3/all_bestMatchDT.Rdata"
# load(outFile)
# all_bestMatchDT[1:5,1:5]
# 
# 
# curr_tad = "ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_lowInf_highInf_chr6_TAD71"
# 
# moreThanOneMatch <- all_matchDT[, c("query_id", "matching_hicds", "matching_exprds")]
# tad_withDup_idx <- which(duplicated(moreThanOneMatch))
# 
# tad_withDup <- all_matchDT[tad_withDup_idx[1], c("query_id")]
# "ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_lowInf_highInf_chr7_TAD58"
# 
# nrow(all_matchDT[all_matchDT$query_id == tad_withDup,])
# nrow(all_bestMatchDT[all_bestMatchDT$query_id == tad_withDup,])
# 
# tad_withDup <- all_matchDT[tad_withDup_idx[2], c("query_id")]
# tad_withDup
# nrow(all_matchDT[all_matchDT$query_id == tad_withDup,])
# nrow(all_bestMatchDT[all_bestMatchDT$query_id == tad_withDup,])
# 
# tad_withDup <- all_matchDT[tad_withDup_idx[3], c("query_id")]
# tad_withDup
# nrow(all_matchDT[all_matchDT$query_id == tad_withDup,])
# nrow(all_bestMatchDT[all_bestMatchDT$query_id == tad_withDup,])
# 
# tad_withDup = "ENCSR444WCZ_A549_40kb_TCGAluad_mutKRAS_mutEGFR_chr6_TAD60" 
# tad_withDup = "GSE105194_spinal_cordGSE105194_cerebellum_40kb_TCGAlgg_IDHwt_IDHmutnc_chr17_TAD80"
# 
# all_tads_withDup <- unique(all_matchDT$query_id[tad_withDup_idx]) # almost all only 1 diff
# 
# 
# for(tad_withDup in all_tads_withDup) {
# cat(nrow(all_matchDT[all_matchDT$query_id == tad_withDup,]), "\n")
#   cat(nrow(all_bestMatchDT[all_bestMatchDT$query_id == tad_withDup,]), "\n")
# }


# xvect <- seq_len(max(queryID_matchDT$all_nMatch, na.rm=TRUE))  
# yvect <- sapply(xvect, function(x){
#   sum(queryID_matchDT$all_nMatch >= x)
# })
# 
# plot(x = xvect,
#      y = yvect,
#      xlab = paste0("# datasets in which matching signif. TAD"), 
#      ylab = paste0("# query TAD"),
#      type="l")

# for DEBUG
# outFile <- "INTERSECT_topTADs_ACROSSDS/top3/all_matchDT.Rdata"
# load(outFile)
# all_matchDT[1:5,1:5]
# 
# outFile <- "INTERSECT_topTADs_ACROSSDS/top3/all_bestMatchDT.Rdata"
# load(outFile)
# all_bestMatchDT[1:5,1:5]

# outFile <- "INTERSECT_topTADs_ACROSSDS/all_matchDT.Rdata""
# load(outFile)

# by(all_matchDT, all_matchDT$query_id, function(subDT) {
#   subDT
# })

# all_matchDT <- all_matchDT[c(1:3,200:205),]

# by(all_matchDT, list(all_matchDT$query_id, all_matchDT$matching_hicds, all_matchDT$matching_exprds), function(subDT) {
#  return(1)
# })
# 
# by(all_matchDT, list(all_matchDT$query_id,all_matchDT$matching_hicds ), function(subDT) {
#   subDT[which.max(subDT$matchingRatio),]
# })

# !!! pintersect removes the match with itself !
# pintersect(query_IR["ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258",
#                     "ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258"], 
#            query_IR["ENCSR346DCU_LNCaP_40kb_TCGAprad_norm_prad_chr1_TAD258",
#                     "ENCSR079VIJ_G401_40kb_TCGAkich_norm_kich_chr1_TAD274"])

# load("INTERSECT_topTADs_ACROSSDS/signifTADs_allDS_data.Rdata")
# head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["geneDT"]])
# head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["signifDT"]])
# head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["posDT"]])
# head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["matchDT"]])
# head(signifTADs_allDS_data[["NCI-H460_40kb/TCGAlusc_norm_lusc"]][["idDT"]])






