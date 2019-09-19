startTime <- Sys.time()
cat(paste0("> Rscript enrich_same_diff_TAD.R\n"))

buildTable <- FALSE

# Rscript enrich_same_diff_TAD.R

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)
require(reshape2)
require(ggpubr)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

library(stringr)
pattern<-"connected_gene=.+?;"

col1 <- get_palette("Dark2", 2)[1]
col2 <- get_palette("Dark2", 2)[2]

hicds="Panc1_rep12_40kb"
exprds="TCGApaad_wt_mutKRAS"
  
outFolder <- file.path("ENRICH_SAME_DIFF_TAD")
dir.create(outFolder, recursive = TRUE)

setDir="/media/electron"
setDir=""
entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
stopifnot(!duplicated(entrez2symb_dt$symbol))
rownames(entrez2symb_dt) <- entrez2symb_dt$symbol

# enhancer data download from https://www.genecards.org/GeneHancer_version_4-4
enhancerFile <- "GeneHancer_version_4-4.csv"
stopifnot(file.exists(enhancerFile))
enhancerDT <- read.delim(enhancerFile, stringsAsFactors = FALSE, sep="\t", header=TRUE)
stopifnot(nrow(enhancerDT) > 0)
i=1

if(buildTable) {
  full_enhancerDT <- foreach(i=1:nrow(enhancerDT), .combine='rbind') %dopar% {
    all_genes <- gsub(".+=(.+);", "\\1", unlist(str_extract_all(enhancerDT$attributes[i],pattern)))
    all_genes <- all_genes[all_genes %in% entrez2symb_dt$symbol]
    stopifnot(all_genes %in% rownames(entrez2symb_dt))
    if(length(all_genes) == 0) return(NULL)
    data.frame(
      geneSymbol = all_genes,
      entrezID = entrez2symb_dt[paste0(all_genes), "entrezID"],
      gene_chromo = entrez2symb_dt[paste0(all_genes), "chromo"],
      gene_start = entrez2symb_dt[paste0(all_genes), "start"],
      gene_end = entrez2symb_dt[paste0(all_genes), "end"],
      enhancer_chromo = enhancerDT$chrom[i],
      enhancer_start = enhancerDT$start[i],
      enhancer_end = enhancerDT$end[i],
      stringsAsFactors = FALSE
    )
  }
  outFile <- file.path(outFolder, "full_enhancerDT.Rdata")
  save(full_enhancerDT, file=outFile, version=2)
  cat(paste0("... written: ", outFile))
} else {
  outFile <- file.path(outFolder, "full_enhancerDT.Rdata")
  full_enhancerDT <- get(load(outFile))
}



mainFolder <- file.path(".")

pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"

all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds



if(buildTable) {
  cat("... start preparing data before matching \n")
  
  # all_hicds=all_hicds[1]
  hicds = all_hicds[1]
  all_enhancer_tad_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "gene_chromo", "gene_start", "gene_end", "gene_region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    stopifnot(!duplicated(g2t_dt$entrezID))
    rownames(g2t_dt) <- g2t_dt$entrezID
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(hicds_file))
    tad_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo","region", "start", "end"))
    # tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
    stopifnot(nrow(tad_dt) > 0 )
    
    
    # all_exprds = all_exprds[[paste0(hicds)]][1]
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat(paste0("... start ", hicds, " - ", exprds, "\n"))
      
      dataset_pipDir <- file.path(pipFolder, hicds, exprds)
      stopifnot(dir.exists(dataset_pipDir))
      
      
      
      geneList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      pipeline_geneList <- get(load(geneList_file))
      
      regionList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
      stopifnot(file.exists(regionList_file))
      regionList <- get(load(regionList_file))
      
      stopifnot(pipeline_geneList %in% g2t_dt$entrezID)
      
      exprds_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% pipeline_geneList,]
      stopifnot(setequal(rownames(exprds_g2t_dt), pipeline_geneList))
      stopifnot(rownames(exprds_g2t_dt) == exprds_g2t_dt$entrezID)
      
      exprds_enhancerDT <- full_enhancerDT[full_enhancerDT$entrezID %in% pipeline_geneList,]
      
      # assign the enhancer to TAD
      exprds_enhancer2tad_dt <- unique(exprds_enhancerDT[,c("enhancer_chromo", "enhancer_start", "enhancer_end")])
      exprds_enhancer2tad_dt$enhancerID <- paste0("enhancer", 1:nrow(exprds_enhancer2tad_dt))
      
      # exprds_enhancer2tad_dt <- exprds_enhancer2tad_dt[1:10000,]
      
      if(nrow(exprds_enhancer2tad_dt) == 0) return(NULL)
      
      i=1
      exprds_enhancer2tad_dt$enhancer_region <- foreach(i = 1:nrow(exprds_enhancer2tad_dt), .combine='c') %dopar% {
        match_region <- tad_dt$region[
          tad_dt$chromo == exprds_enhancer2tad_dt$enhancer_chromo[i] &
            tad_dt$start <= exprds_enhancer2tad_dt$start[i] &
            tad_dt$end >= exprds_enhancer2tad_dt$start[i] 
        ]
        if(length(match_region) == 0) return(NA)
        match_region
      }
      
      
      gene_tad_enhancer_dt <- merge(exprds_enhancerDT[,c("entrezID", "enhancer_chromo", "enhancer_start", "enhancer_end")], exprds_g2t_dt[,c("entrezID", "gene_region")], by=c("entrezID"), all.x=TRUE, all.y=FALSE)
      stopifnot(!is.na(gene_tad_enhancer_dt$gene_region))
      
      gene_tad_enhancer_tad_dt <- merge(gene_tad_enhancer_dt, exprds_enhancer2tad_dt, by=c("enhancer_chromo", "enhancer_start", "enhancer_end"), all.x=TRUE, all.y=FALSE )
      
      gene_tad_enhancer_tad_dt$gene_enhancer_TADs <- ifelse(gene_tad_enhancer_tad_dt$gene_region == gene_tad_enhancer_tad_dt$enhancer_region, "same", "diff")
      gene_tad_enhancer_tad_dt$gene_enhancer_TADs[is.na(gene_tad_enhancer_tad_dt$gene_enhancer_TADs)] <- "diff"
      
      save(gene_tad_enhancer_tad_dt, file=file.path(outFolder, paste0(hicds, "_", exprds, "_", "gene_tad_enhancer_tad_dt.Rdata")), version=2)
      # stop("--ok")
      # 
      stopifnot(!is.na(gene_tad_enhancer_tad_dt$entrezID))
      stopifnot(!is.na(gene_tad_enhancer_tad_dt$enhancerID))
      
      # for each TAD -> count # enhancers same TAD , diffTAD
      
      nSame_dt <- aggregate(gene_enhancer_TADs ~ gene_region, data=gene_tad_enhancer_tad_dt, FUN=function(x) sum(x=="same"))
      colnames(nSame_dt)[colnames(nSame_dt)=="gene_enhancer_TADs"] <- "nEnhancersInSameTAD"

      nDiff_dt <- aggregate(gene_enhancer_TADs ~ gene_region, data=gene_tad_enhancer_tad_dt, FUN=function(x) sum(x=="diff"))
      colnames(nDiff_dt)[colnames(nDiff_dt)=="gene_enhancer_TADs"] <- "nEnhancersInDiffTAD"      
      
      
      all_count_dt <- merge(nSame_dt, nDiff_dt, by="gene_region", all=TRUE)
      stopifnot(!is.na(all_count_dt$nEnhancersInSameTAD))
      stopifnot(!is.na(all_count_dt$nEnhancersInDiffTAD))
      
      colnames(all_count_dt)[colnames(all_count_dt) == "gene_region"] <- "region"
      all_cols <- colnames(all_count_dt)
      all_count_dt$hicds <- hicds
      all_count_dt$exprds <- exprds
      all_count_dt[, c("hicds", "exprds", all_cols)]
      
            
    }# end-foreach iterating over exprds
    
    exprds_dt
    
  }# end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, paste0("all_enhancer_tad_dt.Rdata"))
  save(all_enhancer_tad_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0( "all_enhancer_tad_dt.Rdata"))
  cat("... load data\n")
  all_enhancer_tad_dt <- get(load(outFile))
}

plot_dt <- all_enhancer_tad_dt
plot_dt$region <- NULL
plot_dt_m <- melt(plot_dt, id=c("hicds", "exprds"))



final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

p <- ggdensity(plot_dt_m, 
               title = paste0("# gene enhancers"),
               subtitle=paste0(""),
               x = "value", 
               color = "variable", fill = "variable",
               # add = "mean", rug = TRUE,
               xlab = "# enhancers",
               palette = c(col1, col2))

outFile <- file.path(outFolder, paste0("all_ds_nbr_enhancers_same_diff_TAD", "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))

sub_dt <- plot_dt_m[plot_dt_m$variable == "nEnhancersInDiffTAD",]
sub_dt$value_log10 <- log10(sub_dt$value)

p <- ggdensity(sub_dt, 
               title = paste0("# gene enhancers"),
               subtitle=paste0(""),
               x = "value", 
               color = "variable", fill = "variable",
               # add = "mean", rug = TRUE,
               xlab = "# enhancers (log10)",
               palette = c(col1, col2))

outFile <- file.path(outFolder, paste0("all_ds_nbr_enhancers_only_diff_TAD_log10", "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))

# caution there is 1 link gene -> enhancer

pval_enhancer_DT <- merge(all_enhancer_tad_dt, final_table_DT, by=c("hicds", "exprds", "region"))
pval_enhancer_DT$totEnhancers <- pval_enhancer_DT$nEnhancersInSameTAD + pval_enhancer_DT$nEnhancersInDiffTAD


nDS <- length(unique(paste0(pval_enhancer_DT$hicds, pval_enhancer_DT$exprds)))

all_y <- c("totEnhancers", "nEnhancersInSameTAD", "nEnhancersInDiffTAD")

for(plot_y in all_y) {
  
  
  outFile <- file.path(outFolder, paste0("all_ds_nbr_enhancers_only_diff_TAD_log10", "_density.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    main = paste0("# enhancers/TAD and TAD signif."),
    x = -log10(pval_enhancer_DT[, paste0("adjPvalComb")]),
    y=  (pval_enhancer_DT[, paste0(plot_y)]),
    xlab=paste0("adj. pval comb [-log10]"),
    ylab=paste0("# enhancers (", plot_y, ")")
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile,"\n"))
  
  
}





















# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# signif_column <- "adjPvalComb"
# signifThresh <- 0.01
# signifcol <- paste0(signif_column, "_", signifThresh)
# minOverlapBpRatio <- 0.8
# minIntersectGenes <- 3
# 
# args="norm_vs_tumor"
# args <- commandArgs(trailingOnly = TRUE)
# 
# if(length(args) == 0) {
#   inFolder <- "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2"
#   args <- ""
#   file_prefix <- ""
# } else if(length(args) == 1) {
#   inFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", args[1])
#   args <- args[1]
#   file_prefix <- paste0(args, "_")
# } else {
#   stop("error\n")
# }
# stopifnot(dir.exists(inFolder))
# inFile <- file.path(inFolder, paste0(file_prefix,"conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_df.txt"))
# stopifnot(file.exists(inFile))
# conserved_dt <- read.delim(inFile, header=TRUE, stringsAsFactors = FALSE)
# stopifnot(nrow(conserved_dt) > 0)
# 
# outFolder <- file.path("TAD_MATCHING_ENRICH_ENHANCER_ANNOT", args)
# dir.create(outFolder, recursive = TRUE)
# outfile <- file.path(outFolder, "tad_matching_enrich_enhancer_annot.txt")
# file.remove(outfile)
# 
# mainFolder <- file.path(".")
# stopifnot(dir.exists(mainFolder))
# 
# i=2
# for(i in 1:nrow(conserved_dt)) {
#   
#   # cat("i = ", i, "\n")
#   
#   curr_reg <- conserved_dt$conserved_region[i]
#   
#   curr_tads <- conserved_dt$corresp_tads[i]
#   
#   all_regs <- unlist(strsplit(x=curr_tads, split=","))
#   
#   chromo <- unique(unlist(sapply(all_regs, function(reg) gsub("(.+)_.+", "\\1", basename(reg)))))
#   stopifnot(length(chromo) == 1)
#   
#   reg=all_regs[[1]]
#   all_starts_ends <- sapply(all_regs, function(reg) {
#     hicds <- dirname(dirname(reg))
#     tadfile <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
#     tad <- basename(reg)
#     stopifnot(file.exists(tadfile))
#     tad_dt <- read.delim(tadfile, stringsAsFactors = FALSE, header=F, col.names=c("chromo", "region", "start", "end"))
#     stopifnot(tad %in% tad_dt$region)
#     tad_start <- tad_dt$start[tad_dt$region==tad]
#     stopifnot(length(tad_start) == 1)
#     tad_end <- tad_dt$end[tad_dt$region==tad]
#     stopifnot(length(tad_end) == 1)
#     stopifnot( tad_dt$chromo[tad_dt$region==tad] == chromo)
#     c(tad_start, tad_end)
#   })
#   stopifnot(ncol(all_starts_ends) == length(all_regs))
#   stopifnot(nrow(all_starts_ends) == 2)
#   stopifnot(all_starts_ends[1,]<all_starts_ends[2,])
#   min_start <- min(all_starts_ends[1,])
#   stopifnot(is.numeric(min_start))
#   max_end <- max(all_starts_ends[2,])
#   stopifnot(is.numeric(max_end))
#   
#   mycon <- file(outfile,"a")
#   writeLines(text=paste0("\n*** ", i, ") ", curr_reg, "(", min_start, "-", max_end, ")"), con=mycon)
#   close(mycon)
#   
#   matching_enhancer_dt <- enhancerDT[enhancerDT$chrom == chromo &
#                                      enhancerDT$start >= min_start &
#                                        enhancerDT$start <= max_end,]
#   
#   
#   matching_enhancer_dt$connected_genes <- sapply(matching_enhancer_dt$attributes, function(x)  paste0(unlist(str_extract_all(x,pattern)), collapse=""))
#   matching_enhancer_dt$connected_genes <- gsub("connected_gene=", "", matching_enhancer_dt$connected_genes)
#   
#   matching_enhancer_dt <- matching_enhancer_dt[,c("chrom", "source", "feature.name", "start", "end", "connected_genes")]
#   
#   
#   write.table(matching_enhancer_dt, file = outfile, append=T, quote=F,col.names = FALSE,row.names=FALSE,sep="\t")
#   
# }
#   
# cat(paste0("... written: ", outfile, "\n"))

######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

