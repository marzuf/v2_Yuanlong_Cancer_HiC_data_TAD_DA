startTime <- Sys.time()
cat(paste0("> Rscript enrich_db_annot.R\n"))

# Rscript tad_matching_enrich_db_annot.R
# Rscript tad_matching_enrich_db_annot.R subtypes
# Rscript tad_matching_enrich_db_annot.R norm_vs_tumor
# Rscript tad_matching_enrich_db_annot.R wt_vs_mut

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)


signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3

annot1_name <- "cancermine"
annot1_file <- file.path("db_gene_annot/cancermine_collated.tsv")
annot1_dt <- read.delim(annot1_file, header=TRUE, stringsAsFactors = FALSE)
annot1_dt$gene_entrez_id <- as.character(annot1_dt$gene_entrez_id)
ambiguous_entrez <- annot1_dt$gene_entrez_id[duplicated(annot1_dt$gene_entrez_id)]
annot1_dt <- annot1_dt[!annot1_dt$gene_entrez_id %in% ambiguous_entrez,]
stopifnot(!duplicated(annot1_dt$gene_entrez_id))
annot1_entrez <- setNames(annot1_dt$role, annot1_dt$gene_entrez_id)

annot2_name <- "cigene"
annot2_file <- file.path("db_gene_annot/cigene_human.txt")
annot2_dt <- read.delim(annot2_file, header=TRUE, stringsAsFactors = FALSE)
all_genes1 <- annot2_dt$GeneSymbol
all_genes2 <- unlist(strsplit(annot2_dt$Alias, split="\\|"))
annot2_symbols <- setNames(rep("cancer_initiation", length(c(all_genes1, all_genes2))), c(all_genes1, all_genes2))

annot3_name <- "cgi"
annot3_file <- file.path("db_gene_annot/cgi-bin/cgi_bin_gene_annot.txt")
annot3_dt <- read.delim(annot3_file, header=FALSE, stringsAsFactors = FALSE, col.names = c("symbol", "role", "description") )
stopifnot( !duplicated(annot3_dt$symbol))
annot3_symbols <- setNames(annot3_dt$role, annot3_dt$symbol)


args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  inFolder <- "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH"
  args <- ""
} else if(length(args) == 1) {
  inFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH", args[1])
  args <- args[1]
} else {
  stop("error\n")
}
stopifnot(dir.exists(inFolder))
inFile <- file.path(inFolder, paste0("conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, "_df.txt"))
stopifnot(file.exists(inFile))
conserved_dt <- read.delim(inFile, header=TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(conserved_dt) > 0)

outFolder <- file.path("TAD_MATCHING_ENRICH_DB_ANNOT", args)
dir.create(outFolder, recursive = TRUE)
outfile <- file.path(outFolder, "tad_matching_enrich_db_annot.txt")
file.remove(outfile)

i=2
for(i in 1:nrow(conserved_dt)) {
  
  # cat("i = ", i, "\n")
  
  curr_reg <- conserved_dt$conserved_region[i]
  
  mycon <- file(outfile,"a")
  writeLines(text=paste0("\n*** ", i, ") ", curr_reg), con=mycon)
  close(mycon)
  
  
  
  curr_symbols <- unlist(strsplit(x=conserved_dt$intersect_genes_symbol[i], split=","))
  curr_entrez <- unlist(strsplit(x=conserved_dt$intersect_genes_entrez[i], split=","))
  stopifnot(length(curr_symbols) == length(curr_entrez))
  
  curr_genes <- setNames(curr_symbols, curr_entrez)
  
    if(any(names(curr_genes) %in% names(annot1_entrez))) {
      
      a1_genes <- curr_genes[names(curr_genes) %in% names(annot1_entrez)]
      
      dt <- data.frame(entrezID = names(a1_genes),
                       symbol = as.character(a1_genes),
                       annot = as.character(annot1_entrez[names(a1_genes)]),
                       stringsAsFactors = FALSE)
      stopifnot(!is.na(dt))
      mycon <- file(outfile,"a")
      writeLines(text=paste0("> ", curr_reg, " - ", annot1_name, " annotation:"), con=mycon)
      close(mycon)
      write.table(dt, file = outfile, append=T, quote=F,col.names = FALSE,row.names=FALSE,sep="\t")
      
    }
    
    
    
    if(any(paste0(curr_genes) %in% names(annot2_symbols))) {
      
      a2_genes <- curr_genes[paste0(curr_genes) %in% names(annot2_symbols)]
      
      dt <- data.frame(entrezID = names(a2_genes),
                       symbol = as.character(a2_genes),
                       annot = as.character(annot2_symbols[paste0(a2_genes)]),
                       stringsAsFactors = FALSE)
      stopifnot(!is.na(dt))
      mycon <- file(outfile,"a")
      writeLines(text=paste0("> ", curr_reg, " - ", annot2_name, " annotation:"), con=mycon)
      close(mycon)
      write.table(dt, file = outfile, append=T, quote=F,col.names = FALSE,row.names=FALSE,sep="\t")
      
    }
    
    if(any(paste0(curr_genes) %in% names(annot3_symbols))) {
      
      a3_genes <- curr_genes[paste0(curr_genes) %in% names(annot3_symbols)]
      
      dt <- data.frame(entrezID = names(a3_genes),
                       symbol = as.character(a3_genes),
                       annot = as.character(annot3_symbols[paste0(a3_genes)]),
                       stringsAsFactors = FALSE)
      stopifnot(!is.na(dt))
      
      mycon <- file(outfile,"a")
      writeLines(text=paste0("> ", curr_reg, " - ",  annot3_name, " annotation:"), con=mycon)
      close(mycon)
      write.table(dt, file = outfile, append=T, quote=F,col.names = FALSE,row.names=FALSE,sep="\t")
      
    }
    
  
}
  
cat(paste0("... written: ", outfile, "\n"))

######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

