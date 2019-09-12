startTime <- Sys.time()
cat(paste0("> Rscript create_sameGO_mostSpec_sortNoDup.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript create_sameGO_mostSpec_sortNoDup.R <hicds> <exprds>
# Rscript create_sameGO_mostSpec_sortNoDup.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS

### UPDATE sortNoDup 30.06.2018
# -> sort rows of the gene coord. table to ensure alphabetical order of the genes !!
#    so that after melt the 1st gene will always be the 1st in alphabetical order
# -> add as.character() in apply !!! use "gene1" etc. instead of 1 index in apply 


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
hicds <- args[1]
exprds <- args[2]

### PREPARE GO DATA
ontologyType <- "BP"

## Bimap interface:
x <- org.Hs.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
all_genes_GO_list <- as.list(x[mapped_genes])

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 40))

outFold <- file.path("CREATE_SAME_GO_SORTNODUP", hicds, exprds)
system(paste0("mkdir -p ", outFold))

pipFolder <- file.path(".", "PIPELINE", "OUTPUT_FOLDER")
dataset_pipDir <- file.path(pipFolder, hicds, exprds)
stopifnot(dir.exists(dataset_pipDir))

  
script0_name <- "0_prepGeneData"
pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))

gene_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% pipeline_geneList]
stopifnot(length(gene_GO_list) > 0)

gene_GO_list_filter <- lapply(gene_GO_list, function(x) {
  Filter(function(k) k[["Ontology"]] == ontologyType, x)
})


gene_GO_list_filter_mostSpec <- lapply(gene_GO_list_filter, function(x) {
  if(length(x) == 0) return(NA)
  stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == "BP" )))
  x[[1]][["GOID"]] # in the list, the 1st is the most specific
})


gene_GO_DT <- data.frame(entrezID = as.character(names(unlist(gene_GO_list_filter_mostSpec))),
                         GO = as.character(unlist(gene_GO_list_filter_mostSpec)),
                         stringsAsFactors = FALSE
)

cat("nrow(gene_GO_DT) = ", nrow(gene_GO_DT), "\n")
cat("nrow(na.omit(gene_GO_DT)) = ", nrow(na.omit(gene_GO_DT)), "\n")
cat("length(unique(gene_GO_DT$GO)) = ", length(unique(gene_GO_DT$GO)), "\n")
cat("sum(is.na(gene_GO_DT$GO)) = ", sum(is.na(gene_GO_DT$GO)), "\n")

gene_GO_DT <- gene_GO_DT[!is.na(gene_GO_DT$GO),]

range(table(gene_GO_DT$GO))
hist(table(gene_GO_DT$GO), breaks=100)

all_gos <- unique(gene_GO_DT$GO)


all_GO_pairs <- foreach(go = all_gos, .combine='rbind') %dopar% {
  go_g2g_dt <- gene_GO_DT[gene_GO_DT$GO == go,]
  if(nrow(go_g2g_dt) == 1) return(NULL)
  # UPDATE 30.06.2018 -> ENSURE AS.CHARACTER + ALPHABETICAL ORDER !!!
  go_g2g_dt$entrezID <- as.character(go_g2g_dt$entrezID)
  go_g2g_dt <- go_g2g_dt[order(go_g2g_dt$entrezID),]
  goDT <- as.data.frame(t(combn(go_g2g_dt$entrezID, m=2)))
  colnames(goDT) <- c("gene1", "gene2")
  goDT$GO <- go
  goDT$gene1 <- as.character(goDT$gene1)
  goDT$gene2 <- as.character(goDT$gene2)
  stopifnot(goDT$gene1 < goDT$gene2)
  goDT
}

outFile <- file.path(outFold, "all_GO_pairs.Rdata")
save(all_GO_pairs, file = outFile)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

