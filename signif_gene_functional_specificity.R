startTime <- Sys.time()
cat(paste0("> Rscript signif_genes_functional_specifity.R\n"))

# from Mazandu and Mulder (2014)
# 2.1.2. Assessing functional specificity of a gene
# A given gene can perform several functions or be involved in several processes. In this case, the gene is annotated by a set of GO terms. 
# The specificity of each term assessing its informativeness depends on its position in the GO structure and the deeper the term is in the DAG structure the more specific or informative the term is. 
# This indicates that the closer or the more similar to the leaf term (term without a child term) the more specific or informative the term is. 
# The specificity score of an annotated gene depends on the specificity of terms used to annotate the gene and is the average of specificity scores of terms in the set of strict non-redundant terms annotating the gene. 
# Thus, for a given gene g, its functional specificity score yesS(g) is assessed by measuring how similar its GO annotations are to leaf terms of the GO DAG connected to the set yesXg of strict non-redundant terms annotating g. 
# This functional specificity score is computed as follows:

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
library(ontologyIndex)
data(go)
suppressPackageStartupMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
require(reshape2)

require(topGO)
require(igraph)
require(ggpubr)


myHeightGG <- 7
myWidthGG <- myHeightGG*1.2

buildTable <- TRUE

printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}

options(scipen=100)

# Rscript signif_genes_functional_specifity.R 

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4


myHeightGG <- 7
myWidthGG <- myHeightGG*1.2

limmaCol <- "#00AFBB"
tadCol <- "#FC4E07"

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


mainFolder <- file.path(".")

pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"

all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds


args <- commandArgs(trailingOnly = TRUE)
tads_signifThresh <- 0.01
tads_signifThresh <- args[1]

genes_signifTresh <- 0.05
genes_signifTresh <- args[2]

file_suffix <- paste0("tadPvalThresh", tads_signifThresh, "_genePvalThresh", genes_signifTresh)

outFolder <- file.path("SIGNIF_GENES_FUNCTIONAL_SPECIFICITY", file_suffix)
dir.create(outFolder, recursive = TRUE)

logFile <- file.path(outFolder, "signif_genes_functional_specifity_logFile.txt")
file.remove(logFile)


final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))
signif_column <- "adjPvalComb"
signifcol <- paste0(signif_column, "_", tads_signifThresh)
final_dt[, paste0(signifcol)] <- final_dt[, paste0(signif_column)] <= tads_signifThresh
cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> tads_signifThresh\t=\t", tads_signifThresh, "\n"))


####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start preparing data before matching \n")
  
  hicds = all_hicds[1]
  all_signif_genes_leafDist <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat(paste0("... start ", hicds, " - ", exprds, "\n"))
      
      dataset_pipDir <- file.path(pipFolder, hicds, exprds)
      stopifnot(dir.exists(dataset_pipDir))
      
      
      
      geneList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      pipeline_geneList <- get(load(geneList_file))
      
      stopifnot(pipeline_geneList %in% g2t_dt$entrezID)
      
      exprds_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% pipeline_geneList,]
      
      ### retrieve genes from signif. tads
      signif_tads <- final_dt$region[final_dt$hicds == hicds &
                                       final_dt$exprds == exprds &
                                       final_dt[,paste0(signifcol)]]
      
      
      tads_signif_genes <- exprds_g2t_dt$entrezID[exprds_g2t_dt$region %in% signif_tads]
      stopifnot(tads_signif_genes %in% pipeline_geneList)
      
      ### retrieve limma signif genes
      topTable_DT_file <- file.path(pipFolder,  hicds, exprds, script1_name, "DE_topTable.Rdata")
      stopifnot(file.exists(topTable_DT_file))
      topTable_DT <- get(load(topTable_DT_file))
      topTable_DT$genes <- as.character(topTable_DT$genes)
      stopifnot(names(pipeline_geneList) %in% topTable_DT$genes)
      topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipeline_geneList),]
      topTable_DT$entrezID <- pipeline_geneList[paste0(topTable_DT$genes)]
      stopifnot(!is.na(topTable_DT$entrezID))
      stopifnot(!duplicated(topTable_DT$entrezID))
      stopifnot(topTable_DT$entrezID %in% pipeline_geneList)
      topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipeline_geneList),]
      topTable_DT$entrezID <- pipeline_geneList[paste0(topTable_DT$genes)]
      stopifnot(!is.na(topTable_DT$entrezID))
      stopifnot(!duplicated(topTable_DT$entrezID))
      
      limma_signif_genes <- topTable_DT$entrezID[topTable_DT$adj.P.Val <= genes_signifTresh]
      stopifnot(limma_signif_genes %in% pipeline_geneList)
      
      
      signif_tad_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% tads_signif_genes]
      stopifnot(length(signif_tad_genes_GO_list) > 0)
      signif_tad_genes_GO_list_filter <- lapply(signif_tad_genes_GO_list, function(x) {
        Filter(function(k) k[["Ontology"]] == ontologyType, x)
      })
      
      signif_tad_genes_GO_list_filter_mostSpec <- lapply(signif_tad_genes_GO_list_filter, function(x) {
        if(length(x) == 0) return(NA)
        stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
        x[[1]][["GOID"]] # in the list, the 1st is the most specific
      })
      
      
      # retrieve max distance to leave
      # signif_tad_genes_GO_list_filter_mostSpec_maxDistLeaf <- unlist(sapply(unique(signif_tad_genes_GO_list_filter_mostSpec), function(curr_go){
      #   curr_children <- go$children[[paste0(curr_go)]]
      #   GO_g <- makeGOGraph(ont = tolower(ontologyType)) # AnnotationDbi (GOdata from topGO: no intersect nodes and topTADs_genes)
      #   if(!curr_go %in% GO_g@nodes) return(NA)
      #   if(! any(curr_children %in% GO_g@nodes)) return(NA)
      #   curr_children <- curr_children[curr_children %in% GO_g@nodes]
      #   sub_graph <- inducedGraph(dag=GO_g, startNodes=c(curr_go, curr_children))  
      #   #plot(sub_graph)
      #   # !!! NEED TO REVERSE THE GRAPH !!!
      #   sub_graph_rev <- reverseArch(sub_graph) # topGO
      #   # plot(sub_graph_rev)
      #   sub_graph_rev_igraph <- igraph.from.graphNEL(graphNEL = sub_graph_rev)
      #   stopifnot(is.connected(sub_graph_rev_igraph))
      #   stopifnot(is.directed(sub_graph_rev_igraph))
      #   stopifnot(is.dag(sub_graph_rev_igraph))
      #   # retrieve index of the leaves
      #   leaves_idx <- which(sapply(sapply(V(sub_graph_rev_igraph),
      #                                     function(x) neighbors(sub_graph_rev_igraph, x, mode="out")), length) == 0)
      #   leavesGO <- names(V(sub_graph_rev_igraph)[leaves_idx])
      #   nLeaves <- length(leavesGO)
      #   all_leaves_dist <- distances(sub_graph_rev_igraph, v = leavesGO, to = curr_go)
      #   max_leaf_dist <- max(as.numeric(all_leaves_dist))
      #   max_leaf_dist
      # }))
      signif_tad_genes_GO_list_filter_mostSpec_maxDistLeaf <- foreach(curr_go =  unique(signif_tad_genes_GO_list_filter_mostSpec), .combine='c') %do% {
        curr_children <- go$children[[paste0(curr_go)]]
        GO_g <- makeGOGraph(ont = tolower(ontologyType)) # AnnotationDbi (GOdata from topGO: no intersect nodes and topTADs_genes)
        if(!curr_go %in% GO_g@nodes) return(NA)
        if(! any(curr_children %in% GO_g@nodes)) return(NA)
        curr_children <- curr_children[curr_children %in% GO_g@nodes]
        sub_graph <- inducedGraph(dag=GO_g, startNodes=c(curr_go, curr_children))  
        #plot(sub_graph)
        # !!! NEED TO REVERSE THE GRAPH !!!
        sub_graph_rev <- reverseArch(sub_graph) # topGO
        # plot(sub_graph_rev)
        sub_graph_rev_igraph <- igraph.from.graphNEL(graphNEL = sub_graph_rev)
        stopifnot(is.connected(sub_graph_rev_igraph))
        stopifnot(is.directed(sub_graph_rev_igraph))
        stopifnot(is.dag(sub_graph_rev_igraph))
        # retrieve index of the leaves
        leaves_idx <- which(sapply(sapply(V(sub_graph_rev_igraph),
                                          function(x) neighbors(sub_graph_rev_igraph, x, mode="out")), length) == 0)
        leavesGO <- names(V(sub_graph_rev_igraph)[leaves_idx])
        nLeaves <- length(leavesGO)
        all_leaves_dist <- distances(sub_graph_rev_igraph, v = leavesGO, to = curr_go)
        max_leaf_dist <- max(as.numeric(all_leaves_dist))
        max_leaf_dist
      }
            
      
      
      
      signif_limma_genes_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% limma_signif_genes]
      stopifnot(length(signif_limma_genes_GO_list) > 0)
      signif_limma_genes_GO_list_filter <- lapply(signif_limma_genes_GO_list, function(x) {
        Filter(function(k) k[["Ontology"]] == ontologyType, x)
      })
      signif_limma_genes_GO_list_filter_mostSpec <- lapply(signif_limma_genes_GO_list_filter, function(x) {
        if(length(x) == 0) return(NA)
        stopifnot(  unlist(lapply(x, function(k) k[["Ontology"]] == ontologyType )))
        x[[1]][["GOID"]] # in the list, the 1st is the most specific
      })
      
      
      
      
      # retrieve max distance to leave
      # signif_limma_genes_GO_list_filter_mostSpec_maxDistLeaf <- unlist(sapply(unique(signif_limma_genes_GO_list_filter_mostSpec), function(curr_go){
      #   curr_children <- go$children[[paste0(curr_go)]]
      #   GO_g <- makeGOGraph(ont = tolower(ontologyType)) # AnnotationDbi (GOdata from topGO: no intersect nodes and topTADs_genes)
      #   if(!curr_go %in% GO_g@nodes) return(NA)
      #   if(! any(curr_children %in% GO_g@nodes)) return(NA)
      #   curr_children <- curr_children[curr_children %in% GO_g@nodes]
      #   sub_graph <- inducedGraph(dag=GO_g, startNodes=c(curr_go, curr_children))
      #   #plot(sub_graph)
      #   # !!! NEED TO REVERSE THE GRAPH !!!
      #   sub_graph_rev <- reverseArch(sub_graph) # topGO
      #   # plot(sub_graph_rev)
      #   sub_graph_rev_igraph <- igraph.from.graphNEL(graphNEL = sub_graph_rev)
      #   stopifnot(is.connected(sub_graph_rev_igraph))
      #   stopifnot(is.directed(sub_graph_rev_igraph))
      #   stopifnot(is.dag(sub_graph_rev_igraph))
      #   # retrieve index of the leaves
      #   leaves_idx <- which(sapply(sapply(V(sub_graph_rev_igraph),
      #                                     function(x) neighbors(sub_graph_rev_igraph, x, mode="out")), length) == 0)
      #   leavesGO <- names(V(sub_graph_rev_igraph)[leaves_idx])
      #   nLeaves <- length(leavesGO)
      #   all_leaves_dist <- distances(sub_graph_rev_igraph, v = leavesGO, to = curr_go)
      #   max_leaf_dist <- max(as.numeric(all_leaves_dist))
      #   max_leaf_dist
      # }))

      signif_limma_genes_GO_list_filter_mostSpec_maxDistLeaf <- foreach(curr_go = unique(signif_limma_genes_GO_list_filter_mostSpec), .combine='c') %do% {
        curr_children <- go$children[[paste0(curr_go)]]
        GO_g <- makeGOGraph(ont = tolower(ontologyType)) # AnnotationDbi (GOdata from topGO: no intersect nodes and topTADs_genes)
        if(!curr_go %in% GO_g@nodes) return(NA)
        if(! any(curr_children %in% GO_g@nodes)) return(NA)
        curr_children <- curr_children[curr_children %in% GO_g@nodes]
        sub_graph <- inducedGraph(dag=GO_g, startNodes=c(curr_go, curr_children))
        #plot(sub_graph)
        # !!! NEED TO REVERSE THE GRAPH !!!
        sub_graph_rev <- reverseArch(sub_graph) # topGO
        # plot(sub_graph_rev)
        sub_graph_rev_igraph <- igraph.from.graphNEL(graphNEL = sub_graph_rev)
        stopifnot(is.connected(sub_graph_rev_igraph))
        stopifnot(is.directed(sub_graph_rev_igraph))
        stopifnot(is.dag(sub_graph_rev_igraph))
        # retrieve index of the leaves
        leaves_idx <- which(sapply(sapply(V(sub_graph_rev_igraph),
                                          function(x) neighbors(sub_graph_rev_igraph, x, mode="out")), length) == 0)
        leavesGO <- names(V(sub_graph_rev_igraph)[leaves_idx])
        nLeaves <- length(leavesGO)
        all_leaves_dist <- distances(sub_graph_rev_igraph, v = leavesGO, to = curr_go)
        max_leaf_dist <- max(as.numeric(all_leaves_dist))
        max_leaf_dist
      }
      
      
      
      
      out_list <-  list(    
                 signif_tad_genes_GO_list_filter_mostSpec_maxDistLeaf=signif_tad_genes_GO_list_filter_mostSpec_maxDistLeaf,
                 signif_limma_genes_GO_list_filter_mostSpec_maxDistLeaf=signif_limma_genes_GO_list_filter_mostSpec_maxDistLeaf
                 )
      
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_out_list.Rdata"))
      save(out_list, file=outFile)
      cat(paste0("... written: ", outFile,"\n"))
      
      
      out_list
      
    }# end-foreach iterating over exprds
    
    exprds_dt
    
  }# end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, paste0("all_signif_genes_leafDist.Rdata"))
  save(all_signif_genes_leafDist, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0( "all_signif_genes_leafDist.Rdata"))
  cat("... load data\n")
  all_signif_genes_leafDist <- get(load(outFile))
}


inFile <- "SIGNIF_GENES_FUNCTIONAL_SPECIFICITY/tadPvalThresh0.01_genePvalThresh0.05/Barutcu_MCF-10A_40kb_TCGAbrca_lum_bas_out_list.Rdata"
inlist <- get(load(inFile))
all_signif_genes_leafDist_ul <- list(inlist)


all_signif_genes_leafDist_ul <- unlist(all_signif_genes_leafDist, recursive = FALSE)

plot_dt <- rbind(
  data.frame(
    signif_type="limma",
    ic_values=unlist(lapply(all_signif_genes_leafDist_ul, function(x) x[["signif_limma_genes_GO_list_filter_mostSpec_maxDistLeaf"]])),
    stringsAsFactors=FALSE),
  data.frame(
    signif_type="tad",
    ic_values=unlist(lapply(all_signif_genes_leafDist_ul, function(x) x[["signif_tad_genes_GO_list_filter_mostSpec_maxDistLeaf"]])),
    stringsAsFactors=FALSE)    
)

nLimma <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="limma"]))
nTAD <- sum(!is.na(plot_dt$ic_values[plot_dt$signif_type=="tad"]))

p <- ggdensity(plot_dt, 
               title = paste0("signif. gene max dist leaf (funct. specificity)"),
               subtitle=paste0("# values: TAD=", nTAD,"; limma=", nLimma),
               x = "ic_values", 
               color = "signif_type", fill = "signif_type",
               # add = "mean", 
               rug = TRUE,
               xlab = "max dist. to leaf",
               palette = c(limmaCol, tadCol))

outFile <- file.path(outFolder, paste0("all_gene_maxDistLeaf_values_", file_suffix, "_density.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile,"\n"))


p <- ggviolin(plot_dt, 
              title = paste0("signif. gene max dist leaf (funct. specificity)"),
              subtitle=paste0("# values: TAD=", nTAD,"; limma=", nLimma),
               y = "max dist. to leaf", 
               x="signif_type",
               color = "signif_type", fill = "signif_type",
               # add = "mean", rug = TRUE,
               xlab = "",
               palette = c(limmaCol, tadCol))

outFile <- file.path(outFolder, paste0("all_gene_maxDistLeaf_values_", file_suffix, "_violin.", plotType))
ggsave(p, file = outFile, height=myHeightGG, width=myHeightGG)
cat(paste0("... written: ", outFile,"\n"))





######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))






