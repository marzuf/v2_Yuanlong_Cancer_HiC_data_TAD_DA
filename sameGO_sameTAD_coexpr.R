startTime <- Sys.time()
cat(paste0("> Rscript sameGO_sameTAD_coexpr.R\n"))

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(org.Hs.eg.db, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

options(scipen=100)

# Rscript sameGO_sameTAD_coexpr.R


buildTable <- TRUE

corMethod <- "pearson"

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

outFold <- file.path("SAMEGO_SAMETAD_COEXPR")
dir.create(outFold, recursive = TRUE)

mainFolder <- file.path(".")
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")

all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds


script0_name <- "0_prepGeneData"

all_hicds=all_hicds[1]

####################################################################################################################################### >>> prepare the data
if(buildTable) {
  all_go_coexpr_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    
    tadFile <- file.path("CREATE_SAME_TAD_SORTNODUP", hicds,  "all_TAD_pairs.Rdata")
    stopifnot(file.exists(tadFile))
    all_TAD_pairs <- get(load(tadFile))
    all_TAD_pairs$gene1 <- as.character(all_TAD_pairs$gene1)
    all_TAD_pairs$gene2 <- as.character(all_TAD_pairs$gene2)
    stopifnot(all_TAD_pairs$gene1 < all_TAD_pairs$gene2)
    
    all_TAD_pairs_init <- all_TAD_pairs
    rm(all_TAD_pairs)
    
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, "\n")
      
      # > PREPARE geneList DATA
      dataset_pipDir <- file.path(pipFolder, hicds, exprds)
      stopifnot(dir.exists(dataset_pipDir))
      pipeline_geneList <- eval(parse(text = load(file.path(dataset_pipDir, script0_name, "pipeline_geneList.Rdata"))))
      

      # > PREPARE sameTAD DATA
      # ! need to be done because TAD dt built at hicds level (not hicds/exprds level)
      all_TAD_pairs <- all_TAD_pairs_init[all_TAD_pairs_init$gene1 %in% pipeline_geneList & 
                                            all_TAD_pairs_init$gene2 %in% pipeline_geneList,]
      stopifnot(pipeline_geneList %in% all_TAD_pairs$gene1 | pipeline_geneList %in% all_TAD_pairs$gene2)      
            
      # > PREPARE GO DATA
      gene_GO_list <- all_genes_GO_list[names(all_genes_GO_list) %in% pipeline_geneList]
      stopifnot(length(gene_GO_list) > 0)
      gene_GO_list_filter <- lapply(gene_GO_list, function(x) {
        Filter(function(k) k[["Ontology"]] == ontologyType, x)
      })
      all_go_list <- lapply(gene_GO_list_filter, names)
      all_go_dt <- data.frame(
        go = as.character(unlist(all_go_list)),
        entrezID = rep(names(all_go_list), lengths(all_go_list) ),
        stringsAsFactors = FALSE
      )
      stopifnot(nrow(all_go_dt) > 0)
      stopifnot(setequal(all_go_dt$go[all_go_dt$entrezID==all_go_dt$entrezID[1]], all_go_list[[all_go_dt$entrezID[1]]]))
      stopifnot(setequal(all_go_dt$go[all_go_dt$entrezID==all_go_dt$entrezID[nrow(all_go_dt)]], all_go_list[[all_go_dt$entrezID[nrow(all_go_dt)]]]))
      
      stopifnot(all_go_dt$entrezID %in% all_TAD_pairs$gene1 | all_go_dt$entrezID %in% all_TAD_pairs$gene2)
      stopifnot(all_go_dt$entrezID %in% pipeline_geneList)
      
      # > PREPARE coexpr DATA
      coexprFile <- file.path("CREATE_COEXPR_SORTNODUP",hicds, exprds, corMethod, "coexprDT.Rdata")
      stopifnot(file.exists(coexprFile))
      coexprDT <- get(load(coexprFile))
      coexprDT$gene1 <- as.character(coexprDT$gene1)
      coexprDT$gene2 <- as.character(coexprDT$gene2)
      stopifnot(pipeline_geneList %in% coexprDT$gene1 | pipeline_geneList %in% coexprDT$gene2)
      stopifnot(all_go_dt$entrezID %in% coexprDT$gene1 | all_go_dt$entrezID %in% coexprDT$gene2)
      stopifnot(coexprDT$gene1 < coexprDT$gene2)
      coexprDT <- coexprDT[coexprDT$gene1 %in% pipeline_geneList & 
                                       coexprDT$gene2 %in% pipeline_geneList,]
      i=1
      all_go_coexpr_dt <- foreach(i = 1:nrow(all_go_dt), .combine='rbind') %dopar% {
        
        go <- all_go_dt$go[i]
        go_genes <- all_go_dt$entrezID[all_go_dt$go == go]
        
        stopifnot(go_genes %in% pipeline_geneList)
        
        
        go_coexpr_dt <- coexprDT[coexprDT$gene1 %in% go_genes | coexprDT$gene2 %in% go_genes,]
      
        go_tad_dt <- all_TAD_pairs[all_TAD_pairs$gene1 %in% go_genes | all_TAD_pairs$gene2 %in% go_genes,]
        
        
        merge_dt <- merge(go_coexpr_dt, go_tad_dt, all=TRUE, by=c("gene1", "gene2"))
        
        stopifnot(!is.na(merge_dt$coexpr))
        
        merge_dt$mapTAD <- ifelse(is.na(merge_dt$region), "diffTAD", "sameTAD")
        stopifnot(merge_dt$mapTAD[is.na(merge_dt$region)] == "diffTAD")
        
        mapTAD_meanCoexpr <- aggregate(coexpr ~ mapTAD, data = merge_dt, mean)
        colnames(mapTAD_meanCoexpr)[colnames(mapTAD_meanCoexpr) == "coexpr"] <- "meanCoexpr"
        
        mapTAD_nPairs <- aggregate(coexpr ~ mapTAD, data = merge_dt, length)
        colnames(mapTAD_nPairs)[colnames(mapTAD_nPairs) == "coexpr"] <- "nPairs"
        
        go_dt <- merge(mapTAD_meanCoexpr, mapTAD_nPairs, by="mapTAD")
        otherCols <- colnames(go_dt)
        firstCols <- c("hicds", "exprds", "GO")
        go_dt$hicds <- hicds
        go_dt$exprds <- exprds
        go_dt$GO <- go
        go_dt[,c(firstCols, otherCols)]
      } # end-foreach iterating over GO categories
    all_go_coexpr_dt
    } # end-foreach iterating over exprds
    exprds_dt
  }# end-foreach iterating over hicds

  outFile <- file.path(outFolder, paste0("all_go_coexpr_dt.Rdata"))
  save(all_go_coexpr_dt, file=outFile, version=2)  
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, paste0("all_go_coexpr_dt.Rdata"))
  cat("... load data\n")
  all_go_coexpr_dt <- get(load(outFile))
}

######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

