# library(ontologyIndex)
# data(go)
library(clusterProfiler)

library(ontologySimilarity)
# data(gene_GO_terms)
data(GO_IC)
# 
# beach <- gene_GO_terms[c("LRBA", "LYST", "NBEA", "NBEAL1", "NBEAL2", "NSMAF", "WDFY3", "WDFY4", "WDR81")]
# 
# go$name[beach$LRBA]
# 
# cc <- go$id[go$name == "cellular_component"]
# beach_cc <- lapply(beach, function(x) intersection_with_descendants(go, roots=cc, x)) 
# data.frame(check.names=FALSE, `#terms`=sapply(beach, length), `#CC terms`=sapply(beach_cc, length))
# 
# sim_matrix <- get_sim_grid(
#   ontology=go, 
#   information_content=GO_IC,
#   term_sets=beach)
# 
# library(paintmap)
# paintmap(colour_matrix(sim_matrix))
# 
# get_sim_p_from_ontology(
#   ontology=go,
#   information_content=GO_IC,
#   term_sets=gene_GO_terms,
#   group=names(beach)
# )
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

options(scipen=100)

SSHFS=F

startTime <- Sys.time()

script_name <- "go_specificity_geneLevel_tadLevel.R"

cat("> START ", script_name, "\n")

# Rscript go_specificity_geneLevel_tadLevel.R

buildTable <- TRUE

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- file.path("GO_SPECIFICITY_GENELEVEL_TADLEVEL")
dir.create(outFolder, recursive = TRUE)



logFile <- file.path(outFolder, "go_signif_geneLevel_tadLevel_logFile.txt")
if(buildTable) file.remove(logFile)


printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}

my_GO2TERM <- clusterProfiler:::get_GO2TERM_table()
my_GO2TERM$term <- paste0("GO_", toupper(gsub(" ", "_", my_GO2TERM$Term)))


printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}


go_signif_col <- "p.adjust"
go_signifThresh <- 0.05


txt <- paste0("... go_signif_col\t=\t", go_signif_col, "\n")
printAndLog(txt, logFile)
txt <- paste0("... go_signifThresh\t=\t", go_signifThresh, "\n")
printAndLog(txt, logFile)


inFile <- file.path("GO_SIGNIF_GENELEVEL_TADLEVEL", "all_go_enrich_list.Rdata")
stopifnot(file.exists(inFile))
all_go_enrich_list <- get(load(inFile))

if(buildTable) {
  
  hicds = all_hicds[1]
  all_go_ic_list <- foreach(hicds = all_hicds) %dopar% {
    
  
      
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      
      
      cat(paste0("... start: ", hicds, " - ", exprds, "\n"))
      
      go_signif_limma_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["limma_signif_enrich_resultDT"]]
      
      
      if(!is.null(go_signif_limma_dt)) {
        
        txt <- paste0(hicds, " - ", exprds, " - limma signif.: # annot. GO:\t", nrow(go_signif_limma_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_limma_dt <- go_signif_limma_dt[go_signif_limma_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - limma signif.: # signif. annot. GO:\t", nrow(go_signif_limma_dt), "\n")
        printAndLog(txt, logFile)
        
        signif_limma_go_terms <- rownames(go_signif_limma_dt)
  
        signif_limma_go_ids <- my_GO2TERM$go_id[my_GO2TERM$term %in% signif_limma_go_terms]
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: found GO ids matching GO terms:\t", length(signif_limma_go_ids), "/", length(signif_limma_go_terms), "\n")
        printAndLog(txt, logFile)
        
        signif_limma_go_ic <- GO_IC[names(GO_IC) %in% signif_limma_go_ids]
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: found GO information content matching GO ids:\t", length(signif_limma_go_ic), "/", length(signif_limma_go_ids), "\n")
        printAndLog(txt, logFile)
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: total retrieved data:\t", length(signif_limma_go_terms), " -> ", length(signif_limma_go_ic), "\n")
        printAndLog(txt, logFile)
        
        limma_signif_mean_go_ic <- mean(signif_limma_go_ic)
        
        
      } else {
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: NULL \n")
        printAndLog(txt, logFile)
        
        limma_signif_mean_go_ic <- NULL
        
      }
      
      go_signif_tads_dt <- all_go_enrich_list[[file.path(hicds, exprds)]][["tads_signif_enrich_resultDT"]]
      
      if(!is.null(go_signif_tads_dt)) {
        
        txt <- paste0(hicds, " - ", exprds, " - TAD signif: # annot. GO:\t", nrow(go_signif_tads_dt), "\n")
        printAndLog(txt, logFile)
        go_signif_tads_dt <- go_signif_tads_dt[go_signif_tads_dt[,paste0(go_signif_col)] <= go_signifThresh,]
        txt <- paste0(hicds, " - ", exprds, " - TAD signif: # signif. annot. GO:\t", nrow(go_signif_tads_dt), "\n")
        printAndLog(txt, logFile)
        
        signif_tad_go_terms <- rownames(go_signif_tads_dt)
        
        signif_tad_go_ids <- my_GO2TERM$go_id[my_GO2TERM$term %in% signif_tad_go_terms]
        txt <- paste0(hicds, " - ", exprds, " - ... TAD signif: found GO ids matching GO terms:\t", length(signif_tad_go_ids), "/", length(signif_tad_go_terms), "\n")
        printAndLog(txt, logFile)
        
        signif_tad_go_ic <- GO_IC[names(GO_IC) %in% signif_tad_go_ids]
        txt <- paste0(hicds, " - ", exprds, " - ... TAD signif: found GO information content matching GO ids:\t", length(signif_tad_go_ic), "/", length(signif_tad_go_ids), "\n")
        printAndLog(txt, logFile)
        txt <- paste0(hicds, " - ", exprds, " - ... TAD signif: total retrieved data:\t", length(signif_tad_go_terms), " -> ", length(signif_tad_go_ic), "\n")
        printAndLog(txt, logFile)
        
        tad_signif_mean_go_ic <- mean(signif_tad_go_ic)
        
      } else {
        
        txt <- paste0(hicds, " - ", exprds, " - ... limma signif.: NULL \n")
        printAndLog(txt, logFile)
        
        tad_signif_mean_go_ic <- NULL
        
        
      }
      
      list(
        limma_signif_mean_go_ic = limma_signif_mean_go_ic,
        tad_signif_mean_go_ic = tad_signif_mean_go_ic
      )
      
    } # end-foreach iterating over exprds
    names(exprds_list) <- file.path(hicds, all_exprds[[paste0(hicds)]])
    exprds_list
  } # end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, paste0("all_go_ic_list.Rdata"))
  save(all_go_ic_list, file = outFile, version=2)
  stopifnot(length(all_go_ic_list) == length(all_datasets))
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0( "all_go_ic_list.Rdata"))
  cat("... load data\n")
  all_go_ic_list <- get(load(outFile))
}





######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



