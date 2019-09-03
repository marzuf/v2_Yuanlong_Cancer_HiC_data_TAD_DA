startTime <- Sys.time()
cat(paste0("> Rscript enrich_db_annot.R\n"))

# Rscript enrich_db_annot.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS
# Rscript enrich_db_annot.R


options(scipen=100)


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

buildData <- TRUE

script0_name <- "0_prepGeneData"

outFolder <- file.path("ENRICH_DB_ANNOT")
dir.create(outFolder, recursive = TRUE)

pipFolder <- file.path(".")

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

annot_file <- file.path("db_gene_annot/cancermine_collated.tsv")
annot_dt <- read.delim(annot_file, header=TRUE)
annot_dt$gene_entrez_id <- as.character(annot_dt$gene_entrez_id)
ambiguous_entrez <- annot_dt$gene_entrez_id[duplicated(annot_dt$gene_entrez_id)]
annot_dt <- annot_dt[!annot_dt$gene_entrez_id %in% ambiguous_entrez,]
stopifnot(!duplicated(annot_dt$gene_entrez_id))
annot_entrez <- setNames(annot_dt$role, annot_dt$gene_entrez_id)



hicds="Panc1_rep12_40kb"
exprds="TCGApaad_wt_mutKRAS"

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}


final_dt_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))



if(buildData) {
  
  
  all_annot_result <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    
    exprds_data <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      
      stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
      geneListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneListFile))
      pipeline_geneList <- eval(parse(text = load(geneListFile))) 
      
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
      g2t_DT <- g2t_DT[ g2t_DT$entrezID %in% pipeline_geneList,]
      
      g2t_DT$annot <- annot_entrez[g2t_DT$entrezID]
      
      
      n_tumorSuppressor_dt <- aggregate(annot ~  region, data = g2t_DT, function(x) sum(na.omit(x) == "Tumor_Suppressor"))
      n_tumorSuppressor <- setNames(n_tumorSuppressor_dt$annot, n_tumorSuppressor_dt$region)
      
      n_Driver_dt <- aggregate(annot ~  region, data = g2t_DT, function(x) sum(na.omit(x) == "Driver"))
      n_Driver <- setNames(n_Driver_dt$annot, n_Driver_dt$region)
      
      n_Oncogene_dt <- aggregate(annot ~  region, data = g2t_DT, function(x) sum(na.omit(x) == "Oncogene"))
      n_Oncogene <- setNames(n_Oncogene_dt$annot, n_Oncogene_dt$region)
      
      curr_final <- final_dt[final_dt$hicds == hicds &
                               final_dt$exprds == exprds,
                               ]
      stopifnot(curr_final$region %in% g2t_DT$region)
      
      curr_final <- curr_final[,c("hicds", "exprds", "region", "adjPvalComb","signifFDR_0.1", "signifFDR_0.2")]
      curr_final$nTumorSuppressor <- n_tumorSuppressor[curr_final$region]
      curr_final$nTumorSuppressor[is.na(curr_final$nTumorSuppressor)] <- 0
      curr_final$nDriver <- n_Driver[curr_final$region]
      curr_final$nDriver[is.na(curr_final$nDriver)] <- 0
      curr_final$nOncogene <- n_Oncogene[curr_final$region]
      curr_final$nOncogene[is.na(curr_final$nOncogene)] <- 0
      
      
      curr_final
      
      
    } # end-for iterating over hicds
    exprds_data
    
  } # end-foreach iterating over exprcds
  
  

    
  outFile <- file.path(outFolder, "all_annot_result.Rdata")
  save(all_annot_result, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_annot_result.Rdata")  
  all_annot_result <- get(load(outFile))
}



######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

