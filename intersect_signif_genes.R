
options(scipen=100)

# Rscript intersect_genes.R

script_name <- "intersect_genes.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <-  ifelse(plotType=="png", 600, 10)
axisCex <- 1.4

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

outFolder <- "INTERSECT_SIGNIF_GENES"
dir.create(outFolder, recursive = TRUE)


infile <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(infile))

all_result_dt <- get(load(infile))

pvalThresh <- 0.01

all_result_dt$signifAdjPvalComb <- all_result_dt$adjPvalComb <= pvalThresh

nDS <- length(unique(paste0(all_result_dt$hicds, all_result_dt$exprds)))


all_vars <- c("signifAdjPvalComb", "signifFDR_0.1", "signifFDR_0.2")

nTop <- 20

signif_var=all_vars[1]


logFile <- file.path(outFolder, paste0("conserv_signif_genes.txt"))
file.remove(logFile)

######################################################################### ALL COMPARISON TYPES
all_data <- list()

for(signif_var in all_vars){
  intersectGenes_dt <- all_result_dt$region_genes[all_result_dt[,paste0(signif_var)]]
  list_intersectGenes <- lapply(intersectGenes_dt, function(str_genes) unlist(strsplit(x=str_genes, split=",")))
  
  conserv_signif <- setNames(
    as.numeric(table(unlist(list_intersectGenes)))/nDS,
    names(table(unlist(list_intersectGenes))))  
  
  outFile <- file.path(outFolder, paste0(signif_var, "_all_cmps_ratio_conserv_signif_genes.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(density(conserv_signif*100), 
       main=paste0("conserv. signif. genes - ", signif_var),
       xlab="% gene conserv. across DS")
  mtext(side=3, text=paste0("all cmps - nDS = ", nDS))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  mostConservedNames <- names(conserv_signif[order(conserv_signif, decreasing=TRUE)])[1:20]
  mostConservedPct <- round(conserv_signif[order(conserv_signif, decreasing=TRUE)]*100, 2)[1:20]
  
  sink(logFile, append=TRUE)
  cat(paste0("> ", signif_var, " - all_cmps\n"))
  cat(paste0("", paste0(paste0(mostConservedNames, " (", mostConservedPct, "%)"), collapse = "\n"), "\n"))
  cat("\n")
  sink()
  # all_data[[length(all_data)+1]] <- conserv_signif*100
    all_data <- c(all_data, list(conserv_signif*100))
}

names(all_data) <- all_vars
outFile <- file.path(outFolder, paste0("all_cmps", "_", "all_signif", "_ratio_conserv_signif_genes.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(all_data, 
     plotTit=paste0("conserv. signif. genes"),
     my_xlab ="% gene conserv. across DS")
mtext(side=3, text=paste0("all cmps - nDS = ", nDS))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


######################################################################### ALL COMPARISON TYPES SEPARATELY

all_result_dt$exprds_type <- all_cmps[all_result_dt$exprds]

cmp_types <- unique(all_result_dt$exprds_type)

all_result_dt_save <- all_result_dt

for(cmp_type in cmp_types) {
  
  all_data <- list()
  
  all_result_dt <- all_result_dt_save[all_result_dt_save$exprds_type == cmp_type,]
  
  nDS <- length(unique(paste0(all_result_dt$hicds, all_result_dt$exprds)))
  
  for(signif_var in all_vars){
    intersectGenes_dt <- all_result_dt$region_genes[all_result_dt[,paste0(signif_var)]]
    list_intersectGenes <- lapply(intersectGenes_dt, function(str_genes) unlist(strsplit(x=str_genes, split=",")))
    
    conserv_signif <- setNames(
      as.numeric(table(unlist(list_intersectGenes)))/nDS,
      names(table(unlist(list_intersectGenes))))  
    
    outFile <- file.path(outFolder, paste0(cmp_type, "_", signif_var, "_ratio_conserv_signif_genes.", plotType))
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    plot(density(conserv_signif*100), 
         main=paste0("conserv. signif. genes - ", signif_var),
         xlab="% gene conserv. across DS")
    mtext(side=3, text=paste0(cmp_type, " - nDS = ", nDS))
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    mostConservedNames <- names(conserv_signif[order(conserv_signif, decreasing=TRUE)])[1:20]
    mostConservedPct <- round(conserv_signif[order(conserv_signif, decreasing=TRUE)]*100, 2)[1:20]
    
    sink(logFile, append=TRUE)
    cat(paste0("> ", signif_var, " - ", cmp_type, "\n"))
    cat(paste0("", paste0(paste0(mostConservedNames, " (", mostConservedPct, "%)"), collapse = "\n"), "\n"))
    cat("\n")
    sink()
    # all_data[[length(all_data)+1]] <- conserv_signif*100
    all_data <- c(all_data, list(conserv_signif*100))
  }
  
  names(all_data) <- all_vars
  outFile <- file.path(outFolder, paste0(cmp_type, "_all_signif", "_ratio_conserv_signif_genes.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot_multiDens(all_data, 
                 plotTit=paste0("conserv. signif. genes"),
                 my_xlab ="% gene conserv. across DS")
  mtext(side=3, text=paste0(cmp_type, " - nDS = ", nDS))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
}







