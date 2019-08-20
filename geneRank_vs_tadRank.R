options(scipen=100)

setDir <- ""

# Rscript geneRank_vs_tadRank.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript geneRank_vs_tadRank.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "geneRank_vs_tadRank.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

tiesMeth <- "min"

buildData <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4


pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "GENERANK_VS_TADRANK"
dir.create(outFolder, recursive=TRUE)


FDRthresh_seq <- seq(from=0.1, to=0.5, by=0.1)
pvalThresh_seq <- seq(from=0.01, to=0.05, by = 0.01)

myHeightGG <- length(pvalThresh_seq)*1.2
myWidthGG <- length(FDRthresh_seq)*1.2

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2 | length(args) == 0)
hicds <- args[1]
exprds <- args[2]

if(length(args) == 0) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}


if(buildData) {
  
  
  
  ### BUILD SIGNIF ALONG FDR THRESH
  cat("... start building signif. along FDR thresh\n")
  all_gene_tad_rank_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_gene_tad_rank_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      
      # RETRIEVE PVAL COMBINED DATA
      
      # RETRIEVE COMBINED EMP PVAL      
      comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      # ADJUST THE PVAL
      adj_combPvals <- p.adjust(comb_empPval, method="BH")




      adj_combPvals_rank <- rank(adj_combPvals, ties = tiesMeth) # best rank = smallest value
      stopifnot(adj_combPvals[names(adj_combPvals_rank)[adj_combPvals_rank == 1]] == min(adj_combPvals))
      tad_rank_DT <- data.frame(region=names(adj_combPvals_rank), TAD_rank = adj_combPvals_rank, stringsAsFactors = FALSE)
      
      
      de_file <- file.path(pipOutFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")
      stopifnot(file.exists(de_file))
      de_DT <- eval(parse(text = load(de_file)))
      de_DT$genes <- as.character(de_DT$genes)
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      geneList_file <- file.path(pipOutFolder, hicds, exprds,script0_name,   "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      geneList <- eval(parse(text = load(geneList_file)))
      stopifnot(names(geneList) %in% de_DT$genes)
      stopifnot(geneList %in% g2t_DT$entrezID)
      de_DT <- de_DT[de_DT$genes %in% names(geneList),]
      stopifnot(nrow(de_DT) == length(geneList))
      de_DT$entrezID <- sapply(de_DT$genes, function(x) as.character(geneList[as.character(x)]))
      stopifnot(!is.na(de_DT$entrezID))
      
      gene_DE <- setNames(de_DT$adj.P.Val, de_DT$entrezID)
      gene_DE_rank <- rank(gene_DE, ties = tiesMeth) # best rank = smallest value
      stopifnot(gene_DE[names(gene_DE_rank)[gene_DE_rank == 1]] == min(gene_DE))
      gene_rank_DT <- data.frame(entrezID=names(gene_DE_rank), gene_rank = gene_DE_rank, stringsAsFactors = FALSE)
      
      g2tFile <- file.path(pipFolder, hicds,  "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      stopifnot(geneList %in% g2t_DT$entrezID)
      g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
      
      stopifnot(setequal(g2t_DT$entrezID, gene_rank_DT$entrezID))
      g2t_rank_DT <- merge(gene_rank_DT, g2t_DT[,c("entrezID", "region")], by="entrezID")
      
      stopifnot(setequal(tad_rank_DT$region, g2t_rank_DT$region))
      
      gene_tad_rank_DT <- merge(tad_rank_DT, g2t_rank_DT, by ="region")
      
      gene_tad_rank_DT$hicds <- hicds
      gene_tad_rank_DT$exprds <- exprds
      gene_tad_rank_DT

    } # end-foreach iterating over exprds
    exprds_gene_tad_rank_dt
  } # end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, "all_gene_tad_rank_dt.Rdata")  
  save(all_gene_tad_rank_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else { # end-if build data
  outFile <- file.path(outFolder, "all_gene_tad_rank_dt.Rdata")  
  all_gene_tad_rank_dt <- eval(parse(text = load(outFile)))
}

all_gene_tad_rank_dt$dataset <- paste0( all_gene_tad_rank_dt$hicds, " - ", all_gene_tad_rank_dt$exprds)
all_datasets <- unique(all_gene_tad_rank_dt$dataset) 
nDS <- length(unique( all_datasets ))

x_var = "TAD_rank"
y_var = "gene_rank"
myx <- all_gene_tad_rank_dt[,paste0(x_var)]
myy <- all_gene_tad_rank_dt[,paste0(y_var)]
outFile <- file.path(outFolder, paste0("geneRank", "_vs_", "tadRank", ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x = myx,
  y = myy,
  cex.axis = axisCex,
  cex.lab = axisCex,
  xlab = paste0(gsub("_", " ", x_var)),
  ylab = paste0(gsub("_", " " , y_var)),
  main = paste0(y_var, " vs. ", x_var)
)
mtext(side=3, text = paste0(nDS, " DS\n(n=", nrow(all_gene_tad_rank_dt),")"), font=3)
#curve(1*x, add=TRUE, lty=2, col="black", lwd=1.5)
addCorr(
  x = myx, y = myy,
  legPos = "topleft", bty="n"
)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

for(ds in all_datasets) {
  sub_dt <- all_gene_tad_rank_dt[all_gene_tad_rank_dt$dataset == ds,]
  myx <- sub_dt[,paste0(x_var)]
  myy <- sub_dt[,paste0(y_var)]
  outFile <- file.path(outFolder, paste0(gsub(" - ", "_", ds), "_", "geneRank", "_vs_", "tadRank", ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x = myx,
    y = myy,
    cex.axis = axisCex,
    cex.lab = axisCex,
    xlab = paste0(gsub("_", " ", x_var)),
    ylab = paste0(gsub("_", " " , y_var)),
    main = paste0(y_var, " vs. ", x_var)
  )
  mtext(side=3, text = paste0(ds, "\n(n=", nrow(sub_dt),")"), font=3, line=-1)
#  curve(1*x, add=TRUE, lty=2, col="black", lwd=1.5)
  addCorr(
    x = myx, y = myy,
    legPos = "topleft", bty="n"
  )
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}





##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

