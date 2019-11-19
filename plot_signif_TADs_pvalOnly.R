options(scipen=100)

setDir=""

# Rscript plot_signif_TADs_pvalOnly.R <p_value_threshold> <hicds> <exprds>
# Rscript plot_signif_TADs_pvalOnly.R <p_value_threshold> # to run all datasets in one shot
# Rscript plot_signif_TADs_pvalOnly.R 0.01 K562_40kb TCGAlaml_wt_mutFLT3
# Rscript plot_signif_TADs_pvalOnly.R 0.01 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_signif_TADs_pvalOnly.R 0.01 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_signif_TADs_pvalOnly.R 0.01 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_signif_TADs_pvalOnly.R 0.01 GSE99051_786_O_40kb TCGAkich_norm_kich
# Rscript plot_signif_TADs_pvalOnly.R 0.01 ENCSR TCGAkich_norm_kich

script_name <- "plot_signif_TADs_pvalOnly.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

buildData <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(ggplot2)
require(reshape2)


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("plot_lolliTAD_funct.R")


plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.4

nToPlot <- 10

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script8c_name <- "8cOnlyRatioDownFastSave_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script19_name <- "19onlyFC_SAM_emp_measurement"
script19sameNbr_name <- "19sameNbr_SAM_emp_measurement"
script11same_name <- "11sameNbr_runEmpPvalCombined"

pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "PLOT_SIGNIF_TADS_PVALONLY"
dir.create(outFolder, recursive=TRUE)

twoSidedStouffer <- FALSE

FDRthresh=0.1
pvalThresh=0.01
hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1 | length(args) == 3)
pvalThresh <- as.numeric(args[1])

stopifnot(!is.na(pvalThresh))
stopifnot(pvalThresh >= 0 & pvalThresh <=1 )
hicds <- args[2]
exprds <- args[3]

if(length(args) == 2) {
  all_hicds <- list.files(pipOutFolder)
  all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipOutFolder, x)))
} else{
  all_hicds <- hicds
  all_exprds <- setNames(exprds, hicds)
}

### BUILD SIGNIF ALONG FDR THRESH
cat("... start retrieving FDR signif. TADs\n")

foo1 <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  foo2 <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    ##########> ADJ. COMB PVAL SIGNIF TADS

    # RETRIEVE COMBINED EMP PVAL      
    comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
    stopifnot(file.exists(comb_empPval_file))
    comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
	all_regs <- names(comb_empPval)
    stopifnot(setequal(all_regs, names(comb_empPval)))
    comb_empPval <- comb_empPval[all_regs]
    # ADJUST THE PVAL
    adj_empPval_comb <- p.adjust(comb_empPval, method="BH")
    stopifnot(names(adj_empPval_comb) == all_regs)

    adj_empPval_comb <- sort(adj_empPval_comb, decreasing=FALSE) # => in order to plot the top-ranking; smaller pval = better
    # => SIGNIF TADs FOR THE DESIRED PVAL THRESHOLD
    adjCombPval_signifTADs <- names(adj_empPval_comb)[(adj_empPval_comb <= pvalThresh)]

    if(length(adjCombPval_signifTADs) > 0) {


      # ADDED -> load 1x/dataset
      rnaseqDTFile <- file.path(pipOutFolder, hicds, exprds,script0_name, "rna_rnaseqDT.Rdata")
      initListFile <- file.path(pipOutFolder, hicds, exprds,script0_name, "rna_geneList.Rdata")
      geneListFile <- file.path(pipOutFolder, hicds, exprds,script0_name, "pipeline_geneList.Rdata")
      DE_topTableFile <- file.path(pipOutFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")
      entrezDT_file <- file.path("/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
      gene2tadDT_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
      cat(paste0("... load gene2tad_dt\n"))
      curr_gene2tadDT <-  read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
      cat(paste0("... load entrez2symb_dt\n"))
      curr_entrez2symbDT<- entrez2symb_dt
      cat(paste0("... load DE_topTable\n"))
      curr_DE_topTable <- eval(parse(text = load(DE_topTableFile))) 
      cat(paste0("... load rnaseqDT\n"))
      curr_rnaseqDT <- eval(parse(text = load(rnaseqDTFile))) 
      cat(paste0("... load initList\n"))
      curr_initList <- eval(parse(text = load(initListFile)))	
      cat(paste0("... load geneList\n"))
      curr_geneList <- eval(parse(text = load(geneListFile)))

      nPlotted <- min(c(nToPlot, length(adjCombPval_signifTADs)))
      plot_adjCombPval_signifTADs <-   adjCombPval_signifTADs[1:nPlotted]
      cat("... ", hicds , " - ", exprds, ": signif. TADs, adj. pval. comb. thresh - start plotting (max. ", nToPlot, "; found: ", length(adjCombPval_signifTADs), ")\n")
      plotList <- list()
      for(i_tad in 1:length(plot_adjCombPval_signifTADs)) {
        plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                          hicds = hicds,
                                          all_TADs = plot_adjCombPval_signifTADs[i_tad],

										gene2tadDT=curr_gene2tadDT,entrez2symbDT=curr_entrez2symbDT,
										DE_topTable=curr_DE_topTable,
										rnaseqDT=curr_rnaseqDT, initList=curr_initList, geneList=curr_geneList,

                                          orderByLolli = "startPos")
      } # end-for iterating over TADs to plot
      outFile <- file.path(outFolder, paste0(hicds, exprds, "_adjCombPvalSignifTADs",  "_nToPlot", nToPlot, ".", plotType ))
      mytit <- paste0(hicds, " - ", exprds, " - top ", nToPlot, "\n(adj. comb. pval. <= ", pvalThresh, "; sorted adjPvalComb)")
      all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nToPlot == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      outHeightGG <- min(c(7 * nPlotted/2, 49))
      outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
      ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
      cat("... written: ", outFile, "\n")

          }else{
      # NO SIGNIF TADs FOUND
      cat("... ", hicds , " - ", exprds, ": no signif. adj. pval. comb. TADs found \n")  
    }
  } # end-foreach iterating over exprds
} # end-foreach iterating over hicds





##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

