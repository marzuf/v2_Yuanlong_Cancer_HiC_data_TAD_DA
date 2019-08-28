options(scipen=100)

setDir=""

# Rscript intersect_FDR_combPval_lolliPlot_conserv.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript intersect_FDR_combPval_lolliPlot_conserv.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich
# Rscript intersect_FDR_combPval_lolliPlot_conserv.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"
hicds="ENCSR079VIJ_G401_40kb"
exprds="TCGAkich_norm_kich"


script_name <- "intersect_FDR_combPval_lolliPlot_conserv.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

buildData <- TRUE
separateHeatmap <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

require(reshape2)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("plot_lolliTAD_funct.R")

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script19_name <- "19onlyFC_SAM_emp_measurement"
script19sameNbr_name <- "19sameNbr_SAM_emp_measurement"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "INTERSECT_FDR_COMBPVAL_LOLLIPLOT_CONSERV"
dir.create(outFolder, recursive=TRUE)


pvalThresh <- 0.01
matching_ratio_thresh <- 0.8

fdr <- 0.2
pval <- 0.01


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

final_table_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

final_table_DT[, paste0("signifEmpPval_", pvalThresh)] <- final_table_DT$adjPvalComb <= pvalThresh

inFolder <- "TAD_MATCHING_ACROSS_HICDS"
inFile <- file.path(inFolder, "all_matching_dt.Rdata")
stopifnot(file.exists(inFile))

all_matching_dt <- get(load(inFile))
all_matching_dt$overlapBp_ratio <- all_matching_dt$overlapBP/all_matching_dt$ref_totBp
all_matching_dt$overlapGenes_ratio <- all_matching_dt$nOverlapGenes/all_matching_dt$ref_nGenes
stopifnot(na.omit(all_matching_dt$overlapBp_ratio) >= 0 & na.omit(all_matching_dt$overlapBp_ratio) <= 1)
stopifnot(na.omit(all_matching_dt$overlapGenes_ratio) >= 0 & na.omit(all_matching_dt$overlapGenes_ratio <= 1))


all_matching_dt$overlapBp_overMatchRatio <- all_matching_dt$overlapBp_ratio >= matching_ratio_thresh
all_matching_dt$overlapGenes_overMatchRatio <- all_matching_dt$overlapGenes_ratio >= matching_ratio_thresh

nDS <- length(unique(all_matching_dt$ref_dataset))

nOverThreshBp_dt <- aggregate(overlapBp_overMatchRatio ~ ref_dataset + refID, data = all_matching_dt, FUN=sum)
stopifnot(range(nOverThreshBp_dt$overlapBp_overMatchRatio) <= nDS)

nOverThreshGenes_dt <- aggregate(overlapGenes_overMatchRatio ~ ref_dataset + refID, data = all_matching_dt, FUN=sum)
stopifnot(range(nOverThreshGenes_dt$overlapGenes_overMatchRatio) <= nDS)

nOverThresh_dt <- merge(nOverThreshGenes_dt, nOverThreshBp_dt,by=c("ref_dataset", "refID"), all=TRUE)
nOverThresh_dt$hicds <- dirname(nOverThresh_dt$ref_dataset)
nOverThresh_dt$exprds <- basename(nOverThresh_dt$ref_dataset)
nOverThresh_dt$ref_dataset <- NULL

colnames(nOverThresh_dt)[colnames(nOverThresh_dt) == "refID"] <- "region"

sub_final_dt <- final_table_DT[,c("hicds", "exprds","region", "signifFDR_0.2", "signifEmpPval_0.01")]

plot_dt <- merge(sub_final_dt, nOverThresh_dt, by=c("hicds", "exprds", "region"))

plotAll <- FALSE

############# > PLOT FOR EACH DATASET SEPARATELY

if(plotAll) {
  allDS_intersect_DT <- foreach(hicds = all_hicds) %dopar% {
    exprds_intersect_DT <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      cat("... start ", hicds, " - " , exprds, "\n")
      ds_dt <- plot_dt[plot_dt$hicds == hicds & plot_dt$exprds == exprds,]
      stopifnot(nrow(ds_dt) > 0)
      signif_DT <- ds_dt[ds_dt$signifFDR_0.2 & ds_dt$signifEmpPval_0.01,]
      signif_DT <- signif_DT[order(signif_DT$overlapGenes_overMatchRatio, signif_DT$overlapBp_overMatchRatio, decreasing=TRUE),]
      
      if(nrow(signif_DT) > 10) signif_DT <- signif_DT[1:10,]
      
      if(nrow(signif_DT) == 0) return(NULL)
      
      nPlotted <- nrow(signif_DT)
      
      outHeightGG <- min(c(7 * nPlotted/2, 49))
      outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
      outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
      
      plotList <- list()
      for(i_tad in 1:nPlotted) {
        
        mytit <- paste0( signif_DT$region[i_tad], "\n", ">= ", matching_ratio_thresh*100, "% nGenes: ",  signif_DT$overlapGenes_overMatchRatio[i_tad], "/", nDS, 
                         "; >= ", matching_ratio_thresh*100, "% nBp: ",  signif_DT$overlapBp_overMatchRatio[i_tad], "/", nDS)
        
        
        
        plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = exprds,
                                              hicds = hicds,
                                              all_TADs = signif_DT$region[i_tad],
                                              orderByLolli = "startPos", mytitle=mytit)
      } # end-for iterating over TADs to plot
      outFile <- file.path(outFolder, paste0(hicds, exprds, "_intersect_FDR", fdr, "_adjPvalComb", pval,"_signifTADs", ".", plotType ))
      mytit <- paste0(hicds, " - ", exprds, "\n(FDR<=", fdr, "; pval<=", pval, "; sorted overlap)")
      all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
      outHeightGG <- min(c(7 * nPlotted/2, 49))
      outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
      outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)
      
      ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
      cat("... written: ", outFile, "\n")
      
    } # end-foreach iterating over exprds
  }# end-foreach iterating over hicds => allDS_intersect_DT
  
  
}



############# > PLOT THE OVERALL TOP OVERLAP
all_plot_dt <- plot_dt[plot_dt$signifFDR_0.2 & plot_dt$signifEmpPval_0.01,]

all_plot_dt <- all_plot_dt[order(all_plot_dt$overlapGenes_overMatchRatio, all_plot_dt$overlapBp_overMatchRatio, decreasing=TRUE),]
# nrow(all_plot_dt) # 286


nPlotted <- 10

outHeightGG <- min(c(7 * nPlotted/2, 49))
outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)

plotList <- list()
for(i_tad in 1:10) {
  
  mytit <- paste0( all_plot_dt$hicds[i_tad], " - ", all_plot_dt$exprds[i_tad], "\n",
                   all_plot_dt$region[i_tad], "; ", ">= ", matching_ratio_thresh*100, "% nGenes: ",  all_plot_dt$overlapGenes_overMatchRatio[i_tad], "/", nDS, 
                   "; >= ", matching_ratio_thresh*100, "% nBp: ",  all_plot_dt$overlapBp_overMatchRatio[i_tad], "/", nDS)
  
  
  
  plotList[[i_tad]] <- plot_lolliTAD_ds(exprds = all_plot_dt$exprds[i_tad],
                                        hicds = all_plot_dt$hicds[i_tad],
                                        all_TADs = all_plot_dt$region[i_tad],
                                        orderByLolli = "startPos", mytitle=mytit)
} # end-for iterating over TADs to plot
save(plotList, file="plotList1.Rdata")
outFile <- file.path(outFolder, paste0("all_ds", "_intersect_FDR", fdr, "_adjPvalComb", pval,"_signifTADs", "1-10.", plotType ))
mytit <- paste0("Top 1-10","\n(FDR<=", fdr, "; pval<=", pval, "; sorted overlap)")
all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
outHeightGG <- min(c(7 * nPlotted/2, 49))
outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)

ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
cat("... written: ", outFile, "\n")



plotList <- list()
for(i_tad in 11:20) {
  
  mytit <- paste0( all_plot_dt$hicds[i_tad], " - ", all_plot_dt$exprds[i_tad], "\n",
                   all_plot_dt$region[i_tad], "; ", ">= ", matching_ratio_thresh*100, "% nGenes: ",  all_plot_dt$overlapGenes_overMatchRatio[i_tad], "/", nDS, 
                   "; >= ", matching_ratio_thresh*100, "% nBp: ",  all_plot_dt$overlapBp_overMatchRatio[i_tad], "/", nDS)
  
  
  
  plotList[[i_tad-10]] <- plot_lolliTAD_ds(exprds = all_plot_dt$exprds[i_tad],
                                        hicds = all_plot_dt$hicds[i_tad],
                                        all_TADs = all_plot_dt$region[i_tad],
                                        orderByLolli = "startPos", mytitle=mytit)
} # end-for iterating over TADs to plot
save(plotList, file="plotList2.Rdata")
outFile <- file.path(outFolder, paste0("all_ds", "_intersect_FDR", fdr, "_adjPvalComb", pval,"_signifTADs", "11-20.", plotType ))
mytit <- paste0("Top 11-20","\n(FDR<=", fdr, "; pval<=", pval, "; sorted overlap)")
all_plots <- do.call(grid.arrange, c(plotList,  list(ncol=ifelse(nPlotted == 1, 1, 2), top=textGrob(mytit, gp=gpar(fontsize=20,font=2)))))
outHeightGG <- min(c(7 * nPlotted/2, 49))
outHeightGG <- ifelse(nPlotted < 3, outHeightGG*1.5,outHeightGG)
outWidthGG <- ifelse(nPlotted == 1, 20/2, 20)

ggsave(filename = outFile, all_plots, width=outWidthGG, height = outHeightGG)
cat("... written: ", outFile, "\n")




##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

