SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "GSE36401_mRNA_topTADs_plot_vEntrez_v2.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

library(foreach)
library(doMC)
library(ggpubr)



registerDoMC(ifelse(SSHFS, 20, 2))

# Rscript GSE36401_mRNA_topTADs_plot_vEntrez_v2.R Rao_HCT-116_2017_40kb
# Rscript GSE36401_mRNA_topTADs_plot_vEntrez_v2.R ENCSR504OTV_transverse_colon_40kb_TCG
# Rscript GSE36401_mRNA_topTADs_plot_vEntrez_v2.R GSE105318_DLD1_40kb

plot_countType <- c("", "_log10")

plotType <- "png"
myHeight <- 7
myWidth <- 10

args <- commandArgs(trailingOnly = TRUE)
curr_dataset = "Rao_HCT-116_2017_40kb"
curr_dataset <- args[1]
cond1 <- "MSI"
cond2 <- "MSS"

exprds <- "TCGAcoad_msi_mss"

# HARD CODED, they should be 2 MSI and 4 MSS !
check_nCond1 <- 4
check_nCond2 <- 5

# outFold <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data/GSE36401_Akhtar-Zaidi2012/GSE36401_mRNA_TOPTADs_BOXPLOT_vEntrez_v2", curr_dataset)
outFold <- file.path("GSE36401_mRNA_TOPTADs_BOXPLOT_vEntrez_v2", curr_dataset)

dir.create(outFold, recursive = TRUE)

logFile <- file.path(outFold, "mRNA_topTADs_boxplot_logFile.txt")
system(paste0("rm -f ", logFile))


# /mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/IGV_PLOTS/ENCSR504OTV_transverse_colon_40kb_TCGAcoad_msi_mss_TADcoord.bed
# /mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/IGV_PLOTS/GSE105318_DLD1_40kb_TCGAcoad_msi_mss_TADcoord.bed
# 
topTADfile <- file.path(setDir, "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/IGV_PLOTS", paste0(curr_dataset, "_", exprds, "_TADcoord.bed"))
stopifnot(file.exists(topTADfile))
topTADdt <- read.delim(topTADfile, header=F, sep="\t", col.names=c("chromo", "start", "end", "region"), stringsAsFactors = F)
stopifnot(is.numeric(topTADdt$start))
stopifnot(is.numeric(topTADdt$end))

mRNAfile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_coreg/h3k27ac_msi_mss_data/GSE_data",
                        "GSE36401_Akhtar-Zaidi2012/mRNA_DATA_entrezID/TCGAcrc_msi_mss/mRNA_with_region_DT_vEntrez.Rdata")
stopifnot(file.exists(mRNAfile))
mrnaDT <- eval(parse(text = load(mRNAfile)))
mrnaDT$region <- NULL
stopifnot(is.numeric(mrnaDT$start))
stopifnot(is.numeric(mrnaDT$end))

pipFolder <- file.path(".")
pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")


tadDT_file <- file.path(curr_dataset, "genes2tad", "all_assigned_regions.txt") 
stopifnot(file.exists(tadDT_file))
tadDT <- read.delim(tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c( "chromo", "region", "start", "end"))
stopifnot(is.numeric(tadDT$start))
stopifnot(is.numeric(tadDT$end))

# assigned genes to region
i=1
tad_assignment <- foreach(i = 1:nrow(mrnaDT), .combine='c') %dopar% {
  curr_chromo <- mrnaDT$chromo[i]
  curr_start <- mrnaDT$start[i]
  curr_end <- mrnaDT$end[i]
  if(! curr_chromo %in% tadDT$chromo) return(NA)
  matchTAD <- tadDT$region[which(tadDT$chromo == curr_chromo & tadDT$start <= curr_start & tadDT$end >= curr_start)]
  stopifnot(length(matchTAD) == 1)
  
  if(! grepl("_TAD", matchTAD)) {
    matchTAD2 <- tadDT$region[which(tadDT$chromo == curr_chromo & tadDT$start <= curr_end & tadDT$end >= curr_end)]
    if(grepl("_TAD", matchTAD2)) matchTAD <- matchTAD2
  }
  stopifnot(length(matchTAD) == 1)
  matchTAD
}
stopifnot(length(tad_assignment) == nrow(mrnaDT))

mrnaDT$region <- tad_assignment

gene2tadDT_file <- file.path(curr_dataset, "genes2tad", "all_genes_positions.txt") 
stopifnot(file.exists(gene2tadDT_file))
g2tDT <- read.delim(gene2tadDT_file, stringsAsFactors = FALSE, header=F, col.names=c("entrezID", "chromo", "start", "end", "region"))
g2tDT$entrezID <- as.character(g2tDT$entrezID)

dataset_geneListFile <- file.path(pipOutFolder, 
                                  curr_dataset, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")
stopifnot(file.exists(dataset_geneListFile))
pipeline_geneList <- eval(parse(text = load(dataset_geneListFile)))

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
entrezDT <- read.delim(entrezDT_file, header=T, stringsAsFactors = FALSE)
entrezDT$entrezID <- as.character(entrezDT$entrezID)


cond1_columns <- colnames(mrnaDT)[grepl(paste0("_", cond1), colnames(mrnaDT))]
cond2_columns <- colnames(mrnaDT)[grepl(paste0("_", cond2), colnames(mrnaDT))]

if(exprds == "TCGAcoad_msi_mss") {
  stopifnot(length(cond1_columns) == check_nCond1)
  stopifnot(length(cond2_columns) == check_nCond2)
} else{
  stop("not available\n")
}

top_tads <- unique(as.character(topTADdt$region))

mrnaDT <- mrnaDT[!is.na(mrnaDT$region),]

# my_comparisons <- list( c("MSI", "MSS"))
# p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
#   stat_compare_means(label.y = 2)                   # Add global p-value
tad=top_tads[1]
for(tad in top_tads) {
  
  cat(paste0("*** ", tad, "*** \n"), file = logFile, append=T)
  cat(paste0("*** ", tad, "*** \n"))
  
  curr_df <- mrnaDT[mrnaDT$region == tad, c("entrezID", cond1_columns, cond2_columns)]
  
  if(nrow(curr_df) == 0) {
    cat(paste0("! No probe mapping to topTAD: ", tad, "!\n"), append =T , file = logFile)
    cat("! No probe mapping to topTAD: ", tad, "!\n")
    next
  }
  
  g2t_genes <- g2tDT$entrezID[g2tDT$region == tad]
  dataset_genes <- pipeline_geneList[pipeline_geneList %in% g2t_genes]
  
  
  # take only genes present in TCGA dataset
  cat(paste0("... # genes before subsetting only reference genes:\t", nrow(curr_df), "\n"), file = logFile, append=T)
  cat(paste0("... # genes before subsetting only reference genes:\t", nrow(curr_df), "\n"))
  curr_df <- curr_df[curr_df$entrezID %in% dataset_genes,]
  cat(paste0("... # genes before subsetting after reference genes:\t", nrow(curr_df), "\n"), file = logFile, append=T)
  cat(paste0("... # genes before subsetting after reference genes:\t", nrow(curr_df), "\n"))
  
  if(nrow(curr_df) == 0) {
    cat(paste0("! No probe mapping to topTAD: ", tad, "!\n"), append =T , file = logFile)
    cat("! No probe mapping to topTAD: ", tad, "!\n")
    next
  }
  
  varying_cols <- c(which(colnames(curr_df) %in% cond1_columns), which(colnames(curr_df) %in% cond2_columns))
  curr_df_m <- reshape(curr_df, idvar="entrezID", direction = "long", varying=varying_cols,  timevar="sample", sep="")
  curr_df_m$status <- gsub(".+_(.+)", "\\1", curr_df_m$sample)
  colnames(curr_df_m)[3] <- "mRNA"
  curr_df_m$mRNA_log10 <- log10(curr_df_m$mRNA+0.001)
  
  write.table(curr_df_m, row.names = F, col.names = T, append = T, file = logFile, quote=F)
  

  subTit <- paste0("nEntrezID in TAD = ", length(unique(curr_df_m$entrezID)), "/", nrow(mrnaDT), 
                   "\n# ", cond1, " = ", length(cond1_columns), 
                   " - # ", cond2, " = ", length(cond2_columns))
  
  for(countType in plot_countType) {
    
    p <- ggboxplot(curr_df_m, 
                   x = "status",
                   xlab = paste0("CRC status"), 
                   y = paste0("mRNA", countType),
                   ylab = paste0("mRNA", countType),
                   legend.title="",
                   # legend="right",
                   legend="none",
                   title = paste0(tad, ": mRNA (vEntrez)"),
                   subtitle=subTit,
                   fill = "status", palette = c("steelblue3", "tan2"),
                   add = "jitter")
    p <- p + font("legend.text", size = 14)
    p <- p + font("subtitle", face="bold")
    p <- p + font("title", size = 18, face="bold")
    p <- p + font("xlab", size = 16, face="bold")
    p <- p + font("ylab", size = 16)
    p <- p + font("xy.text", size = 14)
    p <- p + theme(plot.title = element_text(hjust=0.5))
    p <- p + stat_compare_means(label.x=0.5)
    if(SSHFS) p
    
    outFile <- file.path(outFold, paste0(curr_dataset, "_", exprds, "_", tad, "_mRNA", countType, "_vEntrez.", plotType))
    ggsave(plot=p, filename = outFile, height=myHeight, width=myWidth )
    cat(paste0("... written: ", outFile, "\n"))
    
  }
  
}



##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

