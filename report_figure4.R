options(scipen=100)

SSHFS=F

setDir <- "/media/electron"
setDir <- ""

# Rscript report_figure4.R

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "report_figure4.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(RColorBrewer)
require(reshape2)
require(ggplot2)
require(ggpubr)
# require(gplots)
registerDoMC(4)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb <- setNames(entrez2symb_dt$symbol, entrez2symb_dt$entrezID)


buildTable <- TRUE

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.2

myHeightHeat <- myHeight * 1.8
myWidthHeat <- myWidth * 1.8

myWidthGG <- 12
myHeightGG <- 7

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script11same_name <- "11sameNbr_runEmpPvalCombined"

tieMeth <- "min"


mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

outFolder <- file.path("REPORT_FIGURE4")
dir.create(outFolder, recursive = TRUE)


final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

final_table_DT$cmpType <- all_cmps[paste0(final_table_DT$exprds)]
stopifnot(!is.na(final_table_DT$cmpType))



tad_signif_col <- "tad_adjCombPval"
gene_signif_col <- "adj.P.Val"

signif_column <- "adjPvalComb"
signifThresh <- 0.01




tad_pval_thresh <- 0.01
gene_pval_thresh <- 0.05
tad_pval_thresh_log10 <- -log10(tad_pval_thresh)
gene_pval_thresh_log10 <- -log10(gene_pval_thresh)
file_suffix <- paste0("tad_pval_thresh", tad_pval_thresh, "_gene_pval_thresh", gene_pval_thresh)



minOverlapBpRatio <- 0.8
minIntersectGenes <- 3

nRegionLolli <- 10

cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> signifThresh\t=\t", signifThresh, "\n"))
cat(paste0("> minOverlapBpRatio\t=\t", minOverlapBpRatio, "\n"))
cat(paste0("> minIntersectGenes\t=\t", minIntersectGenes, "\n"))
cat(paste0("> nRegionLolli\t=\t", nRegionLolli, "\n"))


tad_plot_list <- list()

data_cmpType = "norm_vs_tumor"

for(data_cmpType in c("norm_vs_tumor", "subtypes", "wt_vs_mut", "")) {
  
  if(data_cmpType == "") {
    nDS <- length(unique(paste0(final_table_DT$hicds,final_table_DT$exprds)))  
  } else {
    nDS <- length(unique(paste0(final_table_DT$hicds[final_table_DT$cmpType==data_cmpType],final_table_DT$exprds[final_table_DT$cmpType==data_cmpType])))
  }
  stopifnot(nDS > 0)
  
  
  inFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", data_cmpType)
  stopifnot(dir.exists(inFolder))
  
  inFile <- file.path(inFolder, paste0(ifelse(data_cmpType=="", data_cmpType, paste0(data_cmpType, "_")), 
                      "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio,
                      "_minInterGenes", minIntersectGenes, ".Rdata"))
  stopifnot(file.exists(inFile))
  conserved_signif_tads <- get(load(inFile))
  
  countConserv <- abs(sort(-lengths(conserved_signif_tads)))
  stopifnot(countConserv > 1)
  
  countConserv_dt <- data.frame(countConserv)
  countConserv_dt$region <- factor(rownames(countConserv_dt), levels=rownames(countConserv_dt))
  countConserv_dt$conservRatio <- countConserv_dt$countConserv/nDS
  stopifnot(countConserv_dt$conservRatio <= 1)
  
  
  
  
  tad_plot_list[[paste0(ifelse(data_cmpType=="", "all", data_cmpType))]] <- ggbarplot(countConserv_dt, 
                                                                                      title=paste0("Conserv. signif. regions"),
                                                                                      subtitle=paste0("# DS = ", nDS, "; n = ", nrow(countConserv_dt)),
                                                                                      x = "region", 
                                                                                      y = "conservRatio",
                                                                                      xlab="",
                                                                                      ylab="conserv. ratio",
                                                                                      col="darkgrey", fill="darkgrey") +
                                                            theme(axis.text.x=element_blank(), 
                                                                  axis.ticks.x = element_blank()) + 
                                                            coord_cartesian(expand = F)
  
  
}




####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start preparing data before matching \n")
  
  hicds = all_hicds[1]
  all_gene_tad_signif_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      
      regionList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
      stopifnot(file.exists(regionList_file))
      all_regs <- get(load(regionList_file))
      
      
      geneList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      pipeline_geneList <- get(load(geneList_file))
      
      stopifnot(pipeline_geneList %in% g2t_dt$entrezID)
      
      exprds_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% pipeline_geneList,]
      
      
      comb_empPval_file <- file.path(pipFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      stopifnot(setequal(all_regs, names(comb_empPval)))
      comb_empPval <- comb_empPval[all_regs]
      # ADJUST THE PVAL
      tad_adjCombPval <- p.adjust(comb_empPval, method="BH")
      stopifnot(names(tad_adjCombPval) == all_regs)
      
      
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
      stopifnot(exprds_g2t_dt$entrezID %in% topTable_DT$entrezID)
      
      
      tad_dt <- data.frame(region=names(tad_adjCombPval), tad_adjCombPval = as.numeric(tad_adjCombPval), stringsAsFactors=FALSE)
      
      tad_gene_dt <- merge(exprds_g2t_dt[,c("entrezID", "region")], tad_dt, by="region", all.x=TRUE, all.y=FALSE)
      
      out_dt <- merge(tad_gene_dt, topTable_DT[,c("entrezID", "logFC", "adj.P.Val")], by="entrezID", all.x=TRUE, all.y=FALSE )
      out_dt <- unique(out_dt)
      stopifnot(!duplicated(out_dt$entrezID))
      stopifnot(exprds_g2t_dt$entrezID %in% out_dt$entrezID)
      
      out_dt$gene_rank <- rank(out_dt$adj.P.Val, ties=tieMeth)
      out_dt$tad_rank <- rank(out_dt$tad_adjCombPval, ties=tieMeth)
      
      out_dt_cols <- colnames(out_dt)
      out_dt$hicds <- hicds
      out_dt$exprds <- exprds
      
      out_dt[, c("hicds", "exprds", out_dt_cols)]
      
    } # end-foreach iterating over exprds
    exprds_dt
  } # end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, paste0("all_gene_tad_signif_dt.Rdata"))
  save(all_gene_tad_signif_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  # outFile <- "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"
  outFile <- file.path(outFolder, paste0( "all_gene_tad_signif_dt.Rdata"))
  cat("... load data\n")
  all_gene_tad_signif_dt <- get(load(outFile))
}


all_gene_tad_signif_dt$gene_signif <- all_gene_tad_signif_dt[,paste0(gene_signif_col)] <= gene_pval_thresh
all_gene_tad_signif_dt$cmpType <- all_cmps[all_gene_tad_signif_dt$exprds]  
stopifnot(!is.na(all_gene_tad_signif_dt$cmpType))

gene_plot_list <- list()

data_cmpType = "norm_vs_tumor"

all_signif_dt <- all_gene_tad_signif_dt[all_gene_tad_signif_dt$gene_signif,]

for(data_cmpType in c("norm_vs_tumor", "subtypes", "wt_vs_mut", "")) {
  

  if(data_cmpType != "") {
    signif_dt <- all_signif_dt[all_signif_dt$cmpType == data_cmpType,]
    nDS <- length(unique(paste0(all_gene_tad_signif_dt$hicds[all_gene_tad_signif_dt$cmpType==data_cmpType],
                                all_gene_tad_signif_dt$exprds[all_gene_tad_signif_dt$cmpType==data_cmpType])))
    
    
    
  } else {
    signif_dt <- all_signif_dt  
    nDS <- length(unique(paste0(all_gene_tad_signif_dt$hicds,all_gene_tad_signif_dt$exprds)))  
  }
  stopifnot(nDS > 0)
  
  conserved_signif_genes <- setNames(as.numeric(table(signif_dt$entrezID)),
                                     as.character(names(table(signif_dt$entrezID))))
  
  conserved_signif_genes <- conserved_signif_genes[conserved_signif_genes > 1]
  
  countConserv <- abs(sort(-(conserved_signif_genes)))
  
  countConserv_dt <- data.frame(countConserv)
  countConserv_dt$gene <- factor(rownames(countConserv_dt), levels=rownames(countConserv_dt))
  countConserv_dt$conservRatio <- countConserv_dt$countConserv/nDS
  stopifnot(countConserv_dt$conservRatio <= 1)

  
  
  gene_plot_list[[paste0(ifelse(data_cmpType=="", "all", data_cmpType))]] <-  ggbarplot(countConserv_dt,  title=paste0("Conserv. signif. genes - ", data_cmpType),
                                                                                        subtitle=paste0("# DS = ", nDS, "; n = ", nrow(countConserv_dt)),
                                                                                        x = "gene", 
                                                                                        y = "conservRatio",
                                                                                        xlab="",
                                                                                        ylab="conserv. ratio",
                                                                                        col="darkgrey", fill="darkgrey") +
    theme(axis.text.x=element_blank(), 
          axis.ticks.x = element_blank()) + 
    coord_cartesian(expand = F)
  
}

outFile <- file.path(outFolder, paste0( "gene_plot_list.Rdata"))
save(gene_plot_list, file=outFile, version=2)
cat("... load data\n")
gene_plot_list <- get(load(outFile))

outFile <- file.path(outFolder, paste0( "tad_plot_list.Rdata"))
cat("... load data\n")
save(tad_plot_list, file=outFile, version=2)
tad_plot_list <- get(load(outFile))



arrAll <- ggarrange(
  ggarrange(
    gene_plot_list[["all"]],
    tad_plot_list[["all"]],
    nrow=2
  ),
  ggarrange(
    gene_plot_list[["norm_vs_tumor"]],
    tad_plot_list[["norm_vs_tumor"]],
    nrow=2
  ),
  ggarrange(
    gene_plot_list[["subtypes"]],
    tad_plot_list[["subtypes"]],
    nrow=2
  ),
  ggarrange(
    gene_plot_list[["wt_vs_mut"]],
    tad_plot_list[["wt_vs_mut"]],
    nrow=2
  ),
  ncol = 4 
)


outFile <- file.path(outFolder, paste0("conserved_genes_regions_all_barplot.", plotType))
ggsave(plot = arrAll, filename = outFile, height=myHeightGG*2, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))

  
######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


