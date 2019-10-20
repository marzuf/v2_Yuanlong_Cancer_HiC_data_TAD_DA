options(scipen=100)

setDir = ""

buildTable <- FALSE

# Rscript report_figure3.R   # 

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "report_figure3.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(ggplot2)
require(ggpubr)
require(ggsci)
require(ggforce)
require(reshape2)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

all_cols[all_cols == "red"] <- "brown3"
all_cols[all_cols == "blue"] <- "darkblue"
all_cols[all_cols == "green"] <- "forestgreen"

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.2

myWidthGG <- 12
myHeightGG <- 8

pipOutFolder <- file.path( "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

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

outFolder <- "REPORT_FIGURE3"
dir.create(outFolder, recursive=TRUE)

pvalGenes001 <- 0.01
pvalGenes005 <- 0.05

pvalTADs001 <- 0.01
pvalTADs005 <- 0.05


# col1 <- get_palette("Dark2", 3)[1]
# col2 <- get_palette("Dark2", 3)[2]
# col3 <- get_palette("Dark2", 3)[3]
col1 <- pal_d3()(3)[1]
col2 <- pal_d3()(3)[2]
col3 <- pal_d3()(3)[3]


final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

if(buildTable) {
  all_signif_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    g2tFile <- file.path( hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      geneList_file <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      geneList <- get(load(geneList_file))
      
      stopifnot(geneList %in% g2t_DT$entrezID)
      
      hicds_g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
      all_regs <- as.character(hicds_g2t_DT$region)
      
      hicds_final_DT <- final_table_DT[final_table_DT$hicds == hicds & final_table_DT$exprds == exprds,]
      
      stopifnot(setequal(all_regs, hicds_final_DT$region))
      
      topTable_DT_file <- file.path(pipOutFolder,  hicds, exprds, script1_name, "DE_topTable.Rdata")
      stopifnot(file.exists(topTable_DT_file))
      topTable_DT <- get(load(topTable_DT_file))
      topTable_DT$genes <- as.character(topTable_DT$genes)
      stopifnot(names(geneList) %in% topTable_DT$genes)
      
      topTable_DT <- topTable_DT[topTable_DT$genes %in% names(geneList),]
      topTable_DT$entrezID <- geneList[paste0(topTable_DT$genes)]
      stopifnot(!is.na(topTable_DT$entrezID))
      stopifnot(!duplicated(topTable_DT$entrezID))
      
      topTable_DT[, paste0("gene_signifAdjPval_", pvalGenes001)] <- topTable_DT$adj.P.Val <= pvalGenes001
      topTable_DT[, paste0("gene_signifAdjPval_", pvalGenes005)] <- topTable_DT$adj.P.Val <= pvalGenes005
      
      hicds_final_DT[, paste0("signifAdjPvalComb_", pvalTADs001)] <- hicds_final_DT$adjPvalComb <= pvalTADs001
      hicds_final_DT[, paste0("signifAdjPvalComb_", pvalTADs005)] <- hicds_final_DT$adjPvalComb <= pvalTADs005
      hicds_final_DT[, paste0("signifAdjPvalComb_", pvalTADs001, "_signifFDR_0.2")] <- hicds_final_DT[paste0("signifAdjPvalComb_", pvalTADs001)] & hicds_final_DT[paste0("signifFDR_0.2")]
      
      out_dt <- merge(hicds_g2t_DT[, c("entrezID", "region")], hicds_final_DT[, c(grepl("signif|region$", colnames(hicds_final_DT)))], by ="region")
      out_dt <- merge(out_dt, topTable_DT[, c(grepl("entrezID|signif", colnames(topTable_DT)))], by="entrezID")
      out_dt$hicds <- hicds
      out_dt$exprds <- exprds
      out_dt
    }
    exprds_dt
  }
  outFile <- file.path(outFolder, "all_signif_dt.Rdata")
  save(all_signif_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_signif_dt.Rdata")
  all_signif_dt <- get(load(outFile))
}

# nSignifFDR01 <- aggregate(signifFDR_0.1 ~ region, data=all_signif_dt, FUN=sum)
# nSignifFDR02 <- aggregate(signifFDR_0.2 ~ region, data=all_signif_dt, FUN=sum)
# nSignifAdjPvalComb001 <- aggregate(signifAdjPvalComb_0.01 ~ region, data=all_signif_dt, FUN=sum)
# nSignifAdjPvalComb005 <- aggregate(signifAdjPvalComb_0.05 ~ region, data=all_signif_dt, FUN=sum)
# 

signif_gene_vars <- c("gene_signifAdjPval_0.05", "gene_signifAdjPval_0.01")
signif_tad_vars <- colnames(all_signif_dt)[! colnames(all_signif_dt) %in% c("hicds", "exprds", "region", "entrezID", signif_gene_vars)]

signif_gene=signif_gene_vars[1]
signif_tad=signif_tad_vars[1]

signif_gene =  "gene_signifAdjPval_0.05"
signif_tad =  "signifAdjPvalComb_0.01_signifFDR_0.2"

signif_gene_vars =  "gene_signifAdjPval_0.05"
signif_tad_vars =  "signifAdjPvalComb_0.01"

signif_gene = signif_gene_vars
signif_tad =  signif_tad_vars
signif_gene_thresh <- gsub("gene_signifAdjPval_", "", signif_gene)
signif_tad_thresh <- gsub("signifAdjPvalComb_", "", signif_tad)
bothsignifcol <- paste0(signif_tad, "_AND_", signif_gene)


signif_geneName <- paste0("limma p-val. <= ", signif_gene_thresh)
signif_tadName <- paste0("TAD adj. p-val <= ", signif_tad_thresh)
bothsignifcolName <- paste0(signif_geneName, " AND\n", signif_tadName)
newVariableNames <- setNames(c(signif_geneName, signif_tadName, bothsignifcolName),
                               c(signif_gene, signif_tad, bothsignifcol))
                               
                                 
                                 


for(signif_gene in signif_gene_vars) {
  
  for(signif_tad in signif_tad_vars) {
    
    cat(" signif_gene = ", signif_gene, "\n" )
    cat(" signif_tad = ", signif_tad, "\n" )
    
    
    orderCol <- signif_gene
    
    
    all_signif_dt[,paste0(bothsignifcol)] <- all_signif_dt[,paste0(signif_gene)] & all_signif_dt[,paste0(signif_tad)]
    
    
    nSignifGenes_tad <- aggregate( as.formula(paste0(signif_tad, " ~ hicds + exprds")), FUN=sum, data=all_signif_dt)
    colnames(nSignifGenes_tad)[3] <- paste0(signif_tad)
    
    nSignifGenes_limma <- aggregate( as.formula(paste0(signif_gene, " ~ hicds + exprds")), FUN=sum, data=all_signif_dt)
    colnames(nSignifGenes_limma)[3] <- paste0(signif_gene)
    
    
    nSignifGenes_tadLimma <- aggregate( as.formula(paste0(paste0(bothsignifcol), " ~ hicds + exprds")), FUN=sum, data=all_signif_dt)
    colnames(nSignifGenes_tadLimma)[3] <- paste0(bothsignifcol)
    
    
    nSignif_tad_limma <- merge(nSignifGenes_limma, nSignifGenes_tad, by=c("hicds", "exprds"))
    
    nSignif_dt <- merge(nSignif_tad_limma, nSignifGenes_tadLimma, by=c("hicds", "exprds"))
    nSignif_dt$dataset <- paste0(nSignif_dt$hicds, "\n", nSignif_dt$exprds)
    nSignif_dt <- nSignif_dt[order(nSignif_dt[,paste0(orderCol)], decreasing = TRUE),]
    ds_levels <- nSignif_dt$dataset
    
    
    nSignif_dt_m <- melt(nSignif_dt, by=c("hicds", "exprds", "dataset"))
    
    nSignif_dt_m$dataset <- factor(nSignif_dt_m$dataset, levels=ds_levels)
    nSignif_dt_m <- nSignif_dt_m[order(as.numeric(nSignif_dt_m$dataset)),]
    nSignif_dt_m$exprds_type_col <- all_cols[all_cmps[nSignif_dt_m$exprds]]
    mycols <- nSignif_dt_m$exprds_type_col[as.character(nSignif_dt_m$variable) == signif_tad] 
    
    nSignif_dt_m$variable <- as.character(nSignif_dt_m$variable)
    nSignif_dt_m$variable <- newVariableNames[nSignif_dt_m$variable]
    nSignif_dt_m$variable <- factor(nSignif_dt_m$variable, levels=as.character(newVariableNames))
    stopifnot(!is.na(nSignif_dt_m$variable))
    
    # nSignif_dt_m$variable <- gsub("_AND_", "_AND\n", nSignif_dt_m$variable)
    
    
    p_var <-  ggplot(nSignif_dt_m, aes(x = dataset, y = value, fill = variable)) + 
      geom_bar(position="dodge", stat="identity") +
      coord_cartesian(expand = FALSE) +
      ggtitle("# signif. genes", subtitle = paste0(signif_tad, " + ", signif_gene))+
      scale_x_discrete(name="")+
      labs(fill="")+
      scale_fill_manual(values=c(col1,col2, col3))+
      scale_y_continuous(name=paste0("# signif. genes"),
                         breaks = scales::pretty_breaks(n = 10))+
      theme( # Increase size of axis lines
        strip.text = element_text(size = 12),
        # top, right, bottom and left
        # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size=16),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.minor.y = element_line(colour = "grey"),
        strip.text.x = element_text(size = 10),
        axis.line.x = element_line(size = .2, color = "black"),
        axis.line.y = element_line(size = .3, color = "black"),
        axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
        axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank(),
        legend.title = element_text(face="bold")
      )
    
    outFile <- file.path(outFolder, paste0("nSignif_genes_", signif_tad, "_and_", signif_gene, "_geneSignifOrder.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    ### zoom
    
    p_var_zoom <- p_var + facet_zoom(y=variable!= paste0(signif_gene), zoom.data = ifelse(variable==paste0(signif_gene), FALSE, NA))
    outFile <- file.path(outFolder, paste0("nSignif_genes_", signif_tad, "_and_", signif_gene, "_zoom_geneSignifOrder.", plotType))
    ggsave(plot = p_var_zoom, filename = outFile, height=myHeightGG, width = myWidthGG*2)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    ### SAME BUT ORDERING ACCORDING TO TAD
    
    orderCol <- signif_tad
    nSignif_dt <- nSignif_dt[order(nSignif_dt[,paste0(orderCol)], decreasing = TRUE),]
    ds_levels <- nSignif_dt$dataset
    
    
    nSignif_dt_m <- melt(nSignif_dt, by=c("hicds", "exprds", "dataset"))
    
    nSignif_dt_m$dataset <- factor(nSignif_dt_m$dataset, levels=ds_levels)
    nSignif_dt_m <- nSignif_dt_m[order(as.numeric(nSignif_dt_m$dataset)),]
    nSignif_dt_m$exprds_type_col <- all_cols[all_cmps[nSignif_dt_m$exprds]]
    mycols <- nSignif_dt_m$exprds_type_col[as.character(nSignif_dt_m$variable) == signif_tad] 
    
    
    nSignif_dt_m$variable <- as.character(nSignif_dt_m$variable)
    nSignif_dt_m$variable <- newVariableNames[nSignif_dt_m$variable]
    nSignif_dt_m$variable <- factor(nSignif_dt_m$variable, levels=as.character(newVariableNames))
    stopifnot(!is.na(nSignif_dt_m$variable))
    
    # nSignif_dt_m$variable <- gsub("_AND_", "_AND\n", nSignif_dt_m$variable)
    
    
    p_var <-  ggplot(nSignif_dt_m, aes(x = dataset, y = value, fill = variable)) + 
      geom_bar(position="dodge", stat="identity") +
      coord_cartesian(expand = FALSE) +
      ggtitle("# signif. genes", subtitle = paste0(signif_tad, " + ", signif_gene))+
      scale_x_discrete(name="")+
      labs(fill="")+
      scale_fill_manual(values=c(col1,col2, col3))+
      scale_y_continuous(name=paste0("# signif. genes"),
                         breaks = scales::pretty_breaks(n = 10))+
      theme( # Increase size of axis lines
        strip.text = element_text(size = 12),
        # top, right, bottom and left
        # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size=16),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.minor.y = element_line(colour = "grey"),
        strip.text.x = element_text(size = 10),
        axis.line.x = element_line(size = .2, color = "black"),
        axis.line.y = element_line(size = .3, color = "black"),
        axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
        axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank(),
        legend.title = element_text(face="bold")
      )
    
    outFile <- file.path(outFolder, paste0("nSignif_genes_", signif_tad, "_and_", signif_gene, "_tadSignifOrder.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    ### zoom
    
    p_var_zoom <- p_var + facet_zoom(y=variable!= paste0(signif_gene), zoom.data = ifelse(variable==paste0(signif_gene), FALSE, NA))
    outFile <- file.path(outFolder, paste0("nSignif_genes_", signif_tad, "_and_", signif_gene, "_zoom_tadSignifOrder.", plotType))
    ggsave(plot = p_var_zoom, filename = outFile, height=myHeightGG, width = myWidthGG*2)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
    
    
    nSignif_dt_m$value_log10 <- log10(nSignif_dt_m$value)
    p_var <-  ggplot(nSignif_dt_m, aes(x = dataset, y = value_log10, fill = variable)) + 
      geom_bar(position="dodge", stat="identity") +
      coord_cartesian(expand = FALSE) +
      ggtitle("# signif. genes", subtitle = paste0(signif_tad, " + ", signif_gene))+
      scale_x_discrete(name="")+
      labs(fill="")+
      scale_fill_manual(values=c(col1,col2, col3))+
      scale_y_continuous(name=paste0("# signif. genes (log10)"),
                         breaks = scales::pretty_breaks(n = 10))+
      theme( # Increase size of axis lines
        strip.text = element_text(size = 12),
        # top, right, bottom and left
        # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size=16),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.minor.y = element_line(colour = "grey"),
        strip.text.x = element_text(size = 10),
        axis.line.x = element_line(size = .2, color = "black"),
        axis.line.y = element_line(size = .3, color = "black"),
        axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
        axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color="black", size=12),
        axis.title.x = element_text(color="black", size=12),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        legend.background =  element_rect(),
        legend.key = element_blank(),
        legend.title = element_text(face="bold")
      )
    
    outFile <- file.path(outFolder, paste0("nSignif_genes_", signif_tad, "_and_", signif_gene, "_log10.", plotType))
    ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
    
    
    
    
    
  }
}








##########################################################################################
##########################################################################################
##########################################################################################

cat("*** DONE - ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))