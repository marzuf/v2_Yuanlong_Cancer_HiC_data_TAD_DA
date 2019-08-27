options(scipen=100)

setDir = ""

buildTable <- TRUE

options(save.defaults = list(version = 2))

# Rscript geneHancer_TADs.R   # 

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "geneHancer_TADs.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

require(stringr)
require(ggplot2)
require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.2

myWidthGG <- 12
myHeightGG <- 8

col1 <- "dodgerblue3"
col2 <- "darkorange3"

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

outFolder <- "GENEHANCER_TADS"
dir.create(outFolder, recursive=TRUE)

logFile <- file.path(outFolder, "genehancer_tads_logfile.txt")
if(buildTable) file.remove(logFile)

# enhancer data download from https://www.genecards.org/GeneHancer_version_4-4
enhancerFile <- "GeneHancer_version_4-4.csv"
stopifnot(file.exists(enhancerFile))

enhancerDT <- read.delim(enhancerFile, stringsAsFactors = FALSE, sep="\t", header=TRUE)
stopifnot(nrow(enhancerDT) > 0)

stopifnot(enhancerDT$feature.name == "Enhancer")

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)


final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

if(buildTable) {
  
  gene2enhancerDT <- foreach(i = 1:nrow(enhancerDT), .combine='rbind') %dopar% {
    cat(paste0("...... gene2enhancer: ", i, "/", nrow(enhancerDT),"\n"))
    i_att <- enhancerDT$attributes[i]
    enhancer_id <- gsub("genehancer_id=(.+);", "\\1", unlist(str_extract_all(string=i_att, pattern = "genehancer_id=(.+?);" )))
    stopifnot(length(enhancer_id) == 1)
    all_genes <- gsub("connected_gene=(.+);", "\\1", unlist(str_extract_all(string=i_att, pattern = "connected_gene=(.+?);" )))
    data.frame(
      enhancer_id = enhancer_id,
      connected_gene = all_genes,
      stringsAsFactors = FALSE
    )
  }
  outFile <- file.path(outFolder, "gene2enhancerDT.Rdata")
  save(gene2enhancerDT, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))


  
  all_enhancer2tad_dt <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    cat(paste0("... start enhancer2tad for: ", hicds,"\n"))
    tadFile <- file.path( hicds, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(tadFile))
    tad_DT <- read.delim(tadFile, header=F, col.names = c("chromo","region", "start", "end"), stringsAsFactors = FALSE)
    tad_DT <- tad_DT[grepl("_TAD", tad_DT$region),]
    enhancer2tad_dt <- foreach(i=1:nrow(enhancerDT), .combine='rbind') %dopar% {
      cat(paste0("...... enhancer2tad: ", i, "/", nrow(enhancerDT),"\n"))
      enhancer_id <- gsub("genehancer_id=(.+);", "\\1", unlist(str_extract_all(string=enhancerDT$attributes[i], pattern = "genehancer_id=(.+?);" )))
      idx_match <- which(tad_DT$chromo == enhancerDT$chrom[i] & tad_DT$start <= enhancerDT$start[i] & tad_DT$end >= enhancerDT$end[i])
      if(length(idx_match) == 0) return(NULL)
      tad_match <- tad_DT$region[idx_match]
      stopifnot( length(tad_match) == 1 )
      
      data.frame(
        hicds=hicds,
        enhancer_id = enhancer_id,
        enhancer_start = enhancerDT$start[i],
        enhancer_end = enhancerDT$end[i],
        enhancer_tad = tad_match,
        tad_start = tad_DT$start[idx_match],
        tad_end = tad_DT$end[idx_match],
        stringsAsFactors = FALSE
      )
    } # end-foreach iterating over enhancers
    txt <- paste0("> # of enhancers within TADs - ", hicds, "\t=\t", nrow(enhancer2tad_dt),"/", nrow(enhancerDT), "\n")
    cat(txt)
    cat(txt, append=TRUE, file = logFile)
    enhancer2tad_dt
  } # end-foreach iterating over hicds
  outFile <- file.path(outFolder, "all_enhancer2tad_dt.Rdata")
  save(all_enhancer2tad_dt, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
 
  
  
  txt <- paste0("> # of enhancers\t=\t", nrow(enhancerDT), "\n")
  cat(txt)
  cat(txt, append=TRUE, file = logFile)
  
  txt <- paste0("> # of genes connected to enhancer\t=\t", nrow(gene2enhancerDT), "\n")
  cat(txt)
  cat(txt, append=TRUE, file = logFile)
  
  entrez_gene2enhancerDT <- gene2enhancerDT[gene2enhancerDT$connected_gene %in% entrez2symb_dt$symbol,]
  txt <- paste0("> # of genes connected to enhancer with EntrezID\t=\t", nrow(entrez_gene2enhancerDT), "\n")
  cat(txt)
  cat(txt, append=TRUE, file = logFile)
  
  
  
} else { # end-if buildTable
  
  outFile <- file.path(outFolder, "all_enhancer2tad_dt.Rdata")
  all_enhancer2tad_dt <- get(load(outFile))
  
  outFile <- file.path(outFolder, "gene2enhancerDT.Rdata")
  gene2enhancerDT <- get(load(outFile))
  
}


# compute number of enhancer per TAD

nEnhancers2tad_dt <- aggregate(enhancer_id ~ hicds + enhancer_tad, data = all_enhancer2tad_dt, FUN=length)
colnames(nEnhancers2tad_dt)[colnames(nEnhancers2tad_dt) == "enhancer_id"] <- "nEnhancers"


all_signif_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {

  curr_e2t_dt <- nEnhancers2tad_dt[nEnhancers2tad_dt$hicds == hicds,]
  stopifnot(nrow(curr_e2t_dt) > 0) 
  stopifnot(!duplicated(curr_e2t_dt$enhancer_tad))
  all_nEnhancers <- setNames(curr_e2t_dt$nEnhancers, curr_e2t_dt$enhancer_tad)
  
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    
    comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
    stopifnot(file.exists(comb_empPval_file))
    comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
    
    # ADJUST THE PVAL
    adj_empPval_comb <- p.adjust(comb_empPval, method="BH")
    
    all_regs <- names(adj_empPval_comb)
    
    nEnhancers <- all_nEnhancers[paste0(all_regs)]
    stopifnot(length(nEnhancers) == length(all_regs))
    nEnhancers[is.na(nEnhancers)] <- 0
    
    
    data.frame(
      hicds=hicds,
      exprds=exprds,
      region = all_regs,
      adj_empPval_comb=adj_empPval_comb[all_regs], 
      nEnhancers = nEnhancers,
      stringsAsFactors = FALSE
      )
    

  } # end-foreach iterating over exprds
  exprds_dt
} # end-foreach iterating over hicds

########################################################################################## densplot

x_var <- "adj_empPval_comb"
y_var <- "nEnhancers"
myx <- all_signif_dt[,paste0(x_var)]
myy <- all_signif_dt[,paste0(y_var)]

nDS <- length(unique(paste0(all_signif_dt$hicds, all_signif_dt$exprds)))

outFile <- file.path(outFolder, paste0(y_var, "_log10_vs_", x_var, "_log10_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x=-log10(myx),
  y=log10(myy),
  xlab = paste0(x_var, " [-log10]"),
  ylab = paste0(y_var, " [log10]"),
  main = paste0(y_var, " (log10) vs. ", x_var, "(log10)"),
  cex.lab=plotCex,
  cex.axis = plotCex
)
addCorr(x = myx, y = myy, bty="n")
mtext(side=3, text = paste0("nDS = ", nDS))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_log10_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x=-log10(myx),
  y=(myy),
  xlab = paste0(x_var, " [-log10]"),
  ylab = paste0(y_var, ""),
  main = paste0(y_var, " vs. ", x_var, "(log10)"),
  cex.lab=plotCex,
  cex.axis = plotCex
)
addCorr(x = myx, y = myy, bty="n")
mtext(side=3, text = paste0("nDS = ", nDS))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0(y_var, "_vs_", x_var, "_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
densplot(
  x=(myx),
  y=(myy),
  xlab = paste0(x_var, ""),
  ylab = paste0(y_var, ""),
  main = paste0(y_var, " vs. ", x_var, ""),
  cex.lab=plotCex,
  cex.axis = plotCex
)
addCorr(x = myx, y = myy, bty="n")
mtext(side=3, text = paste0("nDS = ", nDS))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


########################################################################################## boxplot

final_table_DT[, "signifAdjPvalComb_0.05"] <- final_table_DT$adjPvalComb <= 0.05
final_table_DT[, "signifAdjPvalComb_0.01"] <- final_table_DT$adjPvalComb <= 0.01

colnames(nEnhancers2tad_dt)[colnames(nEnhancers2tad_dt) == "enhancer_tad"] <- "region"

boxplot_dt <- merge(final_table_DT, nEnhancers2tad_dt, by=c("hicds", "region"))

boxplot_dt$dataset <- paste0(boxplot_dt$hicds, "\n", boxplot_dt$exprds)

meanDT <- aggregate(nEnhancers ~ dataset, FUN=mean,data = boxplot_dt)
meanDT <- meanDT[order(meanDT$nEnhancers, decreasing = TRUE),]
ds_levels <- meanDT$dataset

mycols <- all_cols[all_cmps[paste0(gsub(".+\n(.+)", "\\1", ds_levels))]]
stopifnot(!is.na(mycols))

boxplot_dt$dataset <- factor(boxplot_dt$dataset, levels = ds_levels)

signif_vars <- colnames(boxplot_dt)[grepl("signif", colnames(boxplot_dt))]

signif_var = signif_vars[1]

boxplot_dt$nEnhancers_log10 <- log10(boxplot_dt$nEnhancers)

for(signif_var in signif_vars){
  
  
  p_var <-  ggplot(boxplot_dt, aes_string(x = "dataset", y = paste0("nEnhancers"), fill = paste0(signif_var))) + 
    geom_boxplot() +
    # coord_cartesian(expand = FALSE) +
    ggtitle("# enhancers in TADs", subtitle = paste0(signif_var))+
    scale_x_discrete(name="")+
    # labs(fill="")+
    scale_fill_manual(values=c(col1,col2))+
    scale_y_continuous(name=paste0("# enhancers by TAD"),
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
  outFile <- file.path(outFolder, paste0("all_ds_", signif_var, ".", plotType))
  ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  p_var <-  ggplot(boxplot_dt, aes_string(x = "dataset", y = paste0("nEnhancers_log10"), fill = paste0(signif_var))) + 
    geom_boxplot() +
    # coord_cartesian(expand = FALSE) +
    ggtitle("# enhancers in TADs (log10)", subtitle = paste0(signif_var))+
    scale_x_discrete(name="")+
    # labs(fill="")+
    scale_fill_manual(values=c(col1,col2))+
    scale_y_continuous(name=paste0("# enhancers by TAD [log10]"),
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
  outFile <- file.path(outFolder, paste0("all_ds_", signif_var, "_log10.", plotType))
  ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
}




##########################################################################################
##########################################################################################
##########################################################################################

cat("*** DONE - ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))