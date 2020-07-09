
# Rscript pands_corrPval_selectedTADs.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390

outFolder <- "PANDS_CORRPVAL_SELECTEDTADS"
dir.create(outFolder, recursive = TRUE)

hicds_oi <- "ENCSR489OCU_NCI-H460_40kb"
exprds_oi <- "TCGAluad_norm_luad"
tad_oi <- "chr11_TAD390"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
hicds_oi <- args[1]
exprds_oi <- args[2]
tad_oi <- args[3]

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script4_name <- "4_runMeanTADCorr"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10sameNbr_runEmpPvalMeanTADCorr"

setDir <- "/media/electron"
setDir <- ""
buildTable <- TRUE

plotType <- "svg"
myWidth <- myHeight <- 6

require(doMC)
require(foreach)
registerDoMC(40)
require(ggplot2)

matchingNbrGenesThresh <- 0.5 # will look at the corr pval of the TADs if at least > matchingNbrGenesThresh of its genes belong to "genes_oi"

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

### retrieve the genes that belong to the TAD of interest
pipeline_geneList <- get(load(file.path(pipFolder, hicds_oi, exprds_oi, script0_name, "pipeline_geneList.Rdata")))
gene2tad_dt <- read.delim(file.path(hicds_oi, "genes2tad", "all_genes_positions.txt"),
                          stringsAsFactors = FALSE,
                          header=FALSE, 
                          col.names=c("entrezID", "chromosome", "start", "end", "region"))
gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)
stopifnot(pipeline_geneList %in% gene2tad_dt$entrezID)
stopifnot(tad_oi %in% gene2tad_dt$region)
g2t_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% pipeline_geneList,]
genes_oi <- g2t_dt$entrezID[g2t_dt$region == tad_oi]

rm("g2t_dt")
rm("gene2tad_dt")
rm("pipeline_geneList")


### retrieve pairwise correlations for all datasets (only if they are in the same TAD)

if(buildTable) {
  hicds = all_hicds[1]
  all_ds_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat(paste0("... start: ", hicds, " - ", exprds, "\n"))
      
      
      meanCorr <- get(load(file.path(pipFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")))
      
      emp_fcPval <- get(load(file.path(pipFolder, hicds, exprds, script9_name, "emp_pval_meanLogFC.Rdata")))
      adj_fcPval <- p.adjust(emp_fcPval, method="BH")
      
      emp_corrPval <- get(load(file.path(pipFolder, hicds, exprds, script10_name, "emp_pval_meanCorr.Rdata")))
      adj_corrPval <- p.adjust(emp_corrPval, method="BH")
      
      # see if there are matching TADs
      
      # source(file.path("PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R")))
      # samp1 <- get(load(file.path(setDir, sample1_file)))
      # samp2 <- get(load(file.path(setDir, sample2_file)))
      # qq_dt <- eval(parse(text = load(file.path(pipFolder, hicds, exprds, script0_name, "rna_qqnorm_rnaseqDT.Rdata")))) 
      # stopifnot(setequal(colnames(qq_dt), c(samp1,samp2)))
      pipeline_geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
      gene2tad_dt <- read.delim(file.path(hicds, "genes2tad", "all_genes_positions.txt"),
                                stringsAsFactors = FALSE,
                                header=FALSE, 
                                col.names=c("entrezID", "chromosome", "start", "end", "region"))
      gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)
      stopifnot(pipeline_geneList %in% gene2tad_dt$entrezID)
      g2t_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% pipeline_geneList,]
      
      tad_size <- setNames(as.numeric(table(g2t_dt$region)), names(table(g2t_dt$region)))
      
      oi_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% genes_oi,]
      
      if(nrow(oi_g2t_dt) > 0) {
        
        match_tad_nbrOi <- setNames(as.numeric(table(oi_g2t_dt$region)), names(table(oi_g2t_dt$region)))  
        match_tad_size <- tad_size[names(match_tad_nbrOi)]
        stopifnot(!is.na(match_tad_size))
        match_tad_oiRatio <- match_tad_nbrOi/match_tad_size
        stopifnot(match_tad_oiRatio > 0 & match_tad_oiRatio <= 1)
        
        
        if(hicds == hicds_oi & exprds == exprds_oi) {
          stopifnot(length(match_tad_nbrOi) == 1)
          stopifnot( match_tad_oiRatio == 1)
          stopifnot(names(match_tad_nbrOi) == tad_oi)
        }
        
        matchingTADs <- names(match_tad_oiRatio)[match_tad_oiRatio > matchingNbrGenesThresh]
        
        
        if(length(matchingTADs) > 0) {
          
          stopifnot(matchingTADs %in% names(adj_corrPval))
          stopifnot(matchingTADs %in% names(meanCorr))
          
          
          return(data.frame(
            hicds = hicds,
            exprds = exprds,
            matchingTADs = matchingTADs,
            ratioGenesOi = match_tad_oiRatio[matchingTADs],
            meanCorr = meanCorr[matchingTADs],
            adj_meanCorr_empPval = adj_corrPval[matchingTADs],
            adj_meanFC_empPval = adj_fcPval[matchingTADs],
            stringsAsFactors = FALSE
          ))
        } else {
          return(NULL)
        }
        
      } else {
       return(NULL) 
      }
    }
    exprds_dt
  }
  outfile <-file.path(outFolder, paste0(hicds_oi, "_", exprds_oi, "_", tad_oi, "_all_ds_dt.Rdata"))
  save(all_ds_dt, file=outfile, version=2)
  cat(paste0("... written:" , outfile, "\n"))
  
} else {
  outfile <-file.path(outFolder, paste0(hicds_oi, "_", exprds_oi, "_", tad_oi, "_all_ds_dt.Rdata"))
  all_ds_dt <- get(load(outfile))
}

all_ds_dt$adj_meanCorr_empPval_log10 <- -log10(all_ds_dt$adj_meanCorr_empPval)
all_ds_dt$adj_meanFC_empPval_log10 <- -log10(all_ds_dt$adj_meanFC_empPval)


require(ggpubr)

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

all_ds_dt$cmpType <- all_cmps[all_ds_dt$exprds]
stopifnot(!is.na(all_ds_dt$cmpType))
all_ds_dt$cmpCols <- all_cols[all_ds_dt$cmpType]
stopifnot(!is.na(all_ds_dt$cmpCols))

all_ds_dt$oiType <- ifelse(all_ds_dt$hicds == hicds_oi & all_ds_dt$exprds == exprds_oi, "oi", "other")
oiTypes <- c(oi="red", other="black")
all_ds_dt$oiCols <- oiTypes[paste0(all_ds_dt$oiType)]
stopifnot(!is.na(all_ds_dt$oiCols))

dotShapes <- c(oi=8, other=20)

dotSizes <- c(oi=6, other=3)
dotSizesBase <- c(oi=2.5, other=1)


all_ds_dt$oiSizes <- dotSizesBase[paste0(all_ds_dt$oiType)]
all_ds_dt$oiShapes <- dotShapes[paste0(all_ds_dt$oiType)]

plotType <- "svg"
myHeight <- 7
myWidth <- 5

myTit <- paste0(hicds_oi, " - ", exprds_oi, "\n", tad_oi)

for(plot_var in c("meanCorr", "adj_meanCorr_empPval_log10", "adj_meanFC_empPval_log10")) {
  
  p <- ggstripchart(all_ds_dt, y= plot_var, xlab=plot_var, color="oiType", size="oiType",shape="oiType") +
    scale_color_manual(values=oiTypes) + 
    scale_size_manual(values=dotSizes) + 
    scale_shape_manual(values=dotShapes) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
    scale_x_discrete(name=tad_oi)+
    labs(color="", size="")+
    guides(color=F, size=F, shape=F)+
    ggtitle(paste0(myTit, " - ", plot_var), subtitle = paste0("# other matches = ", nrow(all_ds_dt) -1)) +
    theme(
          legend.text=element_text(size=14),
          plot.title=element_text(size=12, face="bold", hjust = 0.5),
          plot.subtitle=element_text(size=14),
          axis.text.x=element_blank(),
          axis.title.x=element_text(size=14, face="bold"),
          axis.title.y=element_text(size=14)
    )
  
  # outfile <- file.path(outFolder, paste0(hicds_oi, "_", exprds_oi, "_", tad_oi, "_", plot_var, "_cmp_selectedTAD_acrossDS_colSelected.", plotType))
  # ggsave(p, filename=outfile)
  # cat(paste0("... written:" , outfile, "\n"))
  
  
  p <- ggstripchart(all_ds_dt, y = plot_var, xlab=plot_var, color="cmpType", shape="oiType", size="oiType") +
    scale_color_manual(values=all_cols) + 
    scale_size_manual(values=dotSizes) + 
    scale_shape_manual(values=dotShapes) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
    scale_x_discrete(name=tad_oi)+
    labs(color="", size= "")+
    guides(size=F, shape=F)+
    ggtitle(paste0(myTit, " - ", plot_var), subtitle = paste0("# other matches = ", nrow(all_ds_dt) -1)) +
    theme(          legend.text=element_text(size=14),
                    plot.title=element_text(size=12, face="bold", hjust = 0.5),
                    plot.subtitle=element_text(size=14),
                    axis.text.x=element_blank(),
                    axis.title.x=element_text(size=14, face="bold"),
                    axis.title.y=element_text(size=14)
                    )
  
  outfile <- file.path(outFolder, paste0(hicds_oi, "_", exprds_oi, "_", tad_oi, "_", plot_var, "_cmp_selectedTAD_acrossDS_colCmpTypes.", plotType))
  ggsave(p, filename=outfile)
  cat(paste0("... written:" , outfile, "\n"))
  
}

myWidth <- 6
plotCex <- 1.2
# outfile <- file.path(outFolder, paste0(hicds_oi, "_", exprds_oi, "_", tad_oi, "_", "empPvalCorr_empPvalFC_selectedTAD_acrossDS_colSelected.", plotType))
# do.call(plotType, list(outfile, height=myWidth, width=myWidth))
# plot(x=all_ds_dt$adj_meanFC_empPval_log10,
#      y=all_ds_dt$adj_meanCorr_empPval_log10,
#      col=all_ds_dt$oiCols,
#      pch=all_ds_dt$oiShapes,
#      cex=all_ds_dt$oiSizes,
#      main=myTit,
#      xlab="adj_meanFC_empPval_log10",
#      ylab="adj_meanCorr_empPval_log10",
#      cex.lab=plotCex,
#      cex.axis=plotCex
#      )
# foo <- dev.off()
# cat(paste0("... written:" , outfile, "\n"))

outfile <- file.path(outFolder, paste0(hicds_oi, "_", exprds_oi, "_", tad_oi, "_", "empPvalCorr_empPvalFC_selectedTAD_acrossDS_colCmpTypes.", plotType))
do.call(plotType, list(outfile, height=myWidth, width=myWidth))
plot(x=all_ds_dt$adj_meanFC_empPval_log10,
     y=all_ds_dt$adj_meanCorr_empPval_log10,
     col=all_ds_dt$cmpCols,
     pch=all_ds_dt$oiShapes,
     cex=all_ds_dt$oiSizes,
     main=myTit,
     xlab="adj_meanFC_empPval_log10",
     ylab="adj_meanCorr_empPval_log10",
     cex.lab=plotCex,
     cex.axis=plotCex
     )
legend("topleft", legend=names(all_cols), col = all_cols, bty="n", pch=16)
foo <- dev.off()
cat(paste0("... written:" , outfile, "\n"))





