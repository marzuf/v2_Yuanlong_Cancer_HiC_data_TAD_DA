
# Rscript check_akr1c_del.R

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(reshape2)
require(foreach)
require(doMC)
registerDoMC(40)


setDir <- "/media/electron"
setDir <- ""
load(file.path(setDir, "/mnt/ndata//databank/TCGA/TCGA_PancanAtlas/cnv/pancan_cnv/data_CNA.RData"))

# > cna[1:5,1:5]
# TCGA-OR-A5J1-01 TCGA-OR-A5J2-01 TCGA-OR-A5J3-01 TCGA-OR-A5J4-01 TCGA-OR-A5J5-01
# ACAP3                 0               0               0               1               0
# ACTRT2                0               0               0               1               0
# 
# Hugo_Symbol Entrez_Gene_Id Cytoband
# ACAP3         ACAP3         116983  1p36.33
# ACTRT2       ACTRT2         140625  1p36.32
# AGRN           AGRN         375790  1p36.33

outFolder <- file.path("CHECK_AKR1C_DEL")
dir.create(outFolder, recursive = TRUE)


# args <- commandArgs(trailingOnly = TRUE)
# hicds <- args[1]
# exprds <- args[2]

plotType <- "png"
plotCex <- 1.4
myHeight <- myWidth <- 400

buildData <- TRUE

tad_signif_thresh <- 0.01
signif_col <- "red"
not_signif_col <- "grey"

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))




all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds <- all_hicds[! (grepl("PERMUTG2T", all_hicds) | grepl("RANDOM", all_hicds))]

all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))

if(buildData){
  
  all_cnv_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    # all_cnv_dt <- foreach(hicds = all_hicds[1], .combine='rbind') %dopar%{
    
    
    # if(!grepl("LG1_40kb", hicds)) return(NULL)
    
    
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      if(!grepl("TCGAluad", exprds)) return(NULL)
      
      
      # exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]][1], .combine='rbind') %do% {
      
      pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds)
      
      
      
      result_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds,]
      
      settingFile <- file.path("PIPELINE/INPUT_FILES/", hicds, paste0("run_settings_", exprds, ".R"))
      source(settingFile)
      samp1 <- get(load(file.path(setDir, sample1_file)))
      samp2 <- get(load(file.path(setDir, sample2_file)))
      
      g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
      g2t_dt <- read.delim(g2t_file, header=FALSE, stringsAsFactors = FALSE,
                           col.names=c("entrezID", "chromo", "start", "end","region"))
      g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
      
      pipeline_geneList <- get(load(file.path(pipFolder, "0_prepGeneData", "pipeline_geneList.Rdata")))
      
      hugo_cnv_pip <- cna_meta$Hugo_Symbol[cna_meta$Entrez_Gene_Id %in% pipeline_geneList]
      stopifnot(length(hugo_cnv_pip) > 0)
      
      stopifnot(hugo_cnv_pip %in% rownames(cna))
      
      tmp_dt <- cna_meta[cna_meta$Hugo_Symbol %in% hugo_cnv_pip,]
      stopifnot(!duplicated(tmp_dt$Hugo_Symbol))
      hugo2entrez <- setNames(tmp_dt$Entrez_Gene_Id, tmp_dt$Hugo_Symbol)
      
      cnv_samp1_dt <- cna[rownames(cna) %in% hugo_cnv_pip , colnames(cna) %in% samp1, drop=FALSE]
      cnv_samp2_dt <- cna[rownames(cna) %in% hugo_cnv_pip , colnames(cna) %in% samp2, drop=FALSE]
      stopifnot(rownames(cnv_samp1_dt) %in% names(hugo2entrez))
      stopifnot(rownames(cnv_samp2_dt) %in% names(hugo2entrez))
      
      cnv_samp1_dt <- as.data.frame(cnv_samp1_dt)
      cnv_samp1_dt$entrez <- hugo2entrez[paste0(rownames(cnv_samp1_dt))]
      cnv_samp1_dt$symbol <- rownames(cnv_samp1_dt)
      
      cnv_samp2_dt <- as.data.frame(cnv_samp2_dt)
      cnv_samp2_dt$entrez <- hugo2entrez[paste0(rownames(cnv_samp2_dt))]
      cnv_samp2_dt$symbol <- rownames(cnv_samp2_dt)
      
      nSamp1 <- ncol(cnv_samp1_dt) - 1
      nSamp2 <- ncol(cnv_samp2_dt) - 1
      

      long_dt_samp1 <- melt(cnv_samp1_dt, id=c("entrez", "symbol"))
      long_dt_samp1$cond <- cond1

            
      if(ncol(cnv_samp1_dt) > 2) { # only entrez
        stopifnot("value" %in% colnames(long_dt_samp1))
        stopifnot(cnv_samp1_dt$entrez %in% g2t_dt$entrezID)
        # long_dt_samp1 <- melt(cnv_samp1_dt, id=c("entrez", "symbol"))
        # mean_dt_samp1 <- aggregate(value ~ entrez, data=long_dt_samp1, FUN=mean)
        # stopifnot(mean_dt_samp1$entrez %in% pipeline_geneList)
        
      } else {
        # stop("--eror\n")
        stopifnot(!"value" %in% colnames(long_dt_samp1))
        stopifnot(!"variable" %in% colnames(long_dt_samp1))
        long_dt_samp1$value <- NA
        long_dt_samp1$variable <- NA
      }
      
      long_dt_samp2 <- melt(cnv_samp2_dt, id=c("entrez", "symbol"))
      long_dt_samp2$cond <- cond2
      
      if(ncol(cnv_samp2_dt) > 2) {
        stopifnot("value" %in% colnames(long_dt_samp2))
        stopifnot(cnv_samp2_dt$entrez %in% g2t_dt$entrezID)
        # mean_dt_samp2 <- aggregate(value ~ entrez, data=long_dt_samp2, FUN=mean)
        # stopifnot(mean_dt_samp2$entrez %in% pipeline_geneList)
        
      } else {
        # stop("--eror\n")
        stopifnot(!"value" %in% colnames(long_dt_samp2))
        stopifnot(!"variable" %in% colnames(long_dt_samp2))
        long_dt_samp2$value <- NA
        long_dt_samp2$variable <- NA
      }

      if(ncol(long_dt_samp1) != ncol(long_dt_samp2)){
        cat(hicds, " - ",exprds, "\n")
        save(long_dt_samp1, file="long_dt_samp1.Rdata", version=2)
             save(long_dt_samp2, file="long_dt_samp2.Rdata", version=2)
      }
      mycols <- c("entrez",  "symbol"    ,    "variable", "value", "cond") 
      stopifnot(mycols %in% colnames(long_dt_samp1))
      stopifnot(mycols %in% colnames(long_dt_samp2))
      
      out_dt <- rbind(long_dt_samp1[,mycols], long_dt_samp2[,mycols])
      out_dt$hicds <- hicds
      out_dt$exprds <- exprds
      out_dt
    }# end-for iterating over exprds
    exprds_dt
  } # end-for iterating over hicds
  
  outFile <- file.path(outFolder, "all_cnv_dt.Rdata")
  save(all_cnv_dt , file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_cnv_dt.Rdata")
  all_cnv_dt <- get(load(outFile))
}  

stop("-ok\n")




stopifnot(!duplicated(file.path(all_cnv_dt$hicds, all_cnv_dt$exprds, all_cnv_dt$region)))

signifThresh <- 0.01

all_cnv_dt$diff_cnv_samp12_abs <- abs(all_cnv_dt$diff_cnv_samp12)

all_cnv_dt$signif <- ifelse(all_cnv_dt$adjPvalComb <= signifThresh, "signif.", "not signif.")

require(ggpubr)
require(ggsci)

my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], c("signif.", "not signif."))

plotTit <- "abs. diff. CNV samp1-samp2"
mySub <- paste0("# DS = ", length(unique(file.path(all_cnv_dt$hicds, all_cnv_dt$exprds))))
legTitle <- ""

p3 <- ggdensity(all_cnv_dt,
                x = "diff_cnv_samp12_abs",
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = "TAD mean abs. diff. CNV mean samp1-samp2",
                # add = "median",                  # Add median line.
                rug = FALSE,                      # Add marginal rug
                color = "signif",
                fill = "signif",
                palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  scale_color_manual(values=my_cols)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(
    # text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  ) 

outFile <- file.path(outFolder, paste0("abs_diff_mean_cnv_density.", "svg"))
ggsave(p3, file=outFile, height=5.5, width=7)
cat(paste0("... written: ", outFile, "\n"))

stop("-ok\n")


























if(exists("meanCNV_reg_signif_dt")) rm("meanCNV_reg_signif_dt")

plotTit <- paste0("all datasets (n=", length(unique(file.path(all_cnv_dt$hicds, all_cnv_dt$exprds))), ")")
subTit <- ""

my_x <- all_cnv_dt$diff_cnv_samp12
my_y <- all_cnv_dt$meanLogFC

outFile <- file.path(outFolder, paste0("allDS_diffCNV_tadLogFC_signifCol.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="l")
plot(
  x=my_x,
  y=my_y,
  col = all_cnv_dt$dotCol,
  cex=0.7,
  pch = 16,
  cex.lab=plotCex,
  cex.main=plotCex,
  cex.axis=plotCex,
  xlab = "Mean TAD diff. CNV samp1-samp2",
  ylab = "TAD meanLogFC",
  sub=subTit,
  main = plotTit
)
points(
  x=my_x[all_cnv_dt$dotCol == signif_col],
  y=my_y[all_cnv_dt$dotCol == signif_col],
  col = all_cnv_dt$dotCol[all_cnv_dt$dotCol == signif_col],
  cex=1.2,
  pch = 16
)
legend("topright",
       legend=c(
         paste0("signif. TADs (p-val<=", tad_signif_thresh, ")")
       ),
       pch = 16,
       # col = c(signif_col, not_signif_col),
       col = c(signif_col),
       bty="n"
)
# mtext(side = 3, text = subTit)
abline(lm(my_y~my_x), lty=2, col="grey")
addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
mtext(side = 3, text = paste0("n=", nrow(all_cnv_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("allDS_diffCNV_tadAdjCombPval_signifCol.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

my_x <- all_cnv_dt$diff_cnv_samp12
my_y <- all_cnv_dt$adjPvalComb

par(bty="l")
plot(
  x=my_x,
  y=my_y,
  col=all_cnv_dt$dotCol,
  cex=0.7,
  pch = 16,
  cex.lab=plotCex,
  cex.main=plotCex,
  cex.axis=plotCex,
  xlab = "Mean TAD diff. CNV samp1-samp2",
  ylab = "TAD adj. comb. p-val.",
  sub=subTit,
  main = plotTit
)

points(
  x=my_x[all_cnv_dt$dotCol == signif_col],
  y=my_y[all_cnv_dt$dotCol == signif_col],
  col=all_cnv_dt$dotCol[all_cnv_dt$dotCol == signif_col],
  cex=1.2,
  pch = 16
  
)

abline(lm(my_y~my_x), lty=2, col="grey")
addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
# mtext(side = 3, text = subTit)
mtext(side = 3, text = paste0("n=", nrow(all_cnv_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("allDS_diffCNV_tadAdjCombPval_log10_signifCol.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

my_x <- all_cnv_dt$diff_cnv_samp12
my_y <- -log10(all_cnv_dt$adjPvalComb)

par(bty="l")
plot(
  x=my_x,
  y=my_y,
  col=all_cnv_dt$dotCol,
  cex=0.7,
  pch = 16,
  cex.lab=plotCex,
  cex.main=plotCex,
  cex.axis=plotCex,
  xlab = "Mean TAD diff. CNV samp1-samp2",
  ylab = "TAD adj. comb. p-val. [-log10]",
  sub=subTit,
  main = plotTit
)

points(
  x=my_x[all_cnv_dt$dotCol == signif_col],
  y=my_y[all_cnv_dt$dotCol == signif_col],
  col=all_cnv_dt$dotCol[all_cnv_dt$dotCol == signif_col],
  cex=1.2,
  pch = 16
  
)

abline(lm(my_y~my_x), lty=2, col="grey")
addCorr(x = my_x, y = my_y, bty="n",  legPos="bottomright")
# mtext(side = 3, text = subTit)
mtext(side = 3, text = paste0("n=", nrow(all_cnv_dt)))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))






