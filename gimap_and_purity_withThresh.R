########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript gimap_and_purity_withThresh.R\n"))

script_name <- "gimap_and_purity_withThresh.R"


####### !!!!!!!!!!!  FOR THE MOMENT I DO NOT RESTRICT 
# the datasets in which I look at expression (I do not ensure that is a DS where the conserved region is signif. DA)
# rationale: I want to look if this set of genes is related to purity, irrespective of DA


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

# Rscript gimap_and_purity_withThresh.R
# Rscript gimap_and_purity_withThresh.R EPIC

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- T

myHeight <- 400 
myWidth <- 400
plotType <- "png"
plotCex <- 1.4
mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")


plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 11
plotCex <- 1.2

corMet <- "pearson"


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
} else{
  purity_ds <- ""
}


gimap_dt <- get(load("../MANUSCRIPT_FIGURES/FIG_4/CONSERVED_REGIONS_VIZGG_V2/fig4C_conserved_region_130_genes_plot_dt.Rdata"))
gimap_symbols <- as.character(gimap_dt$symbol)
gimap_symbols <- gimap_symbols[grepl("GIMAP", gimap_symbols)]
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$symbol))
stopifnot(!duplicated(gff_dt$entrezID))
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
gimap_entrez <- symb2entrez[gimap_symbols]
stopifnot(!is.na(gimap_entrez))




if(purity_ds == "") {
  file_suffix <- ""
  purity_file <- file.path("tcga_purity_aran2015.csv")
  purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
  pm <- purity_metrics[1]
  purity_plot_name <- "aran"
  # all the ranks are between 1 and 0
} else if(purity_ds == "EPIC") {
  file_suffix <- "_EPIC"
  purity_file <- file.path("../OLDER_v2_Yuanlong_Cancer_HiC_data_TAD_DA/EPIC_PURITY/all_epic_purity_data.Rdata")
  epic_purity_data <- get(load(purity_file))
  purity_dt <- as.data.frame(do.call(rbind, c(lapply(epic_purity_data, function(x) x[["infiltration_fraction"]]))))
  all_pm_metrics <- colnames(purity_dt) #"Bcells"      "CAFs"        "CD4_Tcells"  "CD8_Tcells"  "Endothelial" "Macrophages" "NKcells"     "otherCells" 
  pm <- "otherCells"
  purity_dt$Sample.ID <- rownames(purity_dt)
  purity_dt$Sample.ID <- gsub("\\.", "-", purity_dt$Sample.ID)
  purity_plot_name <- "EPIC"
} else {
  stop("--invalid purity_ds\n")
}


outFolder <- file.path("GIMAP_AND_PURITY_WITHTHRESH", purity_ds)
dir.create(outFolder, recursive = TRUE)

cat(paste0("!!! > purity metric: ", pm, "\n"))


all_hicds <- list.files(pipFolder)
all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]
all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds
all_ds <- unlist(sapply(names(all_exprds), function(x) file.path(x, all_exprds[[paste0(x)]])))
names(all_ds) <- NULL
ds=all_ds[3]

cat(paste0(">>> purity metric\t=\t", pm, "\n"))

# all_ds=all_ds[1]

if(buildTable){
  
  all_fpkm_gimap_purity_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
    
    cat(paste0("... start: ", ds, "\n"))
    
    hicds <- dirname(ds)
    exprds <- basename(ds)
    
    settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingFile))
    source(settingFile)
    
    samp1 <- get(load(file.path(setDir, sample1_file)))
    samp2 <- get(load(file.path(setDir, sample2_file)))
    
    pur_samp1 <- samp1[samp1 %in% purity_dt$Sample.ID | paste0(samp1, "A") %in% purity_dt$Sample.ID]
    cat(paste0("For ", cond1, " - available samples:\t", length(pur_samp1), "/", length(samp1), "\n"))
    
    pur_samp2 <- samp2[samp2 %in% purity_dt$Sample.ID | paste0(samp2, "A") %in% purity_dt$Sample.ID]
    cat(paste0("For ", cond2, " - available samples:\t", length(pur_samp2), "/", length(samp2), "\n"))
    
    if(length(pur_samp1) == 0 & length(pur_samp2) == 0) return(NULL)
    
    pur2_samp1 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp1  | purity_dt$Sample.ID %in% paste0(samp1, "A") ]
    stopifnot(length(pur2_samp1) == length(pur_samp1))
    
    pur2_samp2 <- purity_dt$Sample.ID[purity_dt$Sample.ID %in% samp2  | purity_dt$Sample.ID %in% paste0(samp2, "A") ]
    stopifnot(length(pur2_samp2) == length(pur_samp2))
    
    stopifnot(setequal(gsub("A$", "", pur2_samp1), pur_samp1))
    stopifnot(setequal(gsub("A$", "", pur2_samp2), pur_samp2))
    
    fpkm_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")
    stopifnot(file.exists(fpkm_file))
    fpkm_dt <- get(load(fpkm_file))  
    # DO THE SAME AS IN THE MANUSCRIPT_FIGURES SCRIPTS (cf discussion marco)
    # changed here 11.03.2020 -> should be normalized sample-wise
    fpkm_dt2 <- apply(fpkm_dt, 2, function(x)x/sum(x))
    # stopifnot(colSums(fpkm_dt2) == 1)
    stopifnot(abs(colSums(fpkm_dt2) - 1) <= 10^-4)
    # and then multiply by 10^6 to have FPKM
    fpkm_dt2 <- fpkm_dt2*10^6
    fpkm_dt2 <- data.frame(fpkm_dt2, check.names = FALSE)
    stopifnot(dim(fpkm_dt) == dim(fpkm_dt2))
    stopifnot(rownames(fpkm_dt) == rownames(fpkm_dt2))
    stopifnot(colnames(fpkm_dt) == colnames(fpkm_dt2))
    fpkm_dt <- fpkm_dt2
    
    gene_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")
    stopifnot(file.exists(gene_file))
    geneList <- get(load(gene_file))
    stopifnot(names(geneList) %in% rownames(fpkm_dt))
    
    # stopifnot(gimap_entrez %in% geneList) # not true for all if gimap not retained
    
    gimap_entrez_fpkm <- names(geneList)[geneList %in% gimap_entrez]
    stopifnot(gimap_entrez_fpkm %in% rownames(fpkm_dt))
    gimap_fpkm_dt <- fpkm_dt[gimap_entrez_fpkm,]
    stopifnot(nrow(gimap_fpkm_dt) == length(gimap_entrez_fpkm))
    
    stopifnot(pur_samp1 %in% colnames(fpkm_dt)); stopifnot(pur_samp2 %in% colnames(fpkm_dt))
    stopifnot(pur2_samp1 %in% purity_dt$Sample.ID); stopifnot(pur2_samp2 %in% purity_dt$Sample.ID)
    
    all_samps <- sort(c(pur_samp1, pur_samp2))
    stopifnot(!duplicated(all_samps))
    all_samps2 <- sort(c(pur2_samp1, pur2_samp2))
    stopifnot(!duplicated(all_samps2))
    stopifnot(gsub("A$", "", all_samps2) == all_samps)
    
    purity_values <- setNames(purity_dt[purity_dt$Sample.ID %in% pur2_samp1 | purity_dt$Sample.ID %in% pur2_samp2,paste0(pm)],
                              purity_dt[purity_dt$Sample.ID %in% pur2_samp1 | purity_dt$Sample.ID %in% pur2_samp2,paste0("Sample.ID")])
    
    stopifnot(setequal(names(purity_values), all_samps2))
    names(purity_values) <- gsub("A$", "", names(purity_values))
    stopifnot(setequal(names(purity_values), all_samps))

    stopifnot(all_samps %in% colnames(gimap_fpkm_dt))
    
    gimap_dt <- data.frame(t(gimap_fpkm_dt[,all_samps]), check.names=FALSE)
    stopifnot(colnames(gimap_dt) %in% names(pipeline_geneList))
    stopifnot(!duplicated(pipeline_geneList))
    # reback to gff dt entrezID -> needed to add column if missing gmaps
    colnames(gimap_dt) <- pipeline_geneList[colnames(gimap_dt)]
    stopifnot(colnames(gimap_dt) %in% pipeline_geneList)
    
    # if there are some GIMAPs that were not available for this dataset -> add NA
    stopifnot(colnames(gimap_dt) %in% gimap_entrez)
    missed_gimaps <- setdiff(gimap_entrez, colnames(gimap_dt))
    for(mg in missed_gimaps) {
      gimap_dt[mg] <- NA
    }
    gimap_dt$sampID <- rownames(gimap_dt)
    rownames(gimap_dt) <- NULL
    purity_dt <- data.frame(
      dataset=ds,
      sampID = names(purity_values),
      purity = purity_values,
      stringsAsFactors = FALSE
    )
    purity_expr_dt <- merge(purity_dt, gimap_dt, by="sampID", all=TRUE)
    
    stopifnot(ncol(purity_expr_dt) == length(gimap_entrez) + 3)
    
    purity_expr_dt <- purity_expr_dt[,c("dataset", "sampID", "purity", gimap_entrez)]
    purity_expr_dt
    
    
  }
  outFile <- file.path(outFolder, "all_fpkm_gimap_purity_dt.Rdata")
  save(all_fpkm_gimap_purity_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_fpkm_gimap_purity_dt.Rdata")
  all_fpkm_gimap_purity_dt <- get(load(outFile))
}

all_fpkm_gimap_purity_dt <- all_fpkm_gimap_purity_dt[all_fpkm_gimap_purity_dt$purity >= -2,]

my_ylab <- "RNA-seq TPM [log10]" 
my_xlab <- "sample purity"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

ge = gimap_entrez[1]

offset_log10 <- 0.001

                            # 
                            # for(ge in gimap_entrez) {
                            #   gs <- entrez2symb[ge]
                            #   stopifnot(gs == names(ge))
                            #   myTit <- paste0(gs)
                            #   tmp_dt <- all_fpkm_gimap_purity_dt[,c("dataset", ge, "purity")] 
                            #   tmp_dt <- na.omit(tmp_dt)
                            #   tmp_dt$dataset <- as.character(tmp_dt$dataset)
                            #   
                            #   my_x <- tmp_dt[,"purity"]
                            #   my_x <- log10(my_x + offset_log10)
                            #   my_y <- tmp_dt[,ge]
                            #   
                            #   legTxt <- paste("# DS=",length(unique(tmp_dt$dataset)), "/",length(unique(all_fpkm_gimap_purity_dt$dataset)), 
                            #                   "\n# values=", nrow(tmp_dt), "/", nrow(all_fpkm_gimap_purity_dt))  
                            #   
                            #   # outfile <- file.path(outFolder, paste0(gs, "_expr_vs_purity.", plotType))
                            #   # do.call(plotType, list(outfile, height=myWidth, width=myWidth))
                            #   outfile <- file.path(outFolder, paste0(gs, "_expr_vs_purity.", "png"))
                            #   do.call("png", list(outfile, height=400, width=400))
                            #   plot(x=my_x,
                            #        y=my_y,
                            #        ylab=my_ylab,
                            #        xlab=my_xlab,
                            #        pch=16,
                            #        cex=0.7,
                            #        cex.lab=plotCex,
                            #        cex.axis=plotCex,
                            #        cex.main=plotCex,
                            #        main=myTit)
                            #   mtext(text=paste0(purity_plot_name, " - ", corMet,"'s corr."), side=3)
                            #   legend("topright", legend=legTxt, bty="n")
                            #   addCorr(my_x, my_y, bty="n", legPos = "topleft", corMet = corMet)
                            #   foo <- dev.off()
                            #   cat(paste0("... written: ", outfile, "\n"))
                            # }

all_ds_corr_dt <- foreach(ge = gimap_entrez, .combine='rbind') %dopar% {
  gs <- entrez2symb[ge]
  stopifnot(gs == names(ge))
  myTit <- paste0(gs)
  tmp_dt <- all_fpkm_gimap_purity_dt[,c("dataset", ge, "purity")] 
  tmp_dt <- na.omit(tmp_dt)
  tmp_dt$dataset <- as.character(tmp_dt$dataset)
  
  my_x <- tmp_dt[,"purity"]

  my_y <- tmp_dt[,ge]
  my_y <- log10(my_y + offset_log10)
  
  legTxt <- paste("# DS=",length(unique(tmp_dt$dataset)), "/",length(unique(all_fpkm_gimap_purity_dt$dataset)), 
                  "\n# values=", nrow(tmp_dt), "/", nrow(all_fpkm_gimap_purity_dt))  
  
  # outfile <- file.path(outFolder, paste0(gs, "_expr_vs_purity.", plotType))
  # do.call(plotType, list(outfile, height=myWidth, width=myWidth))
            outfile <- file.path(outFolder, paste0(gs, "_expr_vs_purity.", "png"))
            do.call("png", list(outfile, height=400, width=400))
            plot(x=my_x,
                 y=my_y,
                 ylab=my_ylab,
                 xlab=my_xlab,
                 pch=16,
                 cex=0.7,
                 cex.lab=plotCex,
                 cex.axis=plotCex,
                 cex.main=plotCex,
                 main=myTit)
            mtext(text=paste0(purity_plot_name, " - ", corMet,"'s corr."), side=3)
            legend("topright", legend=legTxt, bty="n")
            addCorr(my_x, my_y, bty="n", legPos = "topleft", corMet = corMet)
            foo <- dev.off()
            cat(paste0("... written: ", outfile, "\n"))
  ds=unique(tmp_dt$dataset)[1]
  ds_dt <- foreach(ds=unique(tmp_dt$dataset), .combine = 'rbind') %dopar% {
    ds_tmp_dt <- tmp_dt[tmp_dt$dataset == ds,]
    ds_x <- ds_tmp_dt[,"purity"]
    ds_y <- ds_tmp_dt[,ge]
    data.frame(
      dataset=ds,
      nSamp = nrow(ds_tmp_dt),
      entrezID=ge,
      symbol=gs,
      corrExprPurity=cor(ds_x, ds_y, method=corMet),
      stringsAsFactors = FALSE)
  }
  ds_dt
}
outFile <- file.path(outFolder, "all_ds_corr_dt.Rdata")
save(all_ds_corr_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

stopifnot(colnames(all_fpkm_gimap_purity_dt) == c("dataset", "sampID", "purity", gimap_entrez))
all_fpkm_gimap_purity_dt$meanExpr <- rowMeans( all_fpkm_gimap_purity_dt[,gimap_entrez], na.rm=TRUE)

my_x <- all_fpkm_gimap_purity_dt[,"purity"]
my_y <- all_fpkm_gimap_purity_dt[,"meanExpr"]

my_y <- log10(my_y+offset_log10)

myTit <- paste0("all GIMAPs (meanExpr)")

# outfile <- file.path(outFolder, paste0("allGIMAPs_meanExpr_vs_purity.", plotType))
# do.call(plotType, list(outfile, height=myWidth, width=myWidth))
outfile <- file.path(outFolder, paste0("allGIMAPs_meanExpr_vs_purity.", "png"))
do.call("png", list(outfile, height=400, width=400))
plot(x=my_x,
     y=my_y,
     ylab=my_ylab,
     xlab=my_xlab,
     pch=16,
     cex=0.7,
     cex.lab=plotCex,
     cex.axis=plotCex,
     cex.main=plotCex,
     main=myTit)
mtext(text=paste0(purity_plot_name, " - ", corMet,"'s corr."), side=3)
addCorr(my_x, my_y, bty="n", legPos = "topright", corMet=corMet)
foo <- dev.off()
cat(paste0("... written: ", outfile, "\n"))

all_fpkm_gimap_purity_dt$dataset <- as.character(all_fpkm_gimap_purity_dt$dataset)
ds=unique(all_fpkm_gimap_purity_dt$dataset)[1]
all_ds_meanExpr_corr_dt <- foreach(ds=unique(all_fpkm_gimap_purity_dt$dataset), .combine = 'rbind') %dopar% {
  ds_tmp_dt <- all_fpkm_gimap_purity_dt[all_fpkm_gimap_purity_dt$dataset == ds, c("purity", "meanExpr")]
  ds_tmp_dt <- na.omit(ds_tmp_dt)
  ds_x <- ds_tmp_dt[,"purity"]
  ds_y <- ds_tmp_dt[,"meanExpr"]
    data.frame(
    dataset=ds,
    nSamp = nrow(ds_tmp_dt),
    corrExprPurity=cor(ds_x, ds_y, method=corMet),
    stringsAsFactors = FALSE)
}
outFile <- file.path(outFolder, "all_ds_meanExpr_corr_dt.Rdata")
save(all_ds_meanExpr_corr_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


all_ds_meanExpr_corr_dt$entrezID <- NA
all_ds_meanExpr_corr_dt$symbol <- "meanExpr"
all_ds_meanExpr_corr_dt <- all_ds_meanExpr_corr_dt[,colnames(all_ds_corr_dt)]

plot_dt <- rbind(all_ds_meanExpr_corr_dt, all_ds_corr_dt)

#### HERE THE BOXPLOTS !!!
fontFamily <- "Hershey"
my_box_theme <- theme(
  text = element_text(family=fontFamily),
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

myTit <- paste0("Correlation GIMAP expression and sample purity")
subTit <- paste0("by GIMAP and dataset - ", purity_plot_name, " - ", corMet,"'s corr.")
xlab <- ""
ylab <- "Correlation between GIMAP expression and sample purity"

p_boxplot <- ggplot(data=plot_dt, aes(x=symbol, y=corrExprPurity))+
    geom_boxplot(notch = TRUE, outlier.shape=NA)+
  geom_jitter()+
  ggtitle(myTit, subtitle = paste0(subTit))+
  scale_x_discrete(name=xlab)+
  scale_y_continuous(name=paste0(ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  my_box_theme

outFile <- file.path(outFolder, paste0("all_corr_expr_purity_GIMAPs_byDS_boxplot.", plotType))
ggsave(plot = p_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

plot_dt$cmpType <- all_cmps[basename(plot_dt$dataset)]
stopifnot(!is.na(plot_dt$cmpType))

p_boxplot <- ggplot(data=plot_dt, aes(x=symbol, y=corrExprPurity))+
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  geom_jitter(aes(color=cmpType))+
  ggtitle(myTit, subtitle = paste0(subTit))+
  scale_x_discrete(name=xlab)+
  scale_y_continuous(name=paste0(ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  my_box_theme

outFile <- file.path(outFolder, paste0("all_corr_expr_purity_GIMAPs_byDS_cmpCol_boxplot.", plotType))
ggsave(plot = p_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


