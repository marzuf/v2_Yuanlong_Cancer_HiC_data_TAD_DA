########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript hoxb_and_purity.R\n"))

script_name <- "hoxb_and_purity.R"


####### !!!!!!!!!!!  FOR THE MOMENT I DO NOT RESTRICT 
# the datasets in which I look at expression (I do not ensure that is a DS where the conserved region is signif. DA)
# rationale: I want to look if this set of genes is related to purity, irrespective of DA


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

# Rscript hoxb_and_purity.R
# Rscript hoxb_and_purity.R EPIC

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


# hoxb_dt <- get(load("../MANUSCRIPT_FIGURES/FIG_4/CONSERVED_REGIONS_VIZGG_V2/fig4C_conserved_region_42_genes_plot_dt.Rdata"))
# hoxb_symbols <- as.character(hoxb_dt$symbol)
# hoxb_symbols <- hoxb_symbols[grepl("HOXB", hoxb_symbols)]
hoxb_symbols <- paste0("HOXB",2:7)
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$symbol))
stopifnot(!duplicated(gff_dt$entrezID))
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
hoxb_entrez <- symb2entrez[hoxb_symbols]
stopifnot(!is.na(hoxb_entrez))




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


outFolder <- file.path("HOXB_AND_PURITY", purity_ds)
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
  
  all_fpkm_hoxb_purity_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
    
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
    
    # stopifnot(hoxb_entrez %in% geneList) # not true for all if gimap not retained
    
    hoxb_entrez_fpkm <- names(geneList)[geneList %in% hoxb_entrez]
    stopifnot(hoxb_entrez_fpkm %in% rownames(fpkm_dt))
    hoxb_fpkm_dt <- fpkm_dt[hoxb_entrez_fpkm,]
    stopifnot(nrow(hoxb_fpkm_dt) == length(hoxb_entrez_fpkm))
    
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

    stopifnot(all_samps %in% colnames(hoxb_fpkm_dt))
    
    hoxb_dt <- data.frame(t(hoxb_fpkm_dt[,all_samps]), check.names=FALSE)
    stopifnot(colnames(hoxb_dt) %in% names(pipeline_geneList))
    stopifnot(!duplicated(pipeline_geneList))
    # reback to gff dt entrezID -> needed to add column if missing gmaps
    colnames(hoxb_dt) <- pipeline_geneList[colnames(hoxb_dt)]
    stopifnot(colnames(hoxb_dt) %in% pipeline_geneList)
    
    # if there are some GIMAPs that were not available for this dataset -> add NA
    stopifnot(colnames(hoxb_dt) %in% hoxb_entrez)
    missed_hoxbs <- setdiff(hoxb_entrez, colnames(hoxb_dt))
    for(mg in missed_hoxbs) {
      hoxb_dt[mg] <- NA
    }
    hoxb_dt$sampID <- rownames(hoxb_dt)
    rownames(hoxb_dt) <- NULL
    purity_dt <- data.frame(
      dataset=ds,
      sampID = names(purity_values),
      purity = purity_values,
      stringsAsFactors = FALSE
    )
    purity_expr_dt <- merge(purity_dt, hoxb_dt, by="sampID", all=TRUE)
    
    stopifnot(ncol(purity_expr_dt) == length(hoxb_entrez) + 3)
    
    purity_expr_dt <- purity_expr_dt[,c("dataset", "sampID", "purity", hoxb_entrez)]
    purity_expr_dt
    
    
  }
  outFile <- file.path(outFolder, "all_fpkm_hoxb_purity_dt.Rdata")
  save(all_fpkm_hoxb_purity_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_fpkm_hoxb_purity_dt.Rdata")
  all_fpkm_hoxb_purity_dt <- get(load(outFile))
}

my_ylab <- "RNA-seq TPM [log10]" 
my_xlab <- "sample purity"

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

ge = hoxb_entrez[1]

offset_log10 <- 0.001

                            # 
                            # for(ge in hoxb_entrez) {
                            #   gs <- entrez2symb[ge]
                            #   stopifnot(gs == names(ge))
                            #   myTit <- paste0(gs)
                            #   tmp_dt <- all_fpkm_hoxb_purity_dt[,c("dataset", ge, "purity")] 
                            #   tmp_dt <- na.omit(tmp_dt)
                            #   tmp_dt$dataset <- as.character(tmp_dt$dataset)
                            #   
                            #   my_x <- tmp_dt[,"purity"]
                            #   my_x <- log10(my_x + offset_log10)
                            #   my_y <- tmp_dt[,ge]
                            #   
                            #   legTxt <- paste("# DS=",length(unique(tmp_dt$dataset)), "/",length(unique(all_fpkm_hoxb_purity_dt$dataset)), 
                            #                   "\n# values=", nrow(tmp_dt), "/", nrow(all_fpkm_hoxb_purity_dt))  
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

all_ds_corr_dt <- foreach(ge = hoxb_entrez, .combine='rbind') %dopar% {
  gs <- entrez2symb[ge]
  stopifnot(gs == names(ge))
  myTit <- paste0(gs)
  tmp_dt <- all_fpkm_hoxb_purity_dt[,c("dataset", ge, "purity")] 
  tmp_dt <- na.omit(tmp_dt)
  tmp_dt$dataset <- as.character(tmp_dt$dataset)
  
  my_x <- tmp_dt[,"purity"]

  my_y <- tmp_dt[,ge]
  my_y <- log10(my_y + offset_log10)
  
  legTxt <- paste("# DS=",length(unique(tmp_dt$dataset)), "/",length(unique(all_fpkm_hoxb_purity_dt$dataset)), 
                  "\n# values=", nrow(tmp_dt), "/", nrow(all_fpkm_hoxb_purity_dt))  
  
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

stopifnot(colnames(all_fpkm_hoxb_purity_dt) == c("dataset", "sampID", "purity", hoxb_entrez))
all_fpkm_hoxb_purity_dt$meanExpr <- rowMeans( all_fpkm_hoxb_purity_dt[,hoxb_entrez], na.rm=TRUE)

my_x <- all_fpkm_hoxb_purity_dt[,"purity"]
my_y <- all_fpkm_hoxb_purity_dt[,"meanExpr"]

my_y <- log10(my_y+offset_log10)

myTit <- paste0("all HOXB (meanExpr)")

# outfile <- file.path(outFolder, paste0("allGIMAPs_meanExpr_vs_purity.", plotType))
# do.call(plotType, list(outfile, height=myWidth, width=myWidth))
outfile <- file.path(outFolder, paste0("allHOXBs_meanExpr_vs_purity.", "png"))
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

all_fpkm_hoxb_purity_dt$dataset <- as.character(all_fpkm_hoxb_purity_dt$dataset)
ds=unique(all_fpkm_hoxb_purity_dt$dataset)[1]
all_ds_meanExpr_corr_dt <- foreach(ds=unique(all_fpkm_hoxb_purity_dt$dataset), .combine = 'rbind') %dopar% {
  ds_tmp_dt <- all_fpkm_hoxb_purity_dt[all_fpkm_hoxb_purity_dt$dataset == ds, c("purity", "meanExpr")]
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

myTit <- paste0("Correlation HOXB expression and sample purity")
subTit <- paste0("by HOXB and dataset - ", purity_plot_name, " - ", corMet,"'s corr.")
xlab <- ""
ylab <- "Correlation between HOXB expression and sample purity"

p_boxplot <- ggplot(data=plot_dt, aes(x=symbol, y=corrExprPurity))+
    geom_boxplot(notch = TRUE, outlier.shape=NA)+
  geom_jitter()+
  ggtitle(myTit, subtitle = paste0(subTit))+
  scale_x_discrete(name=xlab)+
  scale_y_continuous(name=paste0(ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  my_box_theme

outFile <- file.path(outFolder, paste0("all_corr_expr_purity_HOXBs_byDS_boxplot.", plotType))
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

outFile <- file.path(outFolder, paste0("all_corr_expr_purity_HOXBs_byDS_cmpCol_boxplot.", plotType))
ggsave(plot = p_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


# load("hoxb_AND_PURITY/all_ds_corr_dt.Rdata")
# load("hoxb_AND_PURITY/all_ds_meanExpr_corr_dt.Rdata")
# #################################################################
# ################################################################# # barplot all agg.
# #################################################################
# 
# 
# all_cors <- do.call(rbind, by(all_fpkmGO_purity_dt, all_fpkmGO_purity_dt$go_id, function(sub_dt){
#   
#   my_x <- sub_dt$geneExpr
#   my_y <- sub_dt$samplePurity
#   mycor <- cor(my_x[!is.na(my_x) & !is.na(my_y)], 
#                my_y[!is.na(my_x) & !is.na(my_y)], method = corMet)
#   
#   my_x <- log10(sub_dt$geneExpr)
#   my_y <- sub_dt$samplePurity
#   mycor_log10 <- cor(my_x[!is.na(my_x) & !is.na(my_y) & is.finite(my_x) & is.finite(my_y)], 
#                      my_y[!is.na(my_x) & !is.na(my_y) & is.finite(my_x) & is.finite(my_y)], method = corMet)
#   
#   curr_go <- unique(sub_dt$go_id)
#   curr_go_rank <- unique(sub_dt$go_rank)
#   
#   data.frame(
#     GO_id=curr_go,
#     GO_rank=curr_go_rank,
#     corr=mycor,
#     corr_log10=mycor_log10,
#     stringsAsFactors = FALSE
#   )
#   
# } ))   
# outFile <- file.path(outFolder, "all_cors.Rdata")
# save(all_cors, file=outFile, version=2)
# cat(paste0("... written: ", outFile, "\n"))
# # load("IMMUNEGO_SAMPLEINF/all_cors.Rdata")
# 
# all_cors <- all_cors[order(all_cors$GO_rank),]
# 
# all_cors$GO_id_labs <- unlist(lapply(strwrap(gsub("_", " ", as.character(all_cors$GO_id)),
#                                              width = strwdth, simplify=FALSE),
#                                      function(x) paste0(x, collapse="\n")))
# all_cors$GO_id <- factor(all_cors$GO_id, levels=all_cors$GO_id)
# all_cors$GO_id_labs <- factor(all_cors$GO_id_labs, levels=all_cors$GO_id_labs)
# stopifnot(!is.na(all_cors$GO_id))
# stopifnot(!is.na(all_cors$GO_id_labs))
# 
# all_cors$corr_abs <- abs(all_cors$corr)
# all_cors$corr_log10_abs <- abs(all_cors$corr_log10)
# 
# my_theme <-     theme(
#   # plot.margin = margin(b = 2, unit = "cm"),
#   text = element_text(family=fontFamily),
#   panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
#   panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
#   panel.background = element_rect(fill = "transparent"),
#   panel.grid.major.x =  element_blank(),
#   panel.grid.minor.x =  element_blank(),
#   axis.line=element_line(),
#   axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
#   axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
#   axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
#   axis.text.x = element_text(size=8, hjust=1, vjust=0.5, angle=90),
#   plot.title = element_text(hjust=0.5, size = 16, face="bold"),
#   plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
#   legend.title = element_text(face="bold")
# )
# 
# 
# var_suff=""
# for(var_suff in c("", "_log10", "_abs", "_log10_abs")){
#   
#   plotTit <- paste0("Corr. GO gene expr. and purity overall DS")
#   subTit <- paste0(gsub("_", " ", var_suff, " - GO signif. thresh <= ", goSignifThresh))
#   
#   bar_p <- ggplot(all_cors, aes_string(x="GO_id_labs",y=paste0("corr", var_suff)))+
#     ggtitle(plotTit, subtitle=subTit)+
#     geom_bar(stat="identity")+
#     labs(x="Ranked GOs", y =paste0(var_suff, " ", corMet, "'s correlation"))+
#     my_theme
#   
#   outFile <- file.path(outFolder, paste0("GO_geneExpr_samplePurity_", corMet, "Corr", "_", var_suff, "_barplot.", plotType))
#   ggsave(bar_p, filename=outFile, height = myHeightGG, width=myWidthGG)
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   
#   
# }
# #################################################################
# ################################################################# # boxplot by DS & boxplot by GO
# #################################################################
# 
# 
# all_cors_byDS <- do.call(rbind, by(all_fpkmGO_purity_dt, all_fpkmGO_purity_dt$go_id, function(sub_dt){
#   
#   sub_dt$dataset <- file.path(sub_dt$hicds, sub_dt$exprds)
#   
#   out_dt <- foreach(ds= unique(as.character(sub_dt$dataset)), .combine='rbind') %dopar% {
#     
#     sub_dt2 <- sub_dt[sub_dt$dataset == ds,]
#     stopifnot(nrow(sub_dt2) > 0)
#     
#     my_x <- sub_dt2$geneExpr
#     my_y <- sub_dt2$samplePurity
#     mycor <- cor(my_x[!is.na(my_x) & !is.na(my_y)], 
#                  my_y[!is.na(my_x) & !is.na(my_y)], method = corMet)
#     
#     my_x <- log10(sub_dt2$geneExpr)
#     my_y <- sub_dt2$samplePurity
#     mycor_log10 <- cor(my_x[!is.na(my_x) & !is.na(my_y) & is.finite(my_x) & is.finite(my_y)], 
#                        my_y[!is.na(my_x) & !is.na(my_y) & is.finite(my_x) & is.finite(my_y)], method = corMet)
#     
#     curr_go <- unique(sub_dt2$go_id)
#     curr_go_rank <- unique(sub_dt2$go_rank)
#     
#     data.frame(
#       dataset = ds,
#       GO_id=curr_go,
#       GO_rank=curr_go_rank,
#       corr=mycor,
#       corr_log10=mycor_log10,
#       stringsAsFactors = FALSE
#     )
#   }
#   out_dt
#   
# } ))   
# rownames(all_cors_byDS) <- NULL
# outFile <- file.path(outFolder, "all_cors_byDS.Rdata")
# save(all_cors_byDS, file=outFile, version=2)
# cat(paste0("... written: ", outFile, "\n"))
# 
# # load("IMMUNEGO_SAMPLEINF/all_cors_byDS.Rdata")
# 
# tmp <- all_cors_byDS[,c("GO_rank", "GO_id")]
# tmp <- unique(tmp)
# tmp$GO_id_labs <- unlist(lapply(strwrap(gsub("_", " ", as.character(tmp$GO_id)),
#                                         width = strwdth, simplify=FALSE),
#                                 function(x) paste0(x, collapse="\n")))
# tmp <- tmp[order(tmp$GO_rank),]
# 
# all_cors_byDS <- all_cors_byDS[order(all_cors_byDS$GO_rank),]
# 
# all_cors_byDS$GO_id_labs <- unlist(lapply(strwrap(gsub("_", " ", as.character(all_cors_byDS$GO_id)),
#                                                   width = strwdth, simplify=FALSE),
#                                           function(x) paste0(x, collapse="\n")))
# all_cors_byDS$GO_id <- factor(all_cors_byDS$GO_id, levels=unique(tmp$GO_id))
# all_cors_byDS$GO_id_labs <- factor(all_cors_byDS$GO_id_labs, levels=tmp$GO_id_labs)
# stopifnot(!is.na(all_cors_byDS$GO_id))
# stopifnot(!is.na(all_cors_byDS$GO_id_labs))
# 
# all_cors_byDS$corr_abs <- abs(all_cors_byDS$corr)
# all_cors_byDS$corr_log10_abs <- abs(all_cors_byDS$corr_log10)
# 
# tmp2 <- aggregate(corr_abs~dataset, FUN=mean, data=all_cors_byDS)
# tmp2 <- tmp2[order(tmp2$corr, decreasing = TRUE),]
# all_cors_byDS$dataset_lab <- paste0(dirname(all_cors_byDS$dataset), "\n", basename(all_cors_byDS$dataset))
# 
# all_cors_byDS$dataset <- factor(all_cors_byDS$dataset, levels=as.character(tmp2$dataset))
# stopifnot(!is.na(all_cors_byDS$dataset))
# 
# all_cors_byDS$dataset_lab <- factor(all_cors_byDS$dataset_lab, levels=paste0(dirname(as.character(tmp2$dataset)),                                                                              "\n", basename(as.character(tmp2$dataset))))
# stopifnot(!is.na(all_cors_byDS$dataset_lab))
# all_cors_byDS$cmpType <- all_cmps[paste0(basename(as.character(all_cors_byDS$dataset)))]
# stopifnot(!is.na(all_cors_byDS$cmpType))
# 
# var_suff=""
# for(var_suff in c("", "_log10", "_abs", "_log10_abs")){
#   
#   plotTit <- paste0("Corr. GO gene expr. and purity by DS")
#   subTit <- paste0(gsub("_", " ", var_suff, " - GO signif. thresh <= ", goSignifThresh))
#   
#   box_p <- ggplot(all_cors_byDS, aes_string(x="GO_id_labs",y=paste0("corr", var_suff)))+
#     ggtitle(plotTit, subtitle=subTit)+
#     geom_boxplot()+
#     labs(x="Ranked GOs", y =paste0(var_suff, " ", corMet, "'s correlation"))+
#     my_theme
#   
#   outFile <- file.path(outFolder, paste0("GO_geneExpr_samplePurity_", corMet, "Corr", "_", var_suff, "_byDS_byGO_boxplot.", plotType))
#   ggsave(box_p, filename=outFile, height = myHeightGG, width=myWidthGG)
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   box_p_cmp <- ggplot(all_cors_byDS, aes_string(x="GO_id_labs",y=paste0("corr", var_suff), color="cmpType"))+
#     ggtitle(plotTit, subtitle=subTit)+
#     geom_boxplot()+
#     labs(x="Ranked GOs", y =paste0(var_suff, " ", corMet, "'s correlation"))+
#     my_theme
#   
#   outFile <- file.path(outFolder, paste0("GO_geneExpr_samplePurity_", corMet, "Corr", "_", var_suff, "_byDS_byGO_byCmp_boxplot.", plotType))
#   ggsave(box_p_cmp, filename=outFile, height = myHeightGG, width=myWidthGG)
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   
#   box_p_ds <- ggplot(all_cors_byDS, aes_string(x="dataset_lab",y=paste0("corr", var_suff)))+
#     ggtitle(plotTit, subtitle=subTit)+
#     geom_boxplot()+
#     labs(x="", y =paste0(var_suff, " ", corMet, "'s correlation"))+
#     my_theme
#   
#   outFile <- file.path(outFolder, paste0("GO_geneExpr_samplePurity_", corMet, "Corr", "_", var_suff, "_byDS_byDS_boxplot.", plotType))
#   ggsave(box_p_ds, filename=outFile, height = myHeightGG, width=myWidthGG)
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   
#   
#   
# }
# 
# 
# 
# 
# #     
# #     
# #     purity_values <- setNames(purity_dt[purity_dt$Sample.ID %in% pur2_samp1 | purity_dt$Sample.ID %in% pur2_samp2,paste0(pm)],
# #                               purity_dt[purity_dt$Sample.ID %in% pur2_samp1 | purity_dt$Sample.ID %in% pur2_samp2,paste0("Sample.ID")])
# #     
# #     names(purity_values) <- gsub("A$", "", names(purity_values))
# #     pur2_samp1 <- gsub("A$", "", pur2_samp1)
# #     pur2_samp2 <- gsub("A$", "", pur2_samp2)
# #     stopifnot(setequal(names(purity_values), c(pur2_samp1, pur2_samp2)))
# #     
# #     stopifnot(setequal(names(purity_values), c(pur_samp1, pur_samp2)))
# #     stopifnot(length(purity_values) == length(c(pur_samp1, pur_samp2)))
# #     
# #     purity_values <- purity_values[c(pur2_samp1, pur2_samp2)]
# #     fpkm_dt <- fpkm_dt[,c(pur_samp1, pur_samp2)]
# #     
# #     fpkm_purity_dt <- fpkm_dt
# #     fpkm_purity_dt[(nrow(fpkm_purity_dt) + 1),] <- purity_values
# #     rownames(fpkm_purity_dt)[nrow(fpkm_purity_dt)] <- "purity"
# #     stopifnot(rownames(fpkm_purity_dt)[nrow(fpkm_purity_dt)] == "purity")
# #     
# #     cat(paste0("... ", hicds, " - ", exprds, " \t ", "corr. for all samples\n"))
# #     
# #     all_partial_corr_dt <- get_partial_corr_dt(fpkmdt_with_pur=fpkm_purity_dt, cormet=corMet, newcol="partial_coexpr")
# #     all_corr_dt <- get_full_corr_dt(fpkmdt=fpkm_dt, cormet=corMet, newcol="coexpr")
# #     
# #     all_dt <- merge(all_partial_corr_dt, all_corr_dt, by=c("gene1", "gene2")) 
# #     
# #     cat(paste0("... ", hicds, " - ", exprds, " \t ", "corr. for samp1\n"))    
# #     if(length(pur_samp1) > 1) {
# #       samp1_partial_corr_dt <- get_partial_corr_dt(fpkmdt_with_pur=fpkm_purity_dt[,c(pur_samp1)], cormet=corMet, newcol="partial_coexpr_samp1")
# #       samp1_corr_dt <- get_full_corr_dt(fpkmdt=fpkm_dt[,c(pur_samp1)], cormet=corMet, newcol="coexpr_samp1")
# #       samp1_dt <- merge(samp1_partial_corr_dt, samp1_corr_dt, by=c("gene1", "gene2")) 
# #       out_dt <- merge(all_dt, samp1_dt, by=c("gene1", "gene2")) 
# #     } else {
# #       out_dt <- all_dt
# #       out_dt$partial_coexpr_samp1 <- out_dt$coexpr_samp1 <- NA
# #     }
# #     
# #     cat(paste0("... ", hicds, " - ", exprds, " \t ", "corr. for samp2\n"))
# #     if(length(pur_samp2) > 1) {
# #       samp2_partial_corr_dt <- get_partial_corr_dt(fpkmdt_with_pur=fpkm_purity_dt[,c(pur_samp2)], cormet=corMet, newcol="partial_coexpr_samp2")
# #       samp2_corr_dt <- get_full_corr_dt(fpkmdt=fpkm_dt[,c(pur_samp2)], cormet=corMet, newcol="coexpr_samp2")
# #       samp2_dt <- merge(samp2_partial_corr_dt, samp2_corr_dt, by=c("gene1", "gene2")) 
# #       out_dt <- merge(out_dt, samp2_dt, by=c("gene1", "gene2")) 
# #     } else {
# #       out_dt$partial_coexpr_samp2 <- out_dt$coexpr_samp2 <- NA
# #     }
# #     
# #     out_dt$hicds <- hicds
# #     out_dt$exprds <- exprds
# #     
# #     out_dt
# #   } # end-iterating ds 
# #   
# #   
# #   outFile <- file.path(outFolder, "coexpr_and_purity_dt.Rdata")
# #   save(coexpr_and_purity_dt, file=outFile, version=2)
# #   cat(paste0("... written: ", outFile, "\n"))
# #   
# #   
# # } else {
# #   outFile <- file.path(outFolder, "coexpr_and_purity_dt.Rdata")
# #   # outFile <- file.path(outFolder, "sub_coexpr_and_purity_dt.Rdata")
# #   cat(paste0("... load coexpr data - ", Sys.time()))
# #   coexpr_and_purity_dt <- get(load(outFile))
# #   cat(paste0(" - ", Sys.time(), " - done\n"))
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # >   purity_file <- file.path("../OLDER_v2_Yuanlong_Cancer_HiC_data_TAD_DA/EPIC_PURITY//all_epic_purity_data.Rdata")
# # epic_purity_data <- get(load(purity_file))
# # [B] 7.6.1810  <ectron:/mnt/etemp/marie/MANUSCRIPT_FIGURES/FIG_2
# # 
# # > str(epic_purity_data[[1]][["infiltration_fraction"]])                                                                                                                                                     
# #  num [1:670, 1:8] 0.003125 0.006525 0.000781 0.007518 0.005236 ...
# #  - attr(*, "dimnames")=List of 2
# #   ..$ : chr [1:670] "TCGA-3C-AAAU-01" "TCGA-3C-AALK-01" "TCGA-4H-AAAK-01" "TCGA-5L-AAT0-01" ...
# #   ..$ : chr [1:8] "Bcells" "CAFs" "CD4_Tcells" "CD8_Tcells" ...
# # > (epic_purity_data[[1]][["infiltration_fraction"]][1:5,1:5])                                                                                                                                               
# #                       Bcells       CAFs CD4_Tcells   CD8_Tcells Endothelial
# # TCGA-3C-AAAU-01 0.0031254992 0.07052641 0.01975444 6.338933e-03  0.01984755
# # TCGA-3C-AALK-01 0.0065250394 0.40216937 0.03330637 2.264790e-06  0.02735698
# # TCGA-4H-AAAK-01 0.0007808537 0.46954227 0.02861602 2.638945e-05  0.01637793
# # TCGA-5L-AAT0-01 0.0075176516 0.26169569 0.04759900 1.154060e-02  0.02780030
# # 
# # 
# # aran_purity_file <- file.path("tcga_purity_aran2015.csv")
# # aran_purity_dt <- read.delim(aran_purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
# # purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
# # pm <- purity_metrics[1]
# # aran_purity_dt$Sample.ID <- gsub("A$", "", aran_purity_dt$Sample.ID)
# # 
# # epic_purity_file <- file.path("EPIC_PURITY/all_epic_purity_data.Rdata")
# # epic_purity_data <- get(load(epic_purity_file))
# # 
# # samples_dt <- do.call(rbind, lapply(1:length(epic_purity_data), function(x) {
# #   tmp <- rownames(epic_purity_data[[x]][["infiltration_fraction"]])
# #   data.frame(
# #     hicds = dirname(names(epic_purity_data)[x]),
# #     exprds = basename(names(epic_purity_data)[x]),
# #     Sample.ID = tmp,
# #     stringsAsFactors = FALSE
# #   )
# # }))
# # 
# # epic_purity_dt <- as.data.frame(do.call(rbind, c(lapply(epic_purity_data, function(x) x[["infiltration_fraction"]]))))
# #                                                  
# # epic_purity_dt$Sample.ID <- rownames(epic_purity_dt)
# # 
# # stopifnot(any(epic_purity_dt$Sample.ID %in% aran_purity_dt$Sample.ID))
# # 
# # 
# # 
# # go_result_dt <- get(load("../MANUSCRIPT_FIGURES/FIG_4/GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
# # nrow(go_result_dt)
# # go_signif_dt <- go_result_dt[go_result_dt$p.adjust <= goSignifThresh,]
# # nrow(go_signif_dt)
# # go_signif_dt$p.adjust_rank <- rank(go_signif_dt$p.adjust, tie=tieMeth)
# # 
# # all_gos <- as.character(go_signif_dt$ID)
# # curr_go = all_gos[1]
# # 
# # go_stat_dt <- foreach(curr_go = all_gos, .combine='rbind') %dopar% {
# #   go_rank <- unique(go_signif_dt$p.adjust_rank[as.character(go_signif_dt$ID) == curr_go])
# #   stopifnot(length(go_rank) == 1)
# #   curr_entrezID <- as.character(go_signif_dt$geneID[as.character(go_signif_dt$ID) == curr_go])
# #   stopifnot(length(curr_entrezID) == 1)
# #   go_entrezID <- unlist(strsplit(x=curr_entrezID, split="/"))
# #   stopifnot(all(go_entrezID %in% tadSignif_inDT$entrezID))
# #   go_symbol <- entrez2symb[paste0(go_entrezID)]
# #   go_inDT <- tadSignif_inDT[tadSignif_inDT$entrezID %in% go_entrezID,]
# #   occByGenes_dt <- aggregate(dataset ~ entrezID, data = go_inDT, FUN=length)
# #   nSignifByGenes_dt <- aggregate(adj.P.Val ~ entrezID, data = go_inDT, FUN=function(x) sum(x<=geneSignifThresh))
# #   occ_signif <- merge(occByGenes_dt, nSignifByGenes_dt, by="entrezID", all=T)
# #   stopifnot(!is.na(occ_signif))
# #   occ_signif$signifRatio <- occ_signif$adj.P.Val/occ_signif$dataset
# #   stopifnot(occ_signif$signifRatio>=0 & occ_signif$signifRatio <= 1)
# #   meanOccByGenes <- mean(occByGenes_dt$dataset)
# #   meanSignifByGenes <- mean(nSignifByGenes_dt$adj.P.Val)
# #   meanRatioSiginf <- mean(occ_signif$signifRatio)
# #   meanFC <- mean(go_inDT$logFC)
# #   median_tadRank <- median(go_inDT$tad_rank)
# #   median_geneRank <- median(go_inDT$gene_rank)
# #   data.frame(
# #     go_id = curr_go,
# #     go_rank=go_rank,
# #     median_tadRank=median_tadRank,
# #     median_geneRank=median_geneRank,
# #     meanFC=meanFC,
# #     meanOccByGenes=meanOccByGenes,
# #     meanSignifByGenes=meanSignifByGenes,
# #     meanRatioSiginf=meanRatioSiginf,
# #     go_genes=paste0(go_symbol, collapse=","),
# #     stringsAsFactors = FALSE
# #   )
# # }
# # 
# # outFile <- file.path(outFolder, "go_stat_dt.Rdata")
# # save(go_stat_dt, file=outFile, version=2)
# # cat(paste0("... written: ", outFile, "\n"))
# # # load("LOOK_CONSERVEDREG_IMMUNEGO/go_stat_dt.Rdata")
# # 
# # go_stat_dt <- go_stat_dt[order(go_stat_dt$go_rank),]
# # plot_dt <- go_stat_dt
# # plot_dt$go_id_lab <- unlist(lapply(strwrap(gsub("_", " ", plot_dt$go_id), 
# #                                            width = strwdth, simplify=FALSE), 
# #                                    function(x) paste0(x, collapse="\n")))
# # plot_dt$go_id <- factor(plot_dt$go_id, levels=plot_dt$go_id)
# # plot_dt$go_id_lab <- factor(plot_dt$go_id_lab, levels=plot_dt$go_id_lab)
# # 
# # 
# # 
