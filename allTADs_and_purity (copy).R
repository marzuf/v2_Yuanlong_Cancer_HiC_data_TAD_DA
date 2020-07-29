########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript allTADs_and_purity.R\n"))

script_name <- "allTADs_and_purity.R"


####### !!!!!!!!!!!  FOR THE MOMENT I DO NOT RESTRICT 
# the datasets in which I look at expression (I do not ensure that is a DS where the conserved region is signif. DA)
# rationale: I want to look if this set of genes is related to purity, irrespective of DA


suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_cols[all_cols=="green"] <- "darkgreen"

do_plot <- function(my_x, my_y, ...) {
  plot(x=my_x,
       y=my_y,
       pch=16,
       cex=0.7,
       cex.axis=1.2,
       cex.main=1.2,
       cex.lab=1.2,
       ...)
  addCorr(my_x, my_y, bty="n")
}

# Rscript allTADs_and_purity_corrected.R

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

buildTable <- T

fontFamily <- "Hershey"

myHeight <- 400 
myWidth <- 400
plotType <- "png"
plotCex <- 1.4
mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

myHeightGG <- 6
myWidthGG <- 9

# plotType <- "svg"
# myHeightGG <- 7
# myWidthGG <- 11
# plotCex <- 1.2

corMet <- "pearson"

transfExpr <- "log10"
logOff <- 0.01


# for plotting
signifThresh <- 0.01
highCorrThresh <- 0.6

# to quickly retrieve tad-level stat.
all_result_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
} else{
  purity_ds <- ""
}

script0_name <- "0_prepGeneData"


if(purity_ds == "") {
  file_suffix <- ""
  purity_file <- file.path("tcga_purity_aran2015.csv")
  purity_dt <- read.delim(purity_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  purity_metrics <- c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")
#  pm <- purity_metrics[1]
#  purity_plot_name <- "aran"
  pm <- "CPE"
  purity_plot_name <- "aran - CPE"
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


  pm <- "CPE"
  purity_plot_name <- "aran - CPE"

  pm <- "ESTIMATE"
  purity_plot_name <- "aran"

stopifnot("Sample.ID" %in% purity_dt)

purity_dt$Sample.ID <- sapply(purity_dt$Sample.ID, function(x) substr(x, 1, 15))


outFolder <- file.path("ALLTADS_AND_PURITY_CORRECTED", purity_ds, pm, transfExpr)
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

mainFolder <- "."

# all_ds=all_ds[1]

ex_hicds <- "ENCSR489OCU_NCI-H460_40kb"
ex_exprds <- "TCGAluad_norm_luad"
ex_TAD <- "chr11_TAD390"

# all_ds = file.path(ex_hicds, ex_exprds)

if(buildTable){
  
  all_ds_corrPurity_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
    
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
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    
    geneList <- get(load(file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")))
    stopifnot(geneList %in% g2t_dt$entrezID)
    g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
    g2t_dt$region <- as.character(g2t_dt$region)

    
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
    
    stopifnot(names(geneList) %in% rownames(fpkm_dt))
    
    ### -> do the same as for the GIMAPs but for all TADs !
    tad = unique(g2t_dt$region)[1]
    all_tads_dt <- foreach(tad = unique(g2t_dt$region), .combine='rbind') %dopar% {
      
      tad_entrez <- g2t_dt$entrezID[g2t_dt$region == tad]
      stopifnot(tad_entrez %in% geneList)

      tad_entrez_fpkm <- names(geneList)[geneList %in% tad_entrez]
      stopifnot(tad_entrez_fpkm %in% rownames(fpkm_dt))
      tad_fpkm_dt <- fpkm_dt[tad_entrez_fpkm,]
      stopifnot(nrow(tad_fpkm_dt) == length(tad_entrez_fpkm))
      
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
      
      stopifnot(all_samps %in% colnames(tad_fpkm_dt))
      
      tad_dt <- data.frame(t(tad_fpkm_dt[,all_samps]), check.names=FALSE)
      
      stopifnot(is.numeric(unlist(tad_dt)))
      
      if(!is.null(transfExpr)) {
        if(grepl("log", transfExpr)) {
          tad_dt_2 <- do.call(transfExpr, list(tad_dt+logOff))
          stopifnot(dim(tad_dt_2)==dim(tad_dt))
          stopifnot(rownames(tad_dt_2) == rownames(tad_dt))
          stopifnot(colnames(tad_dt_2) == colnames(tad_dt))
          tad_dt <- tad_dt_2
        } else {stop("unnknown\n")}
        labTransf <- transfExpr
      } else {
        labTransf <- ""
      }
      
      stopifnot(colnames(tad_dt) %in% names(pipeline_geneList))
      stopifnot(!duplicated(pipeline_geneList))
      # reback to gff dt entrezID -> needed to add column if missing gmaps
      colnames(tad_dt) <- pipeline_geneList[colnames(tad_dt)]
      stopifnot(colnames(tad_dt) %in% pipeline_geneList)
      
      stopifnot(setequal(colnames(tad_dt),  tad_entrez))
      tad_dt$sampID <- rownames(tad_dt)
      rownames(tad_dt) <- NULL
      
      purity_dt <- data.frame(
        dataset=ds,
        sampID = names(purity_values),
        purity = purity_values,
        stringsAsFactors = FALSE
      )
      purity_expr_dt <- merge(purity_dt, tad_dt, by="sampID", all=TRUE)
      stopifnot(ncol(purity_expr_dt) == length(tad_entrez) + 3)
      purity_expr_dt <- purity_expr_dt[,c("dataset", "sampID", "purity", tad_entrez)]
      
      if(hicds==ex_hicds & exprds==ex_exprds & tad == ex_TAD) {
        for(i_g in tad_entrez) {
          outFile <- file.path(outFolder, paste0(ex_hicds, "_", ex_exprds, "_", ex_TAD, "_", i_g, "_expr_corr.", plotType))
          do.call(plotType, list(outFile, height=myHeight, width=myWidth))
          do_plot(my_x=purity_expr_dt[,"purity"],
                  my_y=purity_expr_dt[,paste0(i_g)],
                  main=paste0(hicds, " - ", exprds),
                  xlab="purity", ylab=paste0(i_g, " expr. ", labTransf))
          mtext(side=3, text=paste0(tad))
          foo <- dev.off()
          cat(paste0("... written: ", outFile, "\n"))
        }
      }
      # stop("ok \n")
      # correlation each column with purity 
      
      purity_expr_dt <- na.omit(purity_expr_dt)
      
      stopifnot(!duplicated(purity_expr_dt$sampID))
      
      if(nrow(purity_expr_dt) == 0) return(NULL)
      
      
      all_purityCors <- apply(purity_expr_dt[,tad_entrez],2, function(col) cor(col, purity_expr_dt$purity, method=corMet))
      
      stopifnot(!is.na(all_purityCors))
      
      
      data.frame(
        dataset=ds,
        nSampWithPurity=length(purity_expr_dt$sampID),
        region = tad,
        entrezID= names(all_purityCors),
        purityCorr = as.numeric(all_purityCors),
        stringsAsFactors = FALSE
      )
      
      
      
      
      
    }
    all_tads_dt
  }
  outFile <- file.path(outFolder,"all_ds_corrPurity_dt.Rdata" )
  save(all_ds_corrPurity_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))

} else {
  outFile <- file.path(outFolder,"all_ds_corrPurity_dt.Rdata" )
  all_ds_corrPurity_dt <- get(load(outFile))
  }
      

# do the TADs that are signif. correspond to those that are well correlated with purity ?
all_result_dt$regID <- file.path(all_result_dt$hicds, all_result_dt$exprds, all_result_dt$region)
all_ds_corrPurity_dt$regID <- file.path(all_ds_corrPurity_dt$dataset, all_ds_corrPurity_dt$region)
stopifnot(all_ds_corrPurity_dt$regID %in% all_result_dt$regID)

if(purity_ds == "EPIC") stopifnot(all_result_dt$regID %in% all_ds_corrPurity_dt$regID)

########### => Mean TAD-level
all_ds_meanCorrPurity_dt <- aggregate(purityCorr~regID, data=all_ds_corrPurity_dt, FUN=mean)
stopifnot(!duplicated(all_ds_meanCorrPurity_dt$regID))
all_dt <- merge(all_ds_meanCorrPurity_dt, all_result_dt, by="regID", all.x=TRUE, all.y=FALSE)
stopifnot(!duplicated(all_dt$regID))
all_dt$adjPvalComb_log10 <- -log10(all_dt$adjPvalComb)
all_dt$dotCols <- all_cols[all_cmps[paste0(all_dt$exprds)]]
stopifnot(!is.na(all_dt$dotCols))

plotTit <- paste0("TAD-level mean corr.")

myx_lab <- paste0(transfExpr, " expr. and purity correlation (meanTAD)")

all_vars <- c("adjPvalComb_log10", "meanCorr", "ratioDown")

plot_var =all_vars[1]

for(plot_var in all_vars) {
  my_x <- all_dt$purityCorr
  my_y <- all_dt[,c(plot_var)]
  outFile <- file.path(outFolder, paste0(plot_var, "_exprPurityCorr_meanTAD.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_plot(my_x,my_y, xlab=myx_lab, 
          col = all_dt$dotCols,
          main=paste0("meanTAD - ", plot_var, " vs. purity corr."),
          ylab=paste0(plot_var))
  mtext(side=3, text = paste0(corMet, "'s corr.", " - ", purity_plot_name, " data"))
  legend("bottomleft",all_cols, legend=names(all_cols), bty="n", pch=16, col=all_cols)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

onlysignif_dt <- all_dt[all_dt$adjPvalComb <= signifThresh,]
for(plot_var in all_vars) {
  my_x <- onlysignif_dt$purityCorr
  my_y <- onlysignif_dt[,c(plot_var)]
  outFile <- file.path(outFolder, paste0(plot_var, "_exprPurityCorr_meanTAD_onlySignif.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_plot(my_x,my_y, xlab=myx_lab, 
          col = onlysignif_dt$dotCols,
          main=paste0("meanTAD - ", plot_var, " vs. purity corr."),
          ylab=paste0(plot_var))
  mtext(side=3, text = paste0(corMet, "'s corr. - adjPvalComb TAD <= ", signifThresh, " - ", purity_plot_name, " data"))
  legend("bottomleft",all_cols, legend=names(all_cols), bty="n", pch=16, col=all_cols)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

onlyCorr_dt <- all_dt[abs(all_dt$purityCorr) >= highCorrThresh,]
for(plot_var in all_vars) {
  my_x <- onlyCorr_dt$purityCorr
  my_y <- onlyCorr_dt[,c(plot_var)]
  outFile <- file.path(outFolder, paste0(plot_var, "_exprPurityCorr_meanTAD_highCorr.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  do_plot(my_x,my_y, xlab=myx_lab, 
          col = onlyCorr_dt$dotCols,
          main=paste0("meanTAD - ", plot_var, " vs. purity corr."),
          ylab=paste0(plot_var))
  mtext(side=3, text = paste0(corMet, "'s corr. - corr. expr-purity abs. >= ", highCorrThresh, " - ", purity_plot_name, " data" ))
  legend("bottomleft",all_cols, legend=names(all_cols), bty="n", pch=16, col=all_cols)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

myx_lab <- paste0(transfExpr, " expr. and purity correlation (mean TAD)")

outFile <- file.path(outFolder, paste0("exprPurityCorr_signif_notSignif_density_meanTAD.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.5))
plot_multiDens(list(
  "signif. TADs" =onlysignif_dt$purity,
  "not signif. TADs" = all_dt$purityCorr[all_dt$adjPvalComb > signifThresh]),legPos = "topleft", my_xlab = myx_lab,
  plotTit =paste0(corMet, "'s corr. expr.-purity (mean TAD)"))
mtext(text=paste0("adj. p-val. signif. thresh = ", signifThresh, " - ", purity_plot_name, " data"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


onlyCorr_dt <- all_dt[abs(all_dt$purityCorr) >= highCorrThresh,]
myx_lab <- paste0("adjPvalComb_log10 (meanTAD)")
outFile <- file.path(outFolder, paste0("adjPvalCombLog10_highCorr_lowCorr_density_meanTAD.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.5))
plot_multiDens(list(
  "high corr. TADs" =onlyCorr_dt$adjPvalComb_log10,
  "low corr. TADs" = all_dt$adjPvalComb_log10[all_dt$purityCorr < highCorrThresh]),legPos = "topright",
  plotTit =paste0(corMet, "'s corr. expr.-purity (mean TAD)"), my_xlab = myx_lab)
mtext(text=paste0("abs. high corr. thresh = ", highCorrThresh, " - ", purity_plot_name, " data"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


# dataset boxplot
all_ds_meanCorrPurity_dt$dataset <- dirname(all_ds_meanCorrPurity_dt$regID)
aggds_meanCorrPurit_dt <- aggregate(purityCorr~dataset, data=all_ds_meanCorrPurity_dt, FUN=mean)
aggds_meanCorrPurit_dt <- aggds_meanCorrPurit_dt[order(aggds_meanCorrPurit_dt$purityCorr),]
all_ds_meanCorrPurity_dt$dataset <- factor(all_ds_meanCorrPurity_dt$dataset, levels=aggds_meanCorrPurit_dt$dataset)
myx_lab <- paste0(transfExpr, " expr. and purity correlation (meanTAD)")

plotTit <- paste0("Expr. ", transfExpr, " purity corr. - by DS")
subTit <- paste0(corMet, "'s corr. - ", purity_plot_name, " data")
p <- ggplot(all_ds_meanCorrPurity_dt, aes(x=dataset, y=purityCorr))+
  geom_boxplot()+
  ggtitle(plotTit, subtitle = subTit)+
  # geom_boxplot(notch = TRUE, outlier.shape=NA)+
  # geom_jitter()+
  scale_y_continuous(breaks=scales::pretty_breaks(n = 10))+
  scale_x_discrete(labels=function(x) paste0(dirname(x), "\n", basename(x)))+
  labs(x="",
       y =paste0(myx_lab))+ 
  theme(
    plot.margin=margin(t = 0, r = 0, b = 10, l = 0, unit = "pt"),
    text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.line = element_line(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_text(size=6, hjust=1, vjust=0.5, angle=90, color=all_cols[all_cmps[basename(aggds_meanCorrPurit_dt$dataset)]]),
    # axis.text.x = element_blank(),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold"),
    legend.text = element_text(size=12)
  ) 

outFile <- file.path(outFolder, paste0("all_purity_meanTAD_byDS_boxplot.", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



########### => all gene-level
all_dt <- merge(all_ds_corrPurity_dt, all_result_dt, by="regID", all.x=TRUE, all.y=FALSE)
all_dt$adjPvalComb_log10 <- -log10(all_dt$adjPvalComb)
stopifnot(!duplicated(file.path(all_dt$regID, all_dt$entrezID)))

plotTit <- paste0("Gene-level")
myx_lab <- paste0(transfExpr, " expr. and purity correlation (gene-level)")

myx_lab <- paste0(transfExpr, " expr. and purity correlation (gene-level)")

onlysignif_dt <- all_dt[all_dt$adjPvalComb <= signifThresh,]

outFile <- file.path(outFolder, paste0("exprPurityCorr_signif_notSignif_density_geneLevel.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.5))
plot_multiDens(list(
  "signif. TADs" =onlysignif_dt$purity,
  "not signif. TADs" = all_dt$purity[all_dt$adjPvalComb > signifThresh]),legPos = "topleft",
  plotTit =paste0(corMet, "'s corr. expr.-purity (gene-level)"), my_xlab = myx_lab)
mtext(text=paste0("adj. p-val. signif. thresh = ", signifThresh," - ", purity_plot_name, " data"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

onlyCorr_dt <- all_dt[abs(all_dt$purityCor)  >= highCorrThresh,]

myx_lab <- paste0("adjPvalComb_log10 (gene-level)")

outFile <- file.path(outFolder, paste0("adjPvalCombLog10_highCorr_lowCorr_density_geneLevel.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.5))
plot_multiDens(list(
  "high corr. TADs" =onlyCorr_dt$adjPvalComb_log10,
  "low corr. TADs" = all_dt$adjPvalComb_log10[all_dt$purityCorr < highCorrThresh]),legPos = "topright",
  plotTit =paste0(corMet, "'s corr. expr.-purity (gene-level)"), my_xlab = myx_lab)
mtext(text=paste0("abs. high corr. thresh = ", highCorrThresh, " - ", purity_plot_name, " data"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
entrez2symb <- setNames(gff_dt$symbol,gff_dt$entrezID)

# meanCorr by gene

mean_agg <- aggregate(purityCorr ~ entrezID, data=all_ds_corrPurity_dt, FUN=mean)
mean_agg <- mean_agg[order(abs(mean_agg$purityCorr), decreasing=T),]
mean_agg$symbol <- entrez2symb[mean_agg$entrezID]
head(mean_agg)


onlysignif_dt <- all_dt[all_dt$adjPvalComb <= signifThresh,]
onlysignif_dt$purityCorr_rd <- round(onlysignif_dt$purityCorr, 2)
onlysignif_dt$adjPvalComb_log10_rd <- round(onlysignif_dt$adjPvalComb_log10, 2)
onlysignif_dt <- onlysignif_dt[order(onlysignif_dt$purityCorr), c("regID", "region_genes", "purityCorr_rd","adjPvalComb_log10_rd")]
rownames(onlysignif_dt) <- NULL
head(onlysignif_dt,20)

######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


