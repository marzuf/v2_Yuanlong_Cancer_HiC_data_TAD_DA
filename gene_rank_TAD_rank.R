options(scipen=100)

SSHFS=F

setDir <- "/media/electron"
setDir <- ""

# Rscript gene_rank_TAD_rank.R

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "gene_rank_TAD_rank.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(RColorBrewer)
require(reshape2)
require(ggplot2)
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

outFolder <- file.path("GENE_RANK_TAD_RANK")
dir.create(outFolder, recursive = TRUE)
logFile <- file.path(outFolder, "gene_rank_TAD_rank_logFile.txt")
if(buildTable) file.remove(logFile)

printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}

tad_signif_col <- "tad_adjCombPval"
gene_signif_col <- "adj.P.Val"


tad_pval_thresh <- 0.01
gene_pval_thresh <- 0.05
tad_pval_thresh_log10 <- -log10(tad_pval_thresh)
gene_pval_thresh_log10 <- -log10(gene_pval_thresh)
file_suffix <- paste0("tad_pval_thresh", tad_pval_thresh, "_gene_pval_thresh", gene_pval_thresh)


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
      tad_dt$tad_rank <- rank(tad_dt$tad_adjCombPval, ties=tieMeth)
      
      
      tad_gene_dt <- merge(exprds_g2t_dt[,c("entrezID", "region")], tad_dt, by="region", all.x=TRUE, all.y=FALSE)
      
      out_dt <- merge(tad_gene_dt, topTable_DT[,c("entrezID", "logFC", "adj.P.Val")], by="entrezID", all.x=TRUE, all.y=FALSE )
      out_dt <- unique(out_dt)
      stopifnot(!duplicated(out_dt$entrezID))
      stopifnot(exprds_g2t_dt$entrezID %in% out_dt$entrezID)
      
      out_dt$gene_rank <- rank(out_dt$adj.P.Val, ties=tieMeth)
      # out_dt$tad_rank <- rank(out_dt$tad_adjCombPval, ties=tieMeth)
      
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



###################################################################################### start plotting

all_gene_tad_signif_dt[, paste0(gene_signif_col, "_log10")] <- -log10(all_gene_tad_signif_dt[, paste0(gene_signif_col)])
all_gene_tad_signif_dt[, paste0(tad_signif_col, "_log10")] <- -log10(all_gene_tad_signif_dt[, paste0(tad_signif_col)])

all_y <- c("logFC", paste0(gene_signif_col), paste0(gene_signif_col, "_log10"), "gene_rank")
all_x <- c(paste0(tad_signif_col), paste0(tad_signif_col, "_log10"), "tad_rank")

stopifnot(all_y %in% colnames(all_gene_tad_signif_dt))
stopifnot(all_x %in% colnames(all_gene_tad_signif_dt))

nDS <- length(unique(file.path(all_gene_tad_signif_dt$hicds, all_gene_tad_signif_dt$exprds)))

xvar=all_x[1]
yvar=all_y[1]
for(xvar in all_x) {
  for(yvar in all_y) {
    
    myx <- all_gene_tad_signif_dt[,paste0(xvar)]
    myy <- all_gene_tad_signif_dt[,paste0(yvar)]
    
    outFile <- file.path(outFolder, paste0("allDS_", yvar, "_vs_", xvar, "_densplot.", plotType))    
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = myx,
      y = myy,
      xlab = paste0(xvar),
      ylab = paste0(yvar),
      main = paste0(yvar, " vs. ", xvar),
      cex.lab = axisCex,
      cex.axis = axisCex
    )
    mtext(side=3, text=paste0("all DS; n = ", nDS))
    addCorr(x=myx, y=myy, legPos = "topleft", bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))  
  }  
  

}


xvar <- paste0(tad_signif_col, "_log10")
yvar <- paste0(gene_signif_col, "_log10")

myx <- all_gene_tad_signif_dt[,paste0(xvar)]
myy <- all_gene_tad_signif_dt[,paste0(yvar)]

all_gene_tad_signif_dt$dotCols <- ifelse( all_gene_tad_signif_dt[,paste0(xvar)] >= tad_pval_thresh & all_gene_tad_signif_dt[,paste0(yvar)] >= gene_pval_thresh, "red", 
                                          ifelse( all_gene_tad_signif_dt[,paste0(xvar)] >= tad_pval_thresh & all_gene_tad_signif_dt[,paste0(yvar)] < gene_pval_thresh, "green", "grey"))

outFile <- file.path(outFolder, paste0("allDS_", yvar, "_vs_", xvar, "_withThresh_densplot.", plotType))    
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = myx,
  y = myy,
  xlab = paste0(xvar),
  ylab = paste0(yvar),
  main = paste0(yvar, " vs. ", xvar),
  col = all_gene_tad_signif_dt$dotCols,
  pch=16,
  cex=0.7,
  cex.lab = axisCex,
  cex.axis = axisCex
)
abline(h=gene_pval_thresh, lty=2, col="grey")
abline(v=tad_pval_thresh, lty=2, col="grey")
mtext(side=3, text=paste0("all DS; n = ", nDS))
addCorr(x=myx, y=myy, legPos = "topleft", bty="n")
legend("topright", 
       legend=c(paste0("n=", sum(all_gene_tad_signif_dt$dotCols == "red")), paste0("n=", sum(all_gene_tad_signif_dt$dotCols == "green")), paste0("n=", sum(all_gene_tad_signif_dt$dotCols == "grey"))),
       text.col = c("red", "green", "grey"),
       bty="n"
       )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))  

  

# limmaMissed_dt_nByTAD <- aggregate(entrezID ~ hicds+exprds+region, data=limmaMissed_dt, FUN=length)

###################################################################################### ZOOM PLOT SEPARATELY FOR EACH DATASET
limmaMissed_dt <- all_gene_tad_signif_dt[all_gene_tad_signif_dt[,paste0(tad_signif_col)] <= tad_pval_thresh & all_gene_tad_signif_dt[,paste0(gene_signif_col)] > gene_pval_thresh,]

limmaMissed_dt$gene_symbol <- entrez2symb[paste0(limmaMissed_dt$entrezID)]
stopifnot(!is.na(limmaMissed_dt$gene_symbol))

all_y <- c("logFC", paste0(gene_signif_col), paste0(gene_signif_col, "_log10"), "gene_rank")
all_x <- c(paste0(tad_signif_col), paste0(tad_signif_col, "_log10"), "tad_rank")

stopifnot(all_y %in% colnames(limmaMissed_dt))
stopifnot(all_x %in% colnames(limmaMissed_dt))

nDS <- length(unique(file.path(limmaMissed_dt$hicds, limmaMissed_dt$exprds)))

xvar=all_x[1]
yvar=all_y[1]
for(xvar in all_x) {
  for(yvar in all_y) {
    myx <- limmaMissed_dt[,paste0(xvar)]
    myy <- limmaMissed_dt[,paste0(yvar)]
    outFile <- file.path(outFolder, paste0("allDS_", yvar, "_vs_", xvar, "_subMissed_densplot.", plotType))    
    do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    densplot(
      x = myx,
      y = myy,
      xlab = paste0(xvar),
      ylab = paste0(yvar),
      main = paste0(yvar, " vs. ", xvar),
      cex.lab = axisCex,
      cex.axis = axisCex
    )
    mtext(side=3, text=paste0("all DS; n = ", nDS))
    addCorr(x=myx, y=myy, legPos = "topleft", bty="n")
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))  
  }  
}

###################################################################################### ZOOM PLOT SEPARATELY FOR EACH DATASET - signif. tad vs. signif. gene

yvar <- paste0(gene_signif_col, "_log10")
xvar <- paste0(tad_signif_col, "_log10")
nDS <- length(unique(file.path(limmaMissed_dt$hicds, limmaMissed_dt$exprds)))
myx <- limmaMissed_dt[,paste0(xvar)]
myy <- limmaMissed_dt[,paste0(yvar)]
limmaMissed_dt$ds <- paste0(limmaMissed_dt$hicds,"_",limmaMissed_dt$exprds)
outFile <- file.path(outFolder, paste0("allDS_", yvar, "_vs_", xvar, "_subMissed_label_scatter.", plotType))    
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = myx,
  y = myy,
  xlab = paste0(xvar),
  ylab = paste0(yvar),
  main = paste0(yvar, " vs. ", xvar),
  sub=paste0("(n=",length(myx), ")"),
  col = "grey",
  pch = 16,
  cex=0.6,
  cex.lab = axisCex,
  cex.axis = axisCex
)
text(x=myx, y=myy, labels = limmaMissed_dt$gene_symbol, cex=0.6)
mtext(side=3, text=paste0("all DS; n = ", nDS))
addCorr(x=myx, y=myy, legPos = "topleft", bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))  


yvar <- paste0(gene_signif_col, "_log10")
xvar <- paste0(tad_signif_col, "_log10")
ds = unique(limmaMissed_dt$ds)[1]
for(ds in unique(limmaMissed_dt$ds)) {
  sub_limmaMissed_dt <- limmaMissed_dt[limmaMissed_dt$ds == ds,]
  nDS <- length(unique(file.path(sub_limmaMissed_dt$hicds, sub_limmaMissed_dt$exprds)))
  stopifnot(nDS==1)
  myx <- sub_limmaMissed_dt[,paste0(xvar)]
  myy <- sub_limmaMissed_dt[,paste0(yvar)]
  
  outFile <- file.path(outFolder, paste0(ds, "_", yvar, "_vs_", xvar, "_subMissed_label_scatter.", plotType))    
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = myx,
    y = myy,
    xlab = paste0(xvar),
    ylab = paste0(yvar),
    main = paste0(yvar, " vs. ", xvar),
    sub=paste0("(n=",length(myx), ")"),
    col = "grey",
    pch = 16,
    cex=0.6,
    cex.lab = axisCex,
    cex.axis = axisCex
  )
  text(x=myx, y=myy, labels = sub_limmaMissed_dt$gene_symbol, cex=0.6)
  mtext(side=3, text=paste0(ds))
  # addCorr(x=myx, y=myy, legPos = "topleft", bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
}

###################################################################################### ZOOM PLOT SEPARATELY FOR EACH DATASET - signif. tad vs. signif. gene - SAME BUT ONLY IF SIGN gene FC == sign TAD FC

limmaMissed_dt <- all_gene_tad_signif_dt[all_gene_tad_signif_dt[,paste0(tad_signif_col)] <= tad_pval_thresh & all_gene_tad_signif_dt[,paste0(gene_signif_col)] > gene_pval_thresh,]

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

fc_limmaMissed_dt <- merge(limmaMissed_dt, final_dt[,c("hicds", "exprds","region", "meanLogFC" )], by=c("hicds", "exprds",  "region"), all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(fc_limmaMissed_dt$meanLogFC))
stopifnot(!is.na(fc_limmaMissed_dt$logFC))

fc_limmaMissed_dt <- fc_limmaMissed_dt[sign(fc_limmaMissed_dt$meanLogFC) == sign(fc_limmaMissed_dt$logFC),]

fc_limmaMissed_dt$gene_symbol <- entrez2symb[paste0(fc_limmaMissed_dt$entrezID)]
stopifnot(!is.na(fc_limmaMissed_dt$gene_symbol))


yvar <- paste0(gene_signif_col, "_log10")
xvar <- paste0(tad_signif_col, "_log10")
nDS <- length(unique(file.path(fc_limmaMissed_dt$hicds, fc_limmaMissed_dt$exprds)))
myx <- fc_limmaMissed_dt[,paste0(xvar)]
myy <- fc_limmaMissed_dt[,paste0(yvar)]
fc_limmaMissed_dt$ds <- paste0(fc_limmaMissed_dt$hicds,"_",fc_limmaMissed_dt$exprds)
outFile <- file.path(outFolder, paste0("allDS_", yvar, "_vs_", xvar, "_subMissed_label_scatter_sameFCsign.", plotType))    
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x = myx,
  y = myy,
  xlab = paste0(xvar),
  ylab = paste0(yvar),
  main = paste0(yvar, " vs. ", xvar),
  sub=paste0("(n=",length(myx), ") - gene FC sign == TAD FC sign"),
  # sub = paste0("! only gene FC sign == TAD FC sign !"),
  col = "grey",
  pch = 16,
  cex=0.6,
  cex.lab = axisCex,
  cex.axis = axisCex
)
text(x=myx, y=myy, labels = fc_limmaMissed_dt$gene_symbol, cex=0.6)
mtext(side=3, text=paste0("all DS; n = ", nDS))
addCorr(x=myx, y=myy, legPos = "topleft", bty="n")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))  


yvar <- paste0(gene_signif_col, "_log10")
xvar <- paste0(tad_signif_col, "_log10")
ds = unique(fc_limmaMissed_dt$ds)[1]
for(ds in unique(fc_limmaMissed_dt$ds)) {
  sub_limmaMissed_dt <- fc_limmaMissed_dt[fc_limmaMissed_dt$ds == ds,]
  nDS <- length(unique(file.path(sub_limmaMissed_dt$hicds, sub_limmaMissed_dt$exprds)))
  stopifnot(nDS==1)
  myx <- sub_limmaMissed_dt[,paste0(xvar)]
  myy <- sub_limmaMissed_dt[,paste0(yvar)]
  
  outFile <- file.path(outFolder, paste0(ds, "_", yvar, "_vs_", xvar, "_subMissed_label_scatter_sameFCsign.", plotType))    
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  plot(
    x = myx,
    y = myy,
    xlab = paste0(xvar),
    ylab = paste0(yvar),
    main = paste0(yvar, " vs. ", xvar),
    sub=paste0("(n=",length(myx), ") - gene FC sign == TAD FC sign"),
    # sub = paste0("! only gene FC sign == TAD FC sign !"),
    col = "grey",
    pch = 16,
    cex=0.6,
    cex.lab = axisCex,
    cex.axis = axisCex
  )
  text(x=myx, y=myy, labels = sub_limmaMissed_dt$gene_symbol, cex=0.6)
  mtext(side=3, text=paste0(ds))
  # addCorr(x=myx, y=myy, legPos = "topleft", bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))  
}






###################################################################################### TEXT TABLE

limmaMissed_dt <- all_gene_tad_signif_dt[all_gene_tad_signif_dt[,paste0(tad_signif_col)] <= tad_pval_thresh & all_gene_tad_signif_dt[,paste0(gene_signif_col)] > gene_pval_thresh,]
limmaMissed_dt$dataset <- file.path(limmaMissed_dt$hicds, limmaMissed_dt$exprds)
limmaMissed_dt$cmpType <- all_cmps[paste0(limmaMissed_dt$exprds)]
stopifnot(!is.na(limmaMissed_dt$cmpType))
stopifnot(limmaMissed_dt$entrezID %in% names(entrez2symb))
all_entrez <- unique(as.character(limmaMissed_dt$entrezID))

ds_collapse_dt <- aggregate(dataset ~ entrezID, data=limmaMissed_dt, FUN=function(x) paste0(sort(x), collapse=","))
entrez2ds <- setNames(ds_collapse_dt$dataset, ds_collapse_dt$entrezID)

entrez2nmiss <- setNames(as.numeric(table(limmaMissed_dt$entrezID)), names(table(limmaMissed_dt$entrezID)) )

stopifnot(setequal(all_entrez, names(entrez2nmiss)))
stopifnot(setequal(all_entrez, names(entrez2ds)))

nLimmaMissed_dt <- data.frame(
  entrezID = all_entrez,
  gene_symbol = entrez2symb[paste0(all_entrez)],
  nMissed = entrez2nmiss[paste0(all_entrez)],
  dsMissed = entrez2ds[paste0(all_entrez)],
  stringsAsFactors = FALSE
)

stopifnot(!is.na(nLimmaMissed_dt))

nLimmaMissed_dt <- nLimmaMissed_dt[order(nLimmaMissed_dt$nMissed, decreasing=TRUE),]

outFile <- file.path(outFolder, paste0(file_suffix, "_nLimaMissed_dt.txt"))
write.table(nLimmaMissed_dt, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, append=F, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


all_cmpTypes <- unique(as.character(all_cmps))
ct = all_cmpTypes[1]
for(ct in all_cmpTypes) {
  
  ct_limmaMissed_dt <- limmaMissed_dt[limmaMissed_dt$cmpType == ct,]

  ct_all_entrez <- unique(as.character(ct_limmaMissed_dt$entrezID))
  
  ct_ds_collapse_dt <- aggregate(dataset ~ entrezID, data=ct_limmaMissed_dt, FUN=function(x) paste0(sort(x), collapse=","))
  ct_entrez2ds <- setNames(ct_ds_collapse_dt$dataset, ct_ds_collapse_dt$entrezID)
  
  ct_entrez2nmiss <- setNames(as.numeric(table(ct_limmaMissed_dt$entrezID)), names(table(ct_limmaMissed_dt$entrezID)) )
  
  stopifnot(setequal(ct_all_entrez, names(ct_entrez2nmiss)))
  stopifnot(setequal(ct_all_entrez, names(ct_entrez2ds)))
  
  ct_nLimmaMissed_dt <- data.frame(
    entrezID = ct_all_entrez,
    gene_symbol = entrez2symb[paste0(ct_all_entrez)],
    nMissed = ct_entrez2nmiss[paste0(ct_all_entrez)],
    dsMissed = ct_entrez2ds[paste0(ct_all_entrez)],
    stringsAsFactors = FALSE
  )
  stopifnot(!is.na(ct_nLimmaMissed_dt))
  
  ct_nLimmaMissed_dt <- ct_nLimmaMissed_dt[order(ct_nLimmaMissed_dt$nMissed, decreasing=TRUE),]

  outFile <- file.path(outFolder, paste0(file_suffix, "_nLimaMissed_dt_", ct, ".txt"))
  write.table(ct_nLimmaMissed_dt, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, append=F, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  

  
    
  
}


###################################################################################### TEXT TABLE - same sign only

limmaMissed_dt <- all_gene_tad_signif_dt[all_gene_tad_signif_dt[,paste0(tad_signif_col)] <= tad_pval_thresh & all_gene_tad_signif_dt[,paste0(gene_signif_col)] > gene_pval_thresh,]


final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

limmaMissed_dt <- merge(limmaMissed_dt, final_dt[,c("hicds", "exprds","region", "meanLogFC" )], by=c("hicds", "exprds",  "region"), all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(limmaMissed_dt$meanLogFC))
stopifnot(!is.na(limmaMissed_dt$logFC))

limmaMissed_dt <- limmaMissed_dt[sign(limmaMissed_dt$meanLogFC) == sign(limmaMissed_dt$logFC),]

limmaMissed_dt_nByTAD <- aggregate(entrezID ~ hicds+exprds+region, data=limmaMissed_dt, FUN=length)
head(limmaMissed_dt_nByTAD)

limmaMissed_dt_nByTAD$dataset <- paste0(limmaMissed_dt_nByTAD$hicds, "\n", limmaMissed_dt_nByTAD$exprds)
tmpdt <- aggregate(entrezID~dataset, FUN=mean, data=limmaMissed_dt_nByTAD)
tmpdt <- tmpdt[order(tmpdt$entrezID, decreasing=T),]
ds_lev <- tmpdt$dataset

mycols <- all_cols[all_cmps[gsub(".+\n(.+)", "\\1", ds_lev)]]

limmaMissed_dt_nByTAD$dataset <- factor(limmaMissed_dt_nByTAD$dataset, levels=ds_lev)
p_var <-  ggplot(limmaMissed_dt_nByTAD, aes(x = dataset, y = entrezID)) + 
  geom_boxplot() +
  coord_cartesian(expand = FALSE) +
  ggtitle("# limma missed genes by TAD", subtitle = paste0())+
  scale_x_discrete(name="")+
  # labs(fill="")+
  scale_y_continuous(name=paste0("# missed genes/TAD"),
                     breaks = scales::pretty_breaks(n = 10))+
  # scale_fill_manual(values=c(fcc_auc=fcc_col, coexpr_auc=coexpr_col), labels=c("FCC", "coexpr."))+
  # geom_hline(yintercept=1, linetype=2)+
  theme( # Increase size of axis lines
    # strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    # strip.text.x = element_text(size = 10),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    # legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("all_ds_nbr_limma_missed_genes_byTAD_boxplot.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

stop("--ok\n")



limmaMissed_dt$gene_symbol <- entrez2symb[paste0(limmaMissed_dt$entrezID)]
stopifnot(!is.na(limmaMissed_dt$gene_symbol))


limmaMissed_dt$dataset <- file.path(limmaMissed_dt$hicds, limmaMissed_dt$exprds)
limmaMissed_dt$cmpType <- all_cmps[paste0(limmaMissed_dt$exprds)]
stopifnot(!is.na(limmaMissed_dt$cmpType))
stopifnot(limmaMissed_dt$entrezID %in% names(entrez2symb))
all_entrez <- unique(as.character(limmaMissed_dt$entrezID))

ds_collapse_dt <- aggregate(dataset ~ entrezID, data=limmaMissed_dt, FUN=function(x) paste0(sort(x), collapse=","))
entrez2ds <- setNames(ds_collapse_dt$dataset, ds_collapse_dt$entrezID)

entrez2nmiss <- setNames(as.numeric(table(limmaMissed_dt$entrezID)), names(table(limmaMissed_dt$entrezID)) )

stopifnot(setequal(all_entrez, names(entrez2nmiss)))
stopifnot(setequal(all_entrez, names(entrez2ds)))

nLimmaMissed_dt <- data.frame(
  entrezID = all_entrez,
  gene_symbol = entrez2symb[paste0(all_entrez)],
  nMissed = entrez2nmiss[paste0(all_entrez)],
  dsMissed = entrez2ds[paste0(all_entrez)],
  stringsAsFactors = FALSE
)

stopifnot(!is.na(nLimmaMissed_dt))

nLimmaMissed_dt <- nLimmaMissed_dt[order(nLimmaMissed_dt$nMissed, decreasing=TRUE),]

outFile <- file.path(outFolder, paste0(file_suffix, "_nLimaMissed_dt_sameFCsign.txt"))
write.table(nLimmaMissed_dt, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, append=F, file=outFile)
cat(paste0("... written: ", outFile, "\n"))


all_cmpTypes <- unique(as.character(all_cmps))
ct = all_cmpTypes[1]
for(ct in all_cmpTypes) {
  
  ct_limmaMissed_dt <- limmaMissed_dt[limmaMissed_dt$cmpType == ct,]
  
  ct_all_entrez <- unique(as.character(ct_limmaMissed_dt$entrezID))
  
  ct_ds_collapse_dt <- aggregate(dataset ~ entrezID, data=ct_limmaMissed_dt, FUN=function(x) paste0(sort(x), collapse=","))
  ct_entrez2ds <- setNames(ct_ds_collapse_dt$dataset, ct_ds_collapse_dt$entrezID)
  
  ct_entrez2nmiss <- setNames(as.numeric(table(ct_limmaMissed_dt$entrezID)), names(table(ct_limmaMissed_dt$entrezID)) )
  
  stopifnot(setequal(ct_all_entrez, names(ct_entrez2nmiss)))
  stopifnot(setequal(ct_all_entrez, names(ct_entrez2ds)))
  
  ct_nLimmaMissed_dt <- data.frame(
    entrezID = ct_all_entrez,
    gene_symbol = entrez2symb[paste0(ct_all_entrez)],
    nMissed = ct_entrez2nmiss[paste0(ct_all_entrez)],
    dsMissed = ct_entrez2ds[paste0(ct_all_entrez)],
    stringsAsFactors = FALSE
  )
  stopifnot(!is.na(ct_nLimmaMissed_dt))
  
  ct_nLimmaMissed_dt <- ct_nLimmaMissed_dt[order(ct_nLimmaMissed_dt$nMissed, decreasing=TRUE),]
  
  outFile <- file.path(outFolder, paste0(file_suffix, "_nLimaMissed_dt_", ct, "_sameFCsign.txt"))
  write.table(ct_nLimmaMissed_dt, col.names=TRUE, row.names=FALSE, sep="\t", quote=F, append=F, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  
}




###################################################################################### boxplot missed

limmaMissed_dt <- all_gene_tad_signif_dt[all_gene_tad_signif_dt[,paste0(tad_signif_col)] <= tad_pval_thresh & all_gene_tad_signif_dt[,paste0(gene_signif_col)] > gene_pval_thresh,]

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

limmaMissed_dt <- merge(limmaMissed_dt, final_dt[,c("hicds", "exprds","region", "meanLogFC" )], by=c("hicds", "exprds",  "region"), all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(limmaMissed_dt$meanLogFC))
stopifnot(!is.na(limmaMissed_dt$logFC))

limmaMissed_dt$sameGeneTADsign <- sign(limmaMissed_dt$meanLogFC) == sign(limmaMissed_dt$logFC)

limmaMissed_dt$dataset <- paste0(limmaMissed_dt$hicds, "\n", limmaMissed_dt$exprds)

count_limmaMissed_dt <- aggregate(entrezID ~ dataset + sameGeneTADsign, data=limmaMissed_dt, FUN=length)

tmp <- count_limmaMissed_dt[count_limmaMissed_dt$sameGeneTADsign,]
tmp <- tmp[order(tmp$entrezID, decreasing = TRUE),]
ds_lev <- as.character(tmp$dataset)

count_limmaMissed_dt$dataset <- factor(count_limmaMissed_dt$dataset, levels=ds_lev)
stopifnot(!is.na(count_limmaMissed_dt$dataset))

mycols <- all_cols[all_cmps[gsub(".+\n(.+)", "\\1", ds_lev)]]

# ggbarplot(data=count_limmaMissed_dt, x = "dataset", y="entrezID", fill="sameGeneTADsign")
# 
# 
# boxplot(entrezID~dataset+sameGeneTADsign, data=count_limmaMissed_dt)

p_var <-  ggplot(count_limmaMissed_dt, aes(x = dataset, y = entrezID, fill = sameGeneTADsign)) + 
  geom_bar(position="dodge", stat="identity") +
  coord_cartesian(expand = FALSE) +
  ggtitle("# limma missed genes", subtitle = paste0())+
  scale_x_discrete(name="")+
  # labs(fill="")+
  scale_y_continuous(name=paste0("# genes"),
                     breaks = scales::pretty_breaks(n = 10))+
  # scale_fill_manual(values=c(fcc_auc=fcc_col, coexpr_auc=coexpr_col), labels=c("FCC", "coexpr."))+
  # geom_hline(yintercept=1, linetype=2)+
  theme( # Increase size of axis lines
    # strip.text = element_text(size = 12),
    # top, right, bottom and left
    # plot.margin = unit(c(1, 1, 4.5, 1), "lines"),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    # strip.text.x = element_text(size = 10),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5),
    axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=12),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    # legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0("all_ds_nbr_limma_missed_genes_barplot.", plotType))
ggsave(plot = p_var, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))








######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))



