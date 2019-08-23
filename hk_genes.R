options(scipen=100)

# Rscript hk_genes.R

startTime <- Sys.time()

script_name <- "hk_genes.R"

cat("> START ", script_name," \n")

setDir <- "/media/electron"
setDir <- ""

buildTable <- TRUE

require(foreach)
require(doMC)
registerDoMC(40)


source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

pipFolder <- file.path("PIPELINE","OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))

outFolder <- paste0("HK_GENES")
dir.create(outFolder, recursive=TRUE)

hkGenes_file <- file.path("..", "2_Yuanlong_Cancer_HiC_data_TAD_DA", "HK_genes.txt")
stopifnot(file.exists(hkGenes_file))
hkGenes_DT <- read.delim(hkGenes_file, header=FALSE, stringsAsFactors = FALSE, col.names=c("symbol", "id"))
hkGenes_DT$symbol <- gsub(" ", "", hkGenes_DT$symbol)
hkGenes_DT$id <- gsub(" ", "", hkGenes_DT$id)

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)

cat(paste0("... # of HK genes with positions available:\t", sum(hkGenes_DT$symbol %in% entrez2symb_dt$symbol), "/", length(hkGenes_DT$symbol),
           " (",round(sum(hkGenes_DT$symbol %in% entrez2symb_dt$symbol)/length(hkGenes_DT$symbol)*100, 2), "%)\n" ))

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))


if(buildTable) {

  hicds="GSE105381_HepG2_40kb"
  all_hk_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds="TCGAlihc_norm_lihc"  
    ds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %dopar% {
      cat("... start ", hicds, " - ", exprds, "\n")
      
      g2tFile <- file.path(dirname(dirname(pipFolder)), hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      
      geneList_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      geneList <- eval(parse(text = load(geneList_file)))
      
      stopifnot(geneList %in% g2t_DT$entrezID)
      
      g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
      stopifnot(g2t_DT$entrezID %in% entrez2symb_dt$entrezID)
      g2t_DT$symbol <- sapply(g2t_DT$entrezID, function(x) entrez2symb_dt$symbol[entrez2symb_dt$entrezID == x])
      g2t_DT$is_hk <- as.numeric(g2t_DT$symbol %in% hkGenes_DT$symbol)
      
      all_regs <- as.character(g2t_DT$region)
      
      comb_empPval_file <- file.path(pipFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      stopifnot(setequal(all_regs, names(comb_empPval)))
      comb_empPval <- comb_empPval[all_regs]
      # ADJUST THE PVAL
      adj_empPval_comb <- p.adjust(comb_empPval, method="BH")
      
      nGenes_DT <- aggregate(symbol ~ region, data = g2t_DT, FUN=length)
      nHK_DT <- aggregate(is_hk ~ region, data = g2t_DT, FUN=sum)
      adjEmpPval_DT <- data.frame(region=names(adj_empPval_comb), adjEmpPval=adj_empPval_comb, stringsAsFactors = FALSE)
      
      all_DT <- merge(adjEmpPval_DT, merge(nGenes_DT, nHK_DT, by="region"), by = "region")
      stopifnot(all_DT$is_hk <= all_DT$symbol)
      
      colnames(all_DT)[colnames(all_DT) == "symbol"] <- "nbrTotGenes"
      colnames(all_DT)[colnames(all_DT) == "is_hk"] <- "nbrHkGenes"
      
      all_DT$hkPct <- all_DT$nbrHkGenes/all_DT$nbrTotGenes * 100
      
      all_vars <- c("nbrTotGenes", "adjEmpPval", "hkPct")
      xvar <- "nbrHkGenes"
      myxlab <- "# housekeeping genes in TAD"
      myx <- all_DT[,paste0(xvar)]
      
      for(yvar in all_vars) {
        myy <- all_DT[,paste0(yvar)]
        
        if(yvar == "nbrTotGenes") {
          myylab <- "tot # genes in TAD"
        } else if(yvar == "adjEmpPval") {
          myylab <- "adj. combined empPval [-log10]"
          myy <- -log10(myy)
        } else if(yvar == "hkPct") {
          myylab <- "% hk genes in TAD"
        } else {
          stop("error invalid yvar")
        }
        outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", yvar, "_vs_", xvar, "_",  "densplot.", plotType))
        do.call(plotType, list(outFile, height=myHeight, width=myWidth))
        densplot(
          x=myx,
          y=myy,
          xlab = myxlab,
          ylab = myylab,
          main = paste0(yvar, " vs. ", xvar)
        )
        addCorr(x = myx, y = myy, bty="n")
        mtext(side=3, text = paste0(hicds, " - ", exprds))
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
      }
      all_DT$hicds <- hicds
      all_DT$exprds <- exprds
      all_DT
    }
    ds_dt
  }
  
  outFile <- file.path(outFolder, "all_hk_dt.Rdata")
  save(all_hk_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_hk_dt.Rdata")
  all_hk_dt <- eval(parse(text = load(outFile)))
  
}
########################################################################################################################
########################################################################################################################
########################################################################################################################

totDS <- length(unique(paste0(all_hk_dt$hicds, "_", all_hk_dt$exprds)))

all_vars <- c("nbrTotGenes", "adjEmpPval", "hkPct")

xvar <- "nbrHkGenes"
myxlab <- "# housekeeping genes in TAD"
myx <- all_hk_dt[,paste0(xvar)]

for(yvar in all_vars) {
  myy <- all_hk_dt[,paste0(yvar)]
  
  if(yvar == "nbrTotGenes") {
    myylab <- "tot # genes in TAD"
  } else if(yvar == "adjEmpPval") {
    myylab <- "adj. combined empPval [-log10]"
    myy <- -log10(myy)
  } else if(yvar == "hkPct") {
    myylab <- "% hk genes in TAD"
  } else {
    stop("error invalid yvar")
  }
  
  outFile <- file.path(outFolder, paste0("allDS", "_", yvar, "_vs_", xvar, "_",  "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x=myx,
    y=myy,
    xlab = myxlab,
    ylab = myylab,
    main = paste0(yvar, " vs. ", xvar)
  )
  mtext(side=3, text = paste0("allDS - n=", totDS))
  addCorr(x = myx, y = myy, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}
########################################################################################################################
########################################################################################################################
########################################################################################################################

totDS <- length(unique(paste0(all_hk_dt$hicds, "_", all_hk_dt$exprds)))

all_vars <- c("nbrTotGenes", "hkPct", "nbrHkGenes")

xvar <- "adjEmpPval"
myxlab <- "adj. combined empPval [-log10]"
myx <- all_hk_dt[,paste0(xvar)]
myx <- -log10(myx)

for(yvar in all_vars) {
  myy <- all_hk_dt[,paste0(yvar)]
  
  if(yvar == "nbrTotGenes") {
    myylab <- "tot # genes in TAD"
  } else if(yvar == "nbrHkGenes") {
    myylab <-"# housekeeping genes in TAD"
  } else if(yvar == "hkPct") {
    myylab <- "% hk genes in TAD"
  } else {
    stop("error invalid yvar")
  }
  outFile <- file.path(outFolder, paste0("allDS", "_", yvar, "_vs_", xvar, "_",  "densplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x=myx,
    y=myy,
    xlab = myxlab,
    ylab = myylab,
    main = paste0(yvar, " vs. ", xvar)
  )
  mtext(side=3, text = paste0("allDS - n=", totDS))
  addCorr(x = myx, y = myy, bty="n")
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

########################################################################################################################
########################################################################################################################
########################################################################################################################


signifThresh <- 0.05

all_hk_dt$combEmpPvalSignif <- ifelse(all_hk_dt$adjEmpPval <= signifThresh, "signif.", "not signif.")
ylab <- "% hk genes in TAD"

outFile <- file.path(outFolder, paste0("hkPct", "_bysignif_",  "boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(hkPct ~combEmpPvalSignif, data = all_hk_dt, ylab=ylab, main="% hk genes in signif./not signif. TADs" )
mtext(side=3, text=paste0("allDS - n=", totDS, " - pval thresh.=", signifThresh))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


ylab <- "# housekeeping genes in TAD"
outFile <- file.path(outFolder, paste0("nbrHkGenes", "_bysignif_",  "boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(nbrHkGenes ~combEmpPvalSignif, data = all_hk_dt, ylab=ylab, main="# hk genes in signif./not signif. TADs" )
mtext(side=3, text=paste0("allDS - n=", totDS, " - pval thresh.=", signifThresh))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
