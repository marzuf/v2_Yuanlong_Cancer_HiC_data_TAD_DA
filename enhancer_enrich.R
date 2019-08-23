options(scipen=100)

# Rscript enhancer_enrich.R
# Rscript enhancer_enrich.R Gm12878
# Rscript enhancer_enrich.R Hepg2
# Rscript enhancer_enrich.R K562

startTime <- Sys.time()

script_name <- "enhancer_enrich.R"

startTime <- Sys.time()

cat("> START ", script_name," \n")

setDir <- "/media/electron"
setDir <- ""

buildData <- TRUE

require(foreach)
require(doMC)
registerDoMC(40)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

hmmdata <- commandArgs(trailingOnly = TRUE)[1]
if(is.na(hmmdata)) hmmdata <- "Gm12878"
stopifnot(hmmdata %in% c("Gm12878", "Hepg2", "K562"))

script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

pipFolder <- file.path( "PIPELINE","OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))

outFolder <- file.path("ENHANCER_ENRICH", hmmdata)
dir.create(outFolder, recursive=TRUE)

### PREPARE THE CHROM-HMM data

chromHMM_file <- file.path("..", "2_Yuanlong_Cancer_HiC_data_TAD_DA", "enhancer_data", paste0("wgEncodeAwgSegmentationChromhmm", hmmdata, ".bed"))
stopifnot(file.exists(chromHMM_file))
# bin	chrom	chromStart	chromEnd	name	score	strand	thickStart	thickEnd	itemRgb # ! is 0-based !
chromHMM_DT <- read.delim(chromHMM_file, header=FALSE, stringsAsFactors = FALSE, 
                          col.names=c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb"))

signifThresh <- 0.05


# Tss, TssF	Bright Red	Active Promoter
# PromF	Light Red	Promoter Flanking
# PromP	Purple	Inactive Promoter
# Enh, EnhF	Orange	Candidate Strong enhancer
# EnhWF, EnhW, DNaseU, DNaseD, FaireW	Yellow	Candidate Weak enhancer/DNase
# CtrcfO, Ctcf	Blue	Distal CTCF/Candidate Insulator
# Gen5', Elon, ElonW, Gen3', Pol2, H4K20	Dark Green	Transcription associated
# Low	Light Green	Low activity proximal to active states
# ReprD, Repr, ReprW	Gray	Polycomb repressed
# Quies, Art	Light Gray	Heterochromatin/Repetitive/Copy Number Variation

chromHMMannot <- setNames(c(
  rep("active promoter", 2), "promoter flanking", "inactive promoter", rep("candidate strong enhancer", 2), rep("candidate weak enhancer/DNase", 7),
  rep("distal CTCF/candidate insulator",4), rep("transc. associated", 6), "low activ. prox. to active states", rep("Polycomb repressed", 3),
  rep("Heterochromatin/Repetitive/Copy Number Variation", 2)),
  c("Tss", "TssF", "PromF", "PromP", "Enh", "EnhF", "EnhWF", "EnhW", "DnaseU", "DNaseU", "DNaseD","DnaseD", "FaireW", 
  "Ctrcf0","Ctcf0","CtcfO",  "Ctcf", "Gen5'", "Elon", "ElonW", "Gen3'", "Pol2", "H4K20", "Low", "ReprD", "Repr", "ReprW",
  "Quies", "Art")
)

all_cat <- unique(chromHMMannot)


stopifnot(chromHMM_DT$name %in% names(chromHMMannot))

chromHMM_DT$hmmCategory <- chromHMMannot[chromHMM_DT$name]
chromHMM_DT$start <- chromHMM_DT$start + 1 

chromHMM_DT$hmmCategory <- factor(chromHMM_DT$hmmCategory, levels = all_cat)

chromHMM_DT$hmmBpSize <- chromHMM_DT$end - chromHMM_DT$start + 1

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)


all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))


if(buildData) {

  hicds="GSE105381_HepG2_40kb"
  all_chromHMM_tadCount_DT <- foreach(hicds = all_hicds, .combine='rbind') %do% {
    exprds="TCGAlihc_norm_lihc"  
    ds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, "\n")
      
      g2tFile <- file.path(dirname(dirname(pipFolder)), hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      
      tadposFile <- file.path(dirname(dirname(pipFolder)), hicds, "genes2tad", "all_assigned_regions.txt")
      stopifnot(file.exists(tadposFile))
      tadpos_DT <- read.delim(tadposFile,
                              header=F, col.names = c("chromo",  "region", "start", "end"), stringsAsFactors = FALSE)
      
      geneList_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      geneList <- eval(parse(text = load(geneList_file)))
      stopifnot(geneList %in% g2t_DT$entrezID)
      
      regionList_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_regionList.Rdata")
      stopifnot(file.exists(regionList_file))
      regionList <- eval(parse(text = load(regionList_file)))
      
      g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
      stopifnot(g2t_DT$entrezID %in% entrez2symb_dt$entrezID)
      stopifnot(setequal(g2t_DT$region, regionList))
      all_regs <- regionList
      stopifnot(all_regs %in% tadpos_DT$region)  
      
      comb_empPval_file <- file.path(pipFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      stopifnot(setequal(all_regs, names(comb_empPval)))
      comb_empPval <- comb_empPval[all_regs]
      # ADJUST THE PVAL
      adj_empPval_comb <- p.adjust(comb_empPval, method="BH")
      stopifnot(names(adj_empPval_comb) == all_regs)
      
     all_reg_dt <- foreach(reg = all_regs, .combine='rbind') %dopar% {
        curr_chrom <- tadpos_DT$chromo[tadpos_DT$region == reg]
        stopifnot(length(curr_chrom) == 1)
        curr_start <- tadpos_DT$start[tadpos_DT$region == reg]
        stopifnot(length(curr_start) == 1)
        stopifnot(is.numeric(curr_start) )
        curr_end <- tadpos_DT$end[tadpos_DT$region == reg]
        stopifnot(length(curr_end) == 1)
        stopifnot(is.numeric(curr_end) )
        tad_size <- curr_end - curr_start + 1
        
        # select the tracks # !!! WHAT TO DO WITH OVERLAP BOUNDARY ???
        curr_hmm <- chromHMM_DT[chromHMM_DT$chr == curr_chrom &
                                  chromHMM_DT$start >= curr_start & 
                                  chromHMM_DT$end <= curr_end,]
        if(nrow(curr_hmm) > 0) {
          count_hmm <- aggregate(hmmBpSize ~ hmmCategory, data = curr_hmm, FUN=sum)
        } else {
          count_hmm <- data.frame(hmmCategory=character(0), hmmBpSize = double(0))
        }
        
        missing_cat <- all_cat[!all_cat %in% count_hmm$hmmCategory]
        
        missing_hmm <- data.frame(hmmCategory=missing_cat, hmmBpSize=rep(0, length(missing_cat)))
        
        totCount_DT_bp <- rbind(count_hmm, missing_hmm)
        
        stopifnot(!duplicated(totCount_DT_bp$hmmCategory))
        rownames(totCount_DT_bp) <- totCount_DT_bp$hmmCategory
        totCount_DT_bp$hmmCategory <- NULL
        totCount_DT_bp <- t(totCount_DT_bp)
        colnames(totCount_DT_bp) <- paste0(colnames(totCount_DT_bp), " (bp)")
        rownames(totCount_DT_bp) <- NULL
        totCount_DT_ratio <- totCount_DT_bp/tad_size 
        colnames(totCount_DT_ratio) <- paste0(colnames(totCount_DT_ratio), " (ratio)")
        rownames(totCount_DT_ratio) <- NULL
        
        
        totCount_DT <- cbind(totCount_DT_bp, totCount_DT_ratio)
        totCount_DT <- data.frame(totCount_DT, check.names = FALSE)
        
        totCount_DT$hicds <- hicds
        totCount_DT$exprds <- exprds
        totCount_DT$region <- reg
        
        stopifnot(reg %in% names(adj_empPval_comb))
        
        totCount_DT$region_adjPvalComb <- adj_empPval_comb[paste0(reg)]
        
        totCount_DT
      } # end-foreach iterating over regions
     all_reg_dt
    } # end-foreach iterating over exprds
    ds_dt
  } # end-foreach iterating over hicds
  outFile <- file.path(outFolder, "all_chromHMM_tadCount_DT.Rdata")
  save(all_chromHMM_tadCount_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_chromHMM_tadCount_DT.Rdata")
  all_chromHMM_tadCount_DT <- eval(parse(text = load(outFile)))
}

############################################################

totDS <- length(unique(paste0(all_chromHMM_tadCount_DT$hicds, "_", all_chromHMM_tadCount_DT$exprds)))

xvar <- "region_adjPvalComb"

all_vars <- colnames(all_chromHMM_tadCount_DT)
all_vars <- all_vars[! all_vars %in% c("hicds", "exprds", "region", xvar)]

myxlab <- "adj. empPval. comb [-log10]"
myx <- all_chromHMM_tadCount_DT[,paste0(xvar)]
myx <- -log10(myx)

for(yvar in all_vars) {
  myy <- all_chromHMM_tadCount_DT[,paste0(yvar)]
  myylab <- yvar

  outFile <- file.path(outFolder, paste0("allDS", "_", gsub("\\(|\\)|\\/","_", gsub(" ", "", yvar)), "_vs_", gsub("\\(|\\)|\\/","_", gsub(" ", "", xvar)), "_",  "densplot.", plotType))
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
  
  outFile <- file.path(outFolder, paste0("allDS", "_", gsub("\\(|\\)|\\/","_", gsub(" ", "", yvar)), "_vs_", gsub("\\(|\\)|\\/","_", gsub(" ", "", xvar)), "_",  "densplot_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  densplot(
    x=myx,
    y=log10(myy),
    xlab = myxlab,
    ylab = paste0(myylab, " [log10]"),
    main = paste0(yvar, " vs. ", xvar)
  )
  mtext(side=3, text = paste0("allDS - n=", totDS))
  addCorr(x = myx, y = myy, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}


all_chromHMM_tadCount_DT$combEmpPvalSignif <- ifelse(all_chromHMM_tadCount_DT$region_adjPvalComb <= signifThresh, "signif.", "not signif.")

xvar <- "combEmpPvalSignif"

colnames(all_chromHMM_tadCount_DT) <- gsub("\\(|\\)|\\/","_", gsub(" ", "", colnames(all_chromHMM_tadCount_DT)))

for(yvar in all_vars) {
  myy <- all_chromHMM_tadCount_DT[,paste0( gsub("\\(|\\)|\\/","_", gsub(" ", "", yvar)))]
  myylab <- yvar
  
  outFile <- file.path(outFolder, paste0(gsub("\\(|\\)|\\/","_", gsub(" ", "", yvar)), "_bysignif_",  "boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(as.formula(paste0(gsub("\\(|\\)|\\/","_", gsub(" ", "", yvar)), " ~ ", xvar)), 
          data = all_chromHMM_tadCount_DT, 
          ylab=myylab, main=paste0(yvar, " in signif./not signif. TADs" ))
  mtext(side=3, text=paste0("allDS - n=", totDS, " - pval thresh.=", signifThresh))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  all_chromHMM_tadCount_DT[,paste0(gsub("\\(|\\)|\\/","_", gsub(" ", "", yvar)), "_log10")] <- log10(all_chromHMM_tadCount_DT[,paste0(gsub("\\(|\\)|\\/","_", gsub(" ", "", yvar)))])
  
  outFile <- file.path(outFolder, paste0(gsub("\\(|\\)|\\/","_", gsub(" ", "", yvar)), "_bysignif_",  "boxplot_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  boxplot(as.formula(paste0(gsub("\\(|\\)|\\/","_", gsub(" ", "", yvar), "_log10"), " ~ ", xvar)), 
          data = all_chromHMM_tadCount_DT, 
          ylab=paste0(myylab, " [log10]"), main=paste0(yvar, " in signif./not signif. TADs" ))
  mtext(side=3, text=paste0("allDS - n=", totDS, " - pval thresh.=", signifThresh))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}
# 
# all_chromHMM_tadCount_DT$combEmpPvalSignif <- ifelse(all_chromHMM_tadCount_DT$region_adjPvalComb <= signifThresh, "signif.", "not signif.")
# ylab <- "% hk genes in TAD"
# 
# outFile <- file.path(outFolder, paste0("hkPct", "_bysignif_",  "boxplot.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# boxplot(hkPct ~combEmpPvalSignif, data = all_chromHMM_tadCount_DT, ylab=ylab, main="% hk genes in signif./not signif. TADs" )
# mtext(side=3, text=paste0("allDS - n=", ntotDS, " - pval thresh.=", signifThresh))
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# 
# ylab <- "# housekeeping genes in TAD"
# outFile <- file.path(outFolder, paste0("nbrHkGenes", "_bysignif_",  "boxplot.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# boxplot(nbrHkGenes ~combEmpPvalSignif, data = all_chromHMM_tadCount_DT, ylab=ylab, main="# hk genes in signif./not signif. TADs" )
# mtext(side=3, text=paste0("allDS - n=", ntotDS, " - pval thresh.=", signifThresh))
# foo <- dev.off()
# cat(paste0("... written: ", outFile, "\n"))



##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))
