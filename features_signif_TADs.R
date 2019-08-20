options(scipen=100)

# Rscript features_signif_TADs.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript features_signif_TADs.R   # to run all datasets in one shot

setDir=""

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "features_signif_TADs.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

buildData <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(ggplot2)
require(reshape2)
require(grid)
require(lattice)
require(RColorBrewer)
hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script8c_name <- "8cOnlyRatioDownFastSave_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script11same_name <- "11sameNbr_runEmpPvalCombined"
script19_name <- "19onlyFC_SAM_emp_measurement" 
script19sameNbr_name <- "19sameNbr_SAM_emp_measurement"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))


outFolder <- "FEATURES_SIGNIF_TADS"
dir.create(outFolder, recursive=TRUE)


FDRthresh_seq <- seq(from=0.1, to=0.5, by=0.1)
pvalThresh_seq <- seq(from=0.01, to=0.05, by = 0.01)

myHeightGG <- length(pvalThresh_seq)*1.2
myWidthGG <- length(FDRthresh_seq)*1.2

twoSidedStouffer <- FALSE

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


familyType <- "hgnc_entrezID"
fam_type <- "hgnc_family_short"

if(buildData) {

  ### BUILD RATIO DOWN
  cat("... start building ratiodown and FC \n")
  allDS_rD_FC_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_ratioDown_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      rd_file <- file.path(pipOutFolder, hicds, exprds, script8c_name, "all_obs_ratioDown.Rdata")
      stopifnot(file.exists(rd_file))
      tad_rD <- eval(parse(text = load(rd_file)))
      all_regs <- names(tad_rD)
      
      
      fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
      corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(fc_file))
      stopifnot(file.exists(corr_file))  
      tad_logFC <- eval(parse(text = load(fc_file)))
      tad_meanCorr <- eval(parse(text = load(corr_file)))
      stopifnot(setequal(names(tad_logFC), all_regs))
      stopifnot(setequal(names(tad_meanCorr), all_regs))
      
      family_file <- file.path(pipFolder, "PREP_GENE_FAMILIES_TAD_DATA", hicds,  paste0(familyType, "_family_TAD_DT.Rdata"))
      stopifnot(file.exists(family_file))
      family_DT <- eval(parse(text = load(family_file)))
      agg_fam_DT <- aggregate(formula(paste0(fam_type, " ~ region")), data = family_DT, FUN=function(x)length(unique(x)))
      tad_nFams <- setNames(agg_fam_DT$hgnc_family_short, agg_fam_DT$region)
      stopifnot(is.numeric(tad_nFams))
      
      
      de_file <- file.path(pipOutFolder, hicds, exprds, script1_name, "DE_topTable.Rdata")
      stopifnot(file.exists(de_file))
      de_DT <- eval(parse(text = load(de_file)))
      de_DT$genes <- as.character(de_DT$genes)
      g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      geneList_file <- file.path(pipOutFolder, hicds, exprds,script0_name,   "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      geneList <- eval(parse(text = load(geneList_file)))
      stopifnot(names(geneList) %in% de_DT$genes)
      stopifnot(geneList %in% g2t_DT$entrezID)
      
      g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
      stopifnot(setequal(g2t_DT$entrezID, geneList))
      
      absMaxLogFC <- sapply(all_regs, function(reg) {
        g2t_tad_id <- g2t_DT$entrezID[g2t_DT$region==reg]
        de_tad_id <- names(geneList)[geneList %in% g2t_tad_id]
        stopifnot(de_tad_id %in% de_DT$genes)
        max(abs(de_DT$logFC[de_DT$genes %in% de_tad_id]))
      })
      names(absMaxLogFC) <- all_regs
      sdLogFC <- sapply(all_regs, function(reg) {
        g2t_tad_id <- g2t_DT$entrezID[g2t_DT$region==reg]
        de_tad_id <- names(geneList)[geneList %in% g2t_tad_id]
        stopifnot(de_tad_id %in% de_DT$genes)
        sd(de_DT$logFC[de_DT$genes %in% de_tad_id])
      })    
      names(sdLogFC) <- all_regs
      
      meanLogFC_v1 <- sapply(all_regs, function(reg) {
        g2t_tad_id <- g2t_DT$entrezID[g2t_DT$region==reg]
        de_tad_id <- names(geneList)[geneList %in% g2t_tad_id]
        stopifnot(de_tad_id %in% de_DT$genes)
        mean(de_DT$logFC[de_DT$genes %in% de_tad_id])
      })    
      names(meanLogFC_v1) <- all_regs
      
      
      ### ADD TAD SIZE AND NUMBER OF GENES
      nGenes_tads <- setNames(as.numeric(table(g2t_DT$region)), as.character(names(table(g2t_DT$region))))
      stopifnot(setequal(names(nGenes_tads), all_regs))
      nGenes_tads <- nGenes_tads[all_regs]
      
      
      ### RETRIEVE THE TAD POSITIONS
      tadposFile <- file.path(pipFolder, hicds, "genes2tad", "all_assigned_regions.txt")
      stopifnot(file.exists(tadposFile))
      tadpos_DT <- read.delim(tadposFile, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
      stopifnot(is.numeric(tadpos_DT$start))
      stopifnot(is.numeric(tadpos_DT$end))
      tadpos_DT <- tadpos_DT[grepl("_TAD", tadpos_DT$region),,drop=FALSE] 
      stopifnot(all_regs %in% tadpos_DT$region)
      tadpos_DT <- tadpos_DT[tadpos_DT$region %in% all_regs,]
      stopifnot(setequal(all_regs, tadpos_DT$region))
      tadpos_DT$tad_size <- tadpos_DT$end-tadpos_DT$start+1
      tadSize <- setNames(tadpos_DT$tad_size, tadpos_DT$region)
      stopifnot(setequal(all_regs, names(tadSize)))
      tadSize <- tadSize[all_regs]
      
      stopifnot(names(tad_rD) == all_regs)
      stopifnot(names(tad_rD) == names(sdLogFC))
      stopifnot(names(tad_rD) == names(absMaxLogFC))
      stopifnot(names(tad_rD) == names(meanLogFC_v1))
      stopifnot(names(tad_rD) == names(nGenes_tads))
      stopifnot(names(tad_rD) == names(tadSize))
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        region = all_regs,
        ratioDown = as.numeric(tad_rD[all_regs]),
        meanCorr = as.numeric(tad_meanCorr[all_regs]),
        meanLogFC = as.numeric(tad_logFC[all_regs]),
        meanLogFC_v1 = as.numeric(tad_logFC[all_regs]),
        maxAbsLogFC = as.numeric(absMaxLogFC[all_regs]),
        sdLogFC = as.numeric(sdLogFC[all_regs]),
        nGenes = as.numeric(nGenes_tads[all_regs]),
        TADsize = as.numeric(tadSize[all_regs]),
        nFamilies = as.numeric(tad_nFams[all_regs]),
        stringsAsFactors = FALSE
      )
    } # end-foreach iterating exprds
  } # end-foreach iterating hicds
  
  
  outFile <- file.path(outFolder, "allDS_rD_FC_DT.Rdata")  
  save(allDS_rD_FC_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  stopifnot(nrow(allDS_rD_FC_DT) > 0)
  
  ### BUILD SIGNIF ALONG FDR THRESH
  cat("... start building signif. along FDR thresh\n")
  allDS_signifFDR_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_signifFDR_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      # PREPARE logFC and meanCorr observed data
      fc_file <-  file.path(pipOutFolder, hicds, exprds, script3_name, "all_meanLogFC_TAD.Rdata")
      corr_file <- file.path(pipOutFolder, hicds, exprds, script4_name, "all_meanCorr_TAD.Rdata")
      stopifnot(file.exists(fc_file))
      stopifnot(file.exists(corr_file))  
      tad_logFC <- eval(parse(text = load(fc_file)))
      tad_meanCorr <- eval(parse(text = load(corr_file)))
      all_regs <- names(tad_logFC)
      stopifnot(setequal(all_regs, names(tad_meanCorr)))
      tad_logFC <- tad_logFC[all_regs]
      tad_meanCorr <- tad_meanCorr[all_regs]
      # RETRIEVE FDR DATA FOR MEAN LOGFC
      logFC_FDR_file <- file.path(pipOutFolder, hicds, exprds, script19_name, "empFDR_list.Rdata")
      stopifnot(file.exists(logFC_FDR_file))
      all_FDR <- eval(parse(text = load(logFC_FDR_file)))
      logFC_FDR <- all_FDR[["empFDR_logFC"]]  # the names is the FC threshold, the value is the FDR
      stopifnot(length(logFC_FDR) > 0)
      # RETRIEVE FDR DATA FOR MEAN CORR



      # RETRIEVE FDR DATA FOR MEAN CORR
      # the same for meanCorr
      meanCorr_FDR_file <-  file.path(pipOutFolder, hicds, exprds, script19sameNbr_name, "meanCorr_empFDR.Rdata")
      stopifnot(file.exists(meanCorr_FDR_file))
      all_corr_FDR <- eval(parse(text = load(meanCorr_FDR_file)))
      meanCorr_FDR <- all_corr_FDR[["empFDR"]]  # the names is the meanCorr threshold, the value is the FDR
      stopifnot(length(meanCorr_FDR) > 0)

      # => SIGNIF TADs FOR VARIOUS FDR THRESHOLD
      # for each of the FDR threshold => FC cut-off, meanCorr cut-off => signif TADs
      cutoff_fdr <- FDRthresh_seq[1]
      signif_tads_dt <- foreach(cutoff_fdr = FDRthresh_seq , .combine='rbind') %do% {
        logFC_cut_off <- min(as.numeric(as.character(na.omit(names(logFC_FDR)[logFC_FDR <= cutoff_fdr]))))  # the smallest FC cut-off that leads to desired FDR; if not returns Inf
        meanCorr_cut_off <- min(as.numeric(as.character(na.omit(names(meanCorr_FDR)[meanCorr_FDR <= cutoff_fdr]))))
        stopifnot(names(tad_logFC) == names(tad_meanCorr))
        signif_tads <- ( abs(tad_logFC) >= logFC_cut_off & tad_meanCorr >= meanCorr_cut_off)
        stopifnot(names(signif_tads) == names(tad_meanCorr))
        stopifnot(names(signif_tads) == all_regs)
        data.frame(
          hicds = hicds,
          exprds = exprds,
          region = names(signif_tads),
          thresh_FDR = cutoff_fdr,
          signif_FDR = as.logical(signif_tads), # to avoid having rownames
          stringsAsFactors = FALSE
        )
      } # end-foreach iterating over FDR threshs
    signif_tads_dt
    } # end-foreach iterating over exprds
    exprds_signifFDR_DT
  } # end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, "allDS_signifFDR_DT.Rdata")  
  save(allDS_signifFDR_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  stopifnot(nrow(allDS_signifFDR_DT) > 0 )
  nIDs <- unique(paste0(allDS_signifFDR_DT$hicds, allDS_signifFDR_DT$exprds, allDS_signifFDR_DT$region))
  stopifnot(length(nIDs) == nrow(allDS_rD_FC_DT))
  
  ### BUILD SIGNIF ALONG PVAL THRESH
  cat("... start building signif. along pval thresh\n")
  allDS_signifPval_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_signifPval_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {


      # RETRIEVE COMBINED EMP PVAL      
      comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      all_regs <- names(comb_empPval)
      # ADJUST THE PVAL
      adj_empPval_comb <- p.adjust(comb_empPval, method="BH")
      stopifnot(names(adj_empPval_comb) == all_regs)

      # => SIGNIF TADs FOR VARIOUS PVAL THRESHOLD
      # for each of the threshold of empPvals -> retrieve signif TADs
      cutoff_pval=pvalThresh_seq[1]
      
      signif_tads_dt <- foreach(cutoff_pval = pvalThresh_seq , .combine='rbind') %do% {
        signif_tads_pval <- (adj_empPval_comb <= cutoff_pval)
        stopifnot(names(signif_tads_pval) == names(adj_empPval_comb))
        stopifnot(names(signif_tads_pval) == all_regs)
        data.frame(
          hicds = hicds,
          exprds = exprds,
          region = names(signif_tads_pval),
          thresh_pval = cutoff_pval,
          signif_pval = as.logical(signif_tads_pval), # to avoid having rownames
          stringsAsFactors = FALSE
        )
      } # end-foreach iterating over FDR threshs
      signif_tads_dt
    } # end-foreach iterating over exprds
    exprds_signifPval_DT
  } # end-foreach iterating over hicds
  
  
  outFile <- file.path(outFolder, "allDS_signifPval_DT.Rdata")  
  save(allDS_signifPval_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  stopifnot(nrow(allDS_signifPval_DT) > 0 )
  
  nIDs <- unique(paste0(allDS_signifPval_DT$hicds, allDS_signifPval_DT$exprds, allDS_signifPval_DT$region))
  stopifnot(length(nIDs) == nrow(allDS_rD_FC_DT))
  
  
  
  cat("... merge signifFDR <-> ratioDown\n")
  all_dt_signifFDR <- merge(allDS_signifFDR_DT, allDS_rD_FC_DT, by =c("hicds", "exprds", "region"), all = TRUE)
  stopifnot(nrow(all_dt_signifFDR) == nrow(allDS_signifFDR_DT))
  
  outFile <- file.path(outFolder, "all_dt_signifFDR.Rdata")  
  save(all_dt_signifFDR, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  cat("... merge signifFDR <-> ratioDown\n")
  all_dt_signifPval <- merge(allDS_signifPval_DT, allDS_rD_FC_DT, by =c("hicds", "exprds", "region"), all = TRUE)
  stopifnot(nrow(all_dt_signifPval) == nrow(allDS_signifPval_DT))
  
  outFile <- file.path(outFolder, "all_dt_signifPval.Rdata")  
  save(all_dt_signifPval, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
} else { # end-if build data
  outFile <- file.path(outFolder, "all_dt_signifPval.Rdata")  
  all_dt_signifPval <- eval(parse(text = load(outFile)))
  
  outFile <- file.path(outFolder, "all_dt_signifFDR.Rdata")  
  all_dt_signifFDR <- eval(parse(text = load(outFile)))
  
}

all_dt_signifPval$adjCombPval <- ifelse(all_dt_signifPval$signif_pval, "signif.", "not signif.")
all_dt_signifPval$thresh_pval_label <- paste0("P=", all_dt_signifPval$thresh_pval)

all_dt_signifFDR$FDR <- ifelse(all_dt_signifFDR$signif_FDR, "signif.", "not signif.")
all_dt_signifFDR$thresh_FDR_label <- paste0("FDR=", all_dt_signifFDR$thresh_FDR)


all_vars <- c("ratioDown", "meanCorr", "meanLogFC", "maxAbsLogFC", "sdLogFC", "nGenes", "TADsize", "nFamilies")


### PLOT THE SIGNIF PVAL FEATURES => compare dist Pval signif. vs. not signif.
plot_var=all_vars[1]
for(plot_var in all_vars) {
  outFile <- file.path(outFolder, paste0("allDS_signifPval_", plot_var, "_signif_vs_notsignif_density_lattice.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*3))
  myplot <- densityplot( formula(paste0("~",plot_var, "| thresh_pval_label")), groups = adjCombPval, data = all_dt_signifPval, #auto.key = TRUE,
              # par.strip.text=list(cex=1), # width of the strip bar
              par.strip.text = list(cex = 1, font = 4, col = "brown"),
              layout = c(length(unique(all_dt_signifPval$thresh_pval_label)), 1), # column,row
              scales=list(y=list(relation="free"),
                          x=list(relation="free")
              ),
              auto.key=list(title="", space = "bottom", cex=1.0, columns=length(unique(all_dt_signifPval$adjCombPval))),
              main = paste0(plot_var))
  print(myplot)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

### PLOT THE SIGNIF FDR FEATURES => compare dist FDR signif. vs. not signif.
plot_var=all_vars[1]
for(plot_var in all_vars) {
  outFile <- file.path(outFolder, paste0("allDS_signifFDR_", plot_var, "_signif_vs_notsignif_density_lattice.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth*3))
  myplot <- densityplot( formula(paste0("~",plot_var, "| thresh_FDR_label")), groups = FDR, data = all_dt_signifFDR, #auto.key = TRUE,
                         # par.strip.text=list(cex=1), # width of the strip bar
                         par.strip.text = list(cex = 1, font = 4, col = "brown"),
                         layout = c(length(unique(all_dt_signifFDR$thresh_FDR_label)), 1), # column, row
                         scales=list(y=list(relation="free"),
                                     x=list(relation="free")
                         ),
                         auto.key=list(title="", space = "bottom", cex=1.0, columns=length(unique(all_dt_signifFDR$FDR))),
                         main = paste0(plot_var))
  print(myplot)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}

### CMP SIGNIF INTERSECT NOT INTERSECT
  signifFDR_dt <- all_dt_signifFDR[all_dt_signifFDR$signif_FDR,]
  signifPval_dt <- all_dt_signifPval[all_dt_signifPval$signif_pval,]
  stopifnot(signifFDR_dt$FDR == "signif.")
  stopifnot(signifPval_dt$adjCombPval == "signif.")


  signif_FDR_pval_dt <- merge(signifFDR_dt, signifPval_dt, by = c("hicds", "exprds", "region"), all = TRUE)

  check_dt <- signif_FDR_pval_dt[!is.na(signif_FDR_pval_dt$signif_FDR) & !is.na(signif_FDR_pval_dt$signif_pval),]

  for(var in all_vars) {
    if(var != "nFamilies") stopifnot(check_dt[,paste0(var, ".x")] == check_dt[,paste0(var, ".y")])
    if(var == "nFamilies") stopifnot(na.omit(check_dt[,paste0(var, ".x")]) == na.omit(check_dt[,paste0(var, ".y")]))
    signif_FDR_pval_dt[,paste0(var)] <- ifelse(is.na(signif_FDR_pval_dt[,paste0(var, ".x")]), signif_FDR_pval_dt[,paste0(var, ".y")], signif_FDR_pval_dt[,paste0(var, ".x")] )
  }
  check_dt <- signif_FDR_pval_dt[!is.na(signif_FDR_pval_dt$signif_FDR) & !is.na(signif_FDR_pval_dt$signif_pval),]
  for(var in all_vars) {
    if(var != "nFamilies") stopifnot(check_dt[,paste0(var, ".x")] == check_dt[,paste0(var, ".y")])
    if(var != "nFamilies") stopifnot(check_dt[,paste0(var, ".x")] == check_dt[,paste0(var)])
    if(var == "nFamilies") stopifnot(na.omit(check_dt[,paste0(var, ".x")]) == na.omit(check_dt[,paste0(var, ".y")]))
    if(var == "nFamilies") stopifnot(na.omit(check_dt[,paste0(var, ".x")]) == na.omit(check_dt[,paste0(var)]))
    signif_FDR_pval_dt[,paste0(var, ".x")] <- NULL
    signif_FDR_pval_dt[,paste0(var, ".y")] <- NULL
  }
  stopifnot( !( is.na(signif_FDR_pval_dt$FDR) & is.na(signif_FDR_pval_dt$adjCombPval) ) ) # should never happen that the 2 are NA
  signif_FDR_pval_dt$signif_label <- ifelse( is.na(signif_FDR_pval_dt$FDR), "onlyPval",
                                             ifelse( is.na(signif_FDR_pval_dt$adjCombPval), "onlyFDR", "both"))

  stopifnot(signif_FDR_pval_dt$FDR[signif_FDR_pval_dt$signif_label == "both"] == "signif.")
  stopifnot(signif_FDR_pval_dt$FDR[signif_FDR_pval_dt$signif_label == "onlyFDR"] == "signif.")
  stopifnot(signif_FDR_pval_dt$adjCombPval[signif_FDR_pval_dt$signif_label == "both"] == "signif.")
  stopifnot(signif_FDR_pval_dt$adjCombPval[signif_FDR_pval_dt$signif_label == "onlyPval"] == "signif.")

  # signif_FDR_pval_dt <- signif_FDR_pval_dt[order(signif_FDR_pval_dt$thresh_FDR_label, signif_FDR_pval_dt$thresh_pval_label),]
  # signif_FDR_pval_dt$both_label <- paste0(signif_FDR_pval_dt$thresh_FDR_label, ";", signif_FDR_pval_dt$thresh_pval_label )
  signif_FDR_pval_dt <- signif_FDR_pval_dt[order(signif_FDR_pval_dt$thresh_pval_label, signif_FDR_pval_dt$thresh_FDR_label),]
  signif_FDR_pval_dt$both_label <- paste0(signif_FDR_pval_dt$thresh_pval_label, ";",  signif_FDR_pval_dt$thresh_FDR_label)

  signif_FDR_pval_dt <- signif_FDR_pval_dt[! grepl("^NA", signif_FDR_pval_dt$both_label) & !grepl("NA$", signif_FDR_pval_dt$both_label), ]



  nCases <- length(unique(signif_FDR_pval_dt$both_label))
  nCol <- 5
  nRow <- ceiling(nCases/nCol)

  groupCols <- trellis.par.get("superpose.line")$col[1:length(unique(signif_FDR_pval_dt$signif_label))]
  names(groupCols) <- levels(as.factor(signif_FDR_pval_dt$signif_label))
  

  plot_var=all_vars[1]
  for(plot_var in all_vars) {
    outFile <- file.path(outFolder, paste0("allDS_", plot_var, "_intersect_vs_only_signif_density_lattice.", plotType))
    do.call(plotType, list(outFile, height=myHeight*nRow, width=myWidth*nCol))
    # myplot <- densityplot( formula(paste0("~",plot_var, "| both_label")), groups = signif_label, data = signif_FDR_pval_dt, #auto.key = TRUE,
    #                        # par.strip.text=list(cex=1), # width of the strip bar
    #                        par.strip.text = list(cex = 1, font = 4, col = "brown"),
    #                        # layout = c(length(unique(signif_FDR_pval_dt$both_label))/4, 1), # column, row
    #                        layout = c(nCol, nRow), # column, row
    #                        scales=list(y=list(relation="free"),
    #                                    x=list(relation="free")
    #                        ),
    #                        auto.key=list(title="", space = "bottom", cex=1.0, columns=length(unique(signif_FDR_pval_dt$signif_label))),
    #                        main = paste0(plot_var))
    
    plotDT <- signif_FDR_pval_dt
    if(plot_var=="nFamilies") {
      plotDT <- plotDT[,c(plot_var, "both_label", "signif_label")]
      plotDT <- na.omit(signif_FDR_pval_dt) 
    }
    
    myplot <- densityplot( formula(paste0("~",plot_var, "| both_label")), groups = signif_label, data = plotDT, #auto.key = TRUE, 
                           par.strip.text = list(cex = 1, font = 4, col = "brown"),
                           layout = c(nCol, nRow), # column, row
                           scales=list(y=list(relation="free"),
                                       x=list(relation="free")
                           ),
                           panel = function(x,groups,subscripts, ...) {
                             panel.densityplot(x,groups=groups, subscripts=subscripts)
                             obs_nbr <-  as.numeric(table(groups[subscripts]))
                             obs_cols <- groupCols[names(table(groups[subscripts]))]
                             draw.key(key = list(text = list(#as.character(length(x)), 
                               paste0("n=", as.character(obs_nbr))), col= obs_cols), 
                               draw=TRUE,  
                               vp=viewport(x=0.9,y=0.9))
                           },
                           auto.key=list(title="", space = "bottom", cex=1.0, columns=length(unique(signif_FDR_pval_dt$signif_label))),
                           main = paste0(plot_var))
    
    

    
    print(myplot)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
  }

### CMP SIGNIF INTERSECT NOT INTERSECT
signifFDR_dt <- all_dt_signifFDR
signifPval_dt <- all_dt_signifPval

cat(paste0("... merging signifFDR_dt and signifPval_dt \n"))
signif_FDR_pval_dt <- merge(signifFDR_dt, signifPval_dt, by = c("hicds", "exprds", "region"), all = TRUE)

check_dt <- signif_FDR_pval_dt[!is.na(signif_FDR_pval_dt$signif_FDR) & !is.na(signif_FDR_pval_dt$signif_pval),]

for(var in all_vars) {
  if(var != "nFamilies") stopifnot(check_dt[,paste0(var, ".x")] == check_dt[,paste0(var, ".y")])
  if(var == "nFamilies") stopifnot(na.omit(check_dt[,paste0(var, ".x")] == check_dt[,paste0(var, ".y")]))
  signif_FDR_pval_dt[,paste0(var)] <- ifelse(is.na(signif_FDR_pval_dt[,paste0(var, ".x")]), signif_FDR_pval_dt[,paste0(var, ".y")], signif_FDR_pval_dt[,paste0(var, ".x")] )
}
check_dt <- signif_FDR_pval_dt[!is.na(signif_FDR_pval_dt$signif_FDR) & !is.na(signif_FDR_pval_dt$signif_pval),]
for(var in all_vars) {
  if(var != "nFamilies") stopifnot(check_dt[,paste0(var, ".x")] == check_dt[,paste0(var, ".y")])
  if(var != "nFamilies") stopifnot(check_dt[,paste0(var, ".x")] == check_dt[,paste0(var)])
  if(var == "nFamilies") stopifnot(na.omit(check_dt[,paste0(var, ".x")]) == na.omit(check_dt[,paste0(var, ".y")]))
  if(var == "nFamilies") stopifnot(na.omit(check_dt[,paste0(var, ".x")]) == na.omit(check_dt[,paste0(var)]))
  signif_FDR_pval_dt[,paste0(var, ".x")] <- NULL
  signif_FDR_pval_dt[,paste0(var, ".y")] <- NULL
}
stopifnot( !( is.na(signif_FDR_pval_dt$FDR) & is.na(signif_FDR_pval_dt$adjCombPval) ) ) # should never happen that the 2 are NA
signif_FDR_pval_dt$signif_label <- ifelse( signif_FDR_pval_dt$signif_FDR & signif_FDR_pval_dt$signif_pval, "bothSignif.",
                                              ifelse( ! signif_FDR_pval_dt$signif_FDR & ! signif_FDR_pval_dt$signif_pval, "noneSignif.", 
                                                ifelse( signif_FDR_pval_dt$signif_FDR , "FDRsignif.", 
                                                       ifelse( signif_FDR_pval_dt$signif_pval , "adjCombPvalSsignif.", NA))))
stopifnot(!is.na(signif_FDR_pval_dt$signif_label))
                               
signif_FDR_pval_dt$signif_label <- factor(as.character(signif_FDR_pval_dt$signif_label),
                                          levels = c("bothSignif.","FDRsignif.", "adjCombPvalSsignif.", "noneSignif." ))
            

# signif_FDR_pval_dt <- signif_FDR_pval_dt[order(signif_FDR_pval_dt$thresh_FDR_label, signif_FDR_pval_dt$thresh_pval_label),]
# signif_FDR_pval_dt$both_label <- paste0(signif_FDR_pval_dt$thresh_FDR_label, ";", signif_FDR_pval_dt$thresh_pval_label )

signif_FDR_pval_dt <- signif_FDR_pval_dt[order(signif_FDR_pval_dt$thresh_pval_label, signif_FDR_pval_dt$thresh_FDR_label),]
signif_FDR_pval_dt$both_label <- paste0(signif_FDR_pval_dt$thresh_pval_label, ";",  signif_FDR_pval_dt$thresh_FDR_label)

plot_var=all_vars[1]

nCases <- length(unique(signif_FDR_pval_dt$both_label))
nCol <- 5
nRow <- ceiling(nCases/nCol)

# for(plot_var in all_vars) {
#   outFile <- file.path(outFolder, paste0("allDS_", plot_var, "_signif_by_thresh_density_lattice.", plotType))
#   do.call(plotType, list(outFile, height=myHeight*nRow, width=myWidth*nCol))
#   myplot <- densityplot( formula(paste0("~",plot_var, "| both_label")), groups = signif_label, data = signif_FDR_pval_dt, #auto.key = TRUE, 
#                          # par.strip.text=list(cex=1), # width of the strip bar
#                          par.strip.text = list(cex = 1, font = 4, col = "brown"),
#                          # layout = c(length(unique(signif_FDR_pval_dt$both_label))/4, 1), # column, row
#                          layout = c(nCol, nRow), # column, row
#                          scales=list(y=list(relation="free"),
#                                      x=list(relation="free")
#                          ),
#                          auto.key=list(title="", space = "bottom", cex=1.0, columns=length(unique(signif_FDR_pval_dt$signif_label))),
#                          main = paste0(plot_var))
#   print(myplot)
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))
# }

##########################################################################
##########################################################################
##########################################################################
groupCols <- trellis.par.get("superpose.line")$col[1:length(unique(signif_FDR_pval_dt$signif_label))]
names(groupCols) <- levels(as.factor(signif_FDR_pval_dt$signif_label))


plot_var=all_vars[1]
for(plot_var in all_vars) {
  outFile <- file.path(outFolder, paste0("allDS_", plot_var, "_signif_by_thresh_density_lattice.", plotType))
  do.call(plotType, list(outFile, height=myHeight*nRow, width=myWidth*nCol))
  
  plotDT <- signif_FDR_pval_dt
  
  if(plot_var=="nFamilies") {
    plotDT <- plotDT[,c(plot_var, "both_label", "signif_label")]
    plotDT <- na.omit(signif_FDR_pval_dt) 
  }

  myplot <- densityplot( formula(paste0("~",plot_var, "| both_label")), groups = signif_label, data = plotDT, #auto.key = TRUE, 
               # color=c("red","blue","green", "yellow"),
               # par.strip.text=list(cex=1), # width of the strip bar
               par.strip.text = list(cex = 1, font = 4, col = "brown"),
               # layout = c(length(unique(signif_FDR_pval_dt$both_label))/4, 1), # column, row
               layout = c(nCol, nRow), # column, row
               scales=list(y=list(relation="free"),
                           x=list(relation="free")
               ),
               # str(trellis.par.get(), max.level = 1)
               # str(trellis.par.get("superpose.line"))
               # par.settings = list(superpose.line = list(  col = c("orange", "blue", "red", "yellow"))),
               panel = function(x,groups,subscripts, ...) {
                 panel.densityplot(x,groups=groups, subscripts=subscripts)
                 # foreach panel, groups=the name of the category; subscripts=the index of x corresponding to this panel
                 # cat("groups=", groups, "\n")
                 # cat("subscripts=", subscripts, "\n")
                 # xx <- setNames(as.numeric(table(groups[subscripts])), names(table(groups[subscripts]), "\, "\)
                 # cat("tablepaste:", paste0(names(table(groups[subscripts])),"=", as.numeric(table(groups[subscripts])), "\n"))
                 # cat("cols:", paste0(groupCols[names(table(groups[subscripts]))],"\n"))
                 obs_nbr <-  as.numeric(table(groups[subscripts]))
                 obs_cols <- groupCols[names(table(groups[subscripts]))]
                 # cat("table: ", table(groups[subscripts]), "\n")
                 # n_ind <- tapply(x, groups[subscripts], length) # similar to table
                 # 
                 # draw.key(key = list(text = list(as.character(length(x)), 
                 #                                 as.character(levels(groups)))), draw=TRUE,  vp=viewport(x=0.25,y=0.9))
                 # 
                 draw.key(key = list(text = list(#as.character(length(x)), 
                                                 # paste0("n=", as.character(n_ind))), col= groupCols), 
                          paste0("n=", as.character(obs_nbr))), col= obs_cols), 
                                            draw=TRUE,  
                                            vp=viewport(x=0.9,y=0.9))
               },
               auto.key=list(title="", space = "bottom", cex=1.0, columns=length(unique(signif_FDR_pval_dt$signif_label))),
               main = paste0(plot_var))
  
  print(myplot)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}




##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

