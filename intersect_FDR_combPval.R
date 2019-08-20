options(scipen=100)

# Rscript intersect_FDR_combPval.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript intersect_FDR_combPval.R   # to run all datasets in one shot

setDir=""

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "intersect_FDR_combPval.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

buildData <- TRUE
separateHeatmap <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 10))
require(ggplot2)
require(reshape2)

require(lattice)
require(grid)
require(RColorBrewer)
hm.palette <- colorRampPalette(rev(brewer.pal(9, 'YlOrRd')), space='Lab')

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")

script0_name <- "0_prepGeneData"
script3_name <- "3_runMeanTADLogFC"
script4_name <- "4_runMeanTADCorr"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script19_name <- "19onlyFC_SAM_emp_measurement"
script19sameNbr_name <- "19sameNbr_SAM_emp_measurement"
script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4



pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))



outFolder <- "INTERSECT_FDR_COMBPVAL"
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

if(buildData) {
  allDS_intersect_DT <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds_intersect_DT <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat("... start ", hicds, " - " , exprds, "\n")
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

      # RETRIEVE COMBINED EMP PVAL      
      comb_empPval_file <- file.path(pipOutFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      stopifnot(setequal(all_regs, names(comb_empPval)))
      comb_empPval <- comb_empPval[all_regs]
      # ADJUST THE PVAL
      adj_empPval_comb <- p.adjust(comb_empPval, method="BH")
      stopifnot(names(adj_empPval_comb) == all_regs)


      # => SIGNIF TADs FOR VARIOUS PVAL THRESHOLD
      # for each of the threshold of empPvals -> retrieve signif TADs
      pval_thresh=pvalThresh_seq[1]
      adjCombPval_signifTADs <- lapply(pvalThresh_seq, function(pval_thresh) {
        names(adj_empPval_comb)[adj_empPval_comb <= pval_thresh]
      })
      names(adjCombPval_signifTADs) <- paste0(pvalThresh_seq)
      
      # RETRIEVE FDR DATA FOR LOGFC
      # for each of the FDR threshold -> get FC cut-off and meanCorr cut-off => signif TADs those with abs(logFC) >= cut-off & meanCorr >= cut-off
      logFC_FDR_file <- file.path(pipOutFolder, hicds, exprds, script19_name, "empFDR_list.Rdata")
      stopifnot(file.exists(logFC_FDR_file))
      all_FDR <- eval(parse(text = load(logFC_FDR_file)))
      logFC_FDR <- all_FDR[["empFDR_logFC"]]  # the names is the FC threshold, the value is the FDR
      stopifnot(length(logFC_FDR) > 0)
      
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
      FDR_signifTADs <- lapply(FDRthresh_seq , function(cutoff_fdr) {
        logFC_cut_off <- min(as.numeric(as.character(na.omit(names(logFC_FDR)[logFC_FDR <= cutoff_fdr]))))  # the smallest FC cut-off that leads to desired FDR; if not returns Inf
        meanCorr_cut_off <- min(as.numeric(as.character(na.omit(names(meanCorr_FDR)[meanCorr_FDR <= cutoff_fdr]))))
        stopifnot(names(tad_logFC) == names(tad_meanCorr))
        names(tad_logFC)[ abs(tad_logFC) >= logFC_cut_off & tad_meanCorr >= meanCorr_cut_off]
      })
      names(FDR_signifTADs) <- paste0(FDRthresh_seq)
      
      fdr=FDRthresh_seq[1]
      # => ITERATE OVER VARIOUS FDR AND PVAL THRESH TO GET THE INTERSECT OF SIGNIF TADS
      all_intersectTADs <- foreach(fdr = FDRthresh_seq, .combine='rbind') %do% {
        stopifnot(paste0(fdr) %in% names(FDR_signifTADs))
        pval=pvalThresh_seq[1]
        fdr_tads <- FDR_signifTADs[[paste0(fdr)]]
        stopifnot(fdr_tads %in% all_regs)
        intersectTADs_pval_dt <- foreach(pval = pvalThresh_seq, .combine='rbind') %do% {
          stopifnot(paste0(pval) %in% names(adjCombPval_signifTADs))
          pval_tads <- adjCombPval_signifTADs[[paste0(pval)]]
          stopifnot(pval_tads %in% all_regs)
          intersect_tads <- intersect(pval_tads, fdr_tads)
          FDR_nSignifTADs <- length(fdr_tads)
          adjCombPval_nSignifTADs <- length(fdr_tads)
          intersect_nSignifTADs <- length(intersect_tads)
          data.frame(
            hicds = hicds,
            exprds = exprds,
            FDR_threshold = fdr,
            adjPval_threshold = pval,
            FDR_signifTADs = paste0(fdr_tads, collapse=","),
            adjCombPval_signifTADs = paste0(pval_tads, collapse=","),
            intersect_signifTADs = paste0(intersect_tads, collapse=","),
            FDR_nSignifTADs = length(fdr_tads),
            adjCombPval_nSignifTADs = length(pval_tads),
            intersect_nSignifTADs = length(intersect_tads),
            stringsAsFactors = FALSE
          )
        } # end-foreach iterating over pval cutoffs
        intersectTADs_pval_dt
      } # end-foreach iterating over FDR cutoffs
      all_intersectTADs
    } # end-foreach iterating over exprds
    exprds_intersect_DT
  }# end-foreach iterating over hicds => allDS_intersect_DT
  outFile <- file.path(outFolder, "allDS_intersect_DT.Rdata")  
  save(allDS_intersect_DT, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
} else { # end-if build data
  outFile <- file.path(outFolder, "allDS_intersect_DT.Rdata")  
  allDS_intersect_DT <- eval(parse(text = load(outFile)))
}
allDS_intersect_DT$hicds <- as.character(allDS_intersect_DT$hicds)
allDS_intersect_DT$exprds <- as.character(allDS_intersect_DT$exprds)
allDS_intersect_DT$dataset <- paste0(allDS_intersect_DT$hicds, " - ", allDS_intersect_DT$exprds)

# HEATMAP SEPARATELY FOR EACH DATASET
myxlab <- "FDR thresh."
myylab <- "adj. comb. pval. thresh."

if(separateHeatmap) {
  ds = unique(allDS_intersect_DT$dataset)[1]
  for(ds in unique(allDS_intersect_DT$dataset)) {
    intersectDT <- allDS_intersect_DT[allDS_intersect_DT$dataset == ds, c("FDR_threshold", "adjPval_threshold", "intersect_nSignifTADs", "FDR_nSignifTADs", "adjCombPval_nSignifTADs")]
    stopifnot(nrow(intersectDT) > 0)
    
    intersectDT$FDR_threshold_label <- paste0(intersectDT$FDR_threshold, "\n(n=", intersectDT$FDR_nSignifTADs, ")")
    intersectDT$adjPval_threshold_label <- paste0(intersectDT$adjPval_threshold, "\n(n=", intersectDT$adjCombPval_nSignifTADs, ")")
    
    
    ds_intersect_plot <- ggplot(intersectDT, aes(x=FDR_threshold_label, y=adjPval_threshold_label, fill = intersect_nSignifTADs))+
      geom_tile(color = "white")+
      geom_text(aes(label = intersect_nSignifTADs), color = "black", size = 6, fontface="bold") +
      scale_fill_gradientn(colours = hm.palette(100))+
      guides(fill = FALSE)+ 
      labs(title = paste0(gsub(" - ", "\n", ds)),
           subtitle = paste0("intersect signif. TADs"),
          x = paste0(myxlab),
           y = paste0(myylab)
           )+
      theme(
        axis.text = element_text(size=12),
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(size=18, face="bold", hjust=0.5),
        plot.subtitle = element_text(size=16, face="italic", hjust=0.5),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
      )
      # coord_equal() # make the x and y at the same scale => same as coord_fixed(ratio = 1)
    outFile <- file.path(outFolder, paste0( gsub(" ", "", gsub("-", "_", ds)), "_intersectSignifTADs_heatmap.", plotType))
    ggsave(ds_intersect_plot, filename = outFile, height = myHeightGG, width = myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))  
  }
}
# HEATMAP SEPARATELY FOR AVERAGE VALUES OVER DS
nDS <- length(unique(allDS_intersect_DT$dataset))
avgDT_dt <- allDS_intersect_DT[, c("FDR_threshold", "adjPval_threshold", "intersect_nSignifTADs", "FDR_nSignifTADs", "adjCombPval_nSignifTADs")]
avgDT <- aggregate(.~FDR_threshold + adjPval_threshold, mean,data=avgDT_dt)
avgDT <- round(avgDT, 2)
avgDT$FDR_threshold_label <- paste0(avgDT$FDR_threshold, "\n(n=", avgDT$FDR_nSignifTADs, ")")
avgDT$adjPval_threshold_label <- paste0(avgDT$adjPval_threshold, "\n(n=", avgDT$adjCombPval_nSignifTADs, ")")

avg_intersect_plot <- ggplot(avgDT, aes(x=FDR_threshold_label, y=adjPval_threshold_label, fill = intersect_nSignifTADs))+
  geom_tile(color = "white")+
  geom_text(aes(label = intersect_nSignifTADs), color = "black", size = 6, fontface="bold") +
  scale_fill_gradientn(colours = hm.palette(100))+
  guides(fill = FALSE)+ 
  labs(title = paste0("Avg. over all ds (n=", nDS, ")"),
       subtitle = paste0("intersect signif. TADs"),
       x = paste0(myxlab),
       y = paste0(myylab)
  )+
  theme(
    axis.text = element_text(size=12),
    axis.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=18, face="bold", hjust=0.5),
    plot.subtitle = element_text(size=16, face="italic", hjust=0.5),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

outFile <- file.path(outFolder, paste0("avgAllDS_intersectSignifTADs_heatmap.", plotType))
ggsave(avg_intersect_plot, filename = outFile, height = myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

nDS <- length(unique(allDS_intersect_DT$dataset))

# HEATMAP FOR POOLED INTERSECT => this is just the sum !!!
# all_fdr <- unique(allDS_intersect_DT$FDR_threshold)
# fdr = all_fdr[1]  
# pooled_dt <- foreach(fdr = all_fdr, .combine='rbind') %do% {
#   sub_fdr_dt <- allDS_intersect_DT[allDS_intersect_DT$FDR_threshold == fdr,]
#   stopifnot(nrow(sub_fdr_dt) > 0)
#   all_pval <- unique(sub_fdr_dt$adjPval_threshold)
#   pval=all_pval[1]
#   pval_dt <- foreach(pval = all_pval, .combine='rbind') %do% {
#     sub_dt <- sub_fdr_dt[sub_fdr_dt$adjPval_threshold ==pval,]
#     stopifnot(nrow(sub_dt) > 0)
#     all_fdr_tads <- unlist(sapply(1:nrow(sub_dt), function(i) {
#               sub_dt
#               hic=sub_dt$hicds[i]
#               expr=sub_dt$exprds[i]
#               tads=sub_dt$FDR_signifTADs[i]
#               as.character(sapply(unlist(strsplit(tads, ",")), function(x) paste(hic, expr, x,sep="-")))  
#     }))
#     all_pval_tads <- unlist(sapply(1:nrow(sub_dt), function(i) {
#       sub_dt
#       hic=sub_dt$hicds[i]
#       expr=sub_dt$exprds[i]
#       tads=sub_dt$adjCombPval_signifTADs[i]
#       as.character(sapply(unlist(strsplit(tads, ",")), function(x) paste(hic, expr, x,sep="-")))  
#     }))
#     all_intersect_tads <- intersect(all_fdr_tads, all_pval_tads)
#     data.frame(
#       FDR_threshold = fdr,
#       adjPval_threshold = pval,
#       FDR_signifTADs = paste0(all_fdr_tads, collapse=","),
#       adjCombPval_signifTADs = paste0(all_pval_tads, collapse=","),
#       intersect_signifTADs = paste0(all_intersect_tads, collapse=","),
#       FDR_nSignifTADs = length(all_fdr_tads),
#       adjCombPval_nSignifTADs = length(all_pval_tads),
#       intersect_nSignifTADs = length(all_intersect_tads),
#       stringsAsFactors = FALSE
#     )
#   } # end-for each iterating over pval
#   pval_dt
# }
sumDT_dt <- allDS_intersect_DT[, c("FDR_threshold", "adjPval_threshold", "intersect_nSignifTADs", "FDR_nSignifTADs", "adjCombPval_nSignifTADs")]
sumDT <- aggregate(.~FDR_threshold + adjPval_threshold, sum,data=sumDT_dt)
# pooled_dt[pooled_dt$FDR_threshold==0.1 & 
#             pooled_dt$adjPval_threshold == 0.01, c("FDR_nSignifTADs", "adjCombPval_nSignifTADs")]
# sumDT[sumDT$FDR_threshold==0.1 & 
#         sumDT$adjPval_threshold == 0.01, c("FDR_nSignifTADs", "adjCombPval_nSignifTADs")]
sumDT <- round(sumDT, 2)
sumDT$FDR_threshold_label <- paste0(sumDT$FDR_threshold, "\n(n=", sumDT$FDR_nSignifTADs, ")")
sumDT$adjPval_threshold_label <- paste0(sumDT$adjPval_threshold, "\n(n=", sumDT$adjCombPval_nSignifTADs, ")")

avg_intersect_plot <- ggplot(sumDT, aes(x=FDR_threshold_label, y=adjPval_threshold_label, fill = intersect_nSignifTADs))+
  geom_tile(color = "white")+
  geom_text(aes(label = intersect_nSignifTADs), color = "black", size = 6, fontface="bold") +
  scale_fill_gradientn(colours = hm.palette(100))+
  guides(fill = FALSE)+ 
  labs(title = paste0("Sum over all ds (n=", nDS, ")"),
       subtitle = paste0("intersect signif. TADs"),
       x = paste0(myxlab),
       y = paste0(myylab)
  )+
  theme(
    axis.text = element_text(size=12),
    axis.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=18, face="bold", hjust=0.5),
    plot.subtitle = element_text(size=16, face="italic", hjust=0.5),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

outFile <- file.path(outFolder, paste0("sumAllDS_intersectSignifTADs_heatmap.", plotType))
ggsave(avg_intersect_plot, filename = outFile, height = myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))




# HEATMAP FOR THE AVERAGE ALL DS


avgDT_dt <- allDS_intersect_DT[, c("FDR_threshold", "adjPval_threshold", "intersect_nSignifTADs", "FDR_nSignifTADs", "adjCombPval_nSignifTADs")]
avgDT <- aggregate(.~FDR_threshold + adjPval_threshold, mean,data=avgDT_dt)
avgDT <- round(avgDT, 2)
avgDT$FDR_threshold_label <- paste0(avgDT$FDR_threshold, "\n(n=", avgDT$FDR_nSignifTADs, ")")
avgDT$adjPval_threshold_label <- paste0(avgDT$adjPval_threshold, "\n(n=", avgDT$adjCombPval_nSignifTADs, ")")

avg_intersect_plot <- ggplot(avgDT, aes(x=FDR_threshold_label, y=adjPval_threshold_label, fill = intersect_nSignifTADs))+
  geom_tile(color = "white")+
  geom_text(aes(label = intersect_nSignifTADs), color = "black", size = 6, fontface="bold") +
  scale_fill_gradientn(colours = hm.palette(100))+
  guides(fill = FALSE)+ 
  labs(title = paste0("Avg. over all ds (n=", nDS, ")"),
       subtitle = paste0("intersect signif. TADs"),
       x = paste0(myxlab),
       y = paste0(myylab)
  )+
  theme(
    axis.text = element_text(size=12),
    axis.title = element_text(size=14, face="bold"),
    plot.title = element_text(size=18, face="bold", hjust=0.5),
    plot.subtitle = element_text(size=16, face="italic", hjust=0.5),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

outFile <- file.path(outFolder, paste0("avgAllDS_intersectSignifTADs_heatmap.", plotType))
ggsave(avg_intersect_plot, filename = outFile, height = myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


########################################################################################## LATTICE PLOT

lattice_dt <- allDS_intersect_DT[, c("FDR_threshold", "adjPval_threshold", "intersect_nSignifTADs", "FDR_nSignifTADs", "adjCombPval_nSignifTADs")]

# lattice_dt$thresholdComb <- interaction(lattice_dt$FDR_threshold, lattice_dt$adjPval_threshold)
# lattice_dt$thresholdComb_label <- paste0("FDR thresh = ", lattice_dt$FDR_threshold, ";\nadj. comb. Pval. thresh = ", lattice_dt$adjPval_threshold)
lattice_dt$thresholdComb_label <- paste0("F=", lattice_dt$FDR_threshold, ";P=", lattice_dt$adjPval_threshold)
lattice_dt$thresholdComb <- interaction(lattice_dt$adjPval_threshold, lattice_dt$FDR_threshold)
lattice_dt$thresholdComb_label <- paste0("Pval<=", lattice_dt$adjPval_threshold, "; FDR<=", lattice_dt$FDR_threshold)


melt_lattice_dt <- melt(lattice_dt, id=c("FDR_threshold", "adjPval_threshold", "thresholdComb", "thresholdComb_label"))
colnames(melt_lattice_dt) [colnames(melt_lattice_dt) == "value"] <- "nSignifTADs"
melt_lattice_dt$variable <- gsub("_nSignifTADs", "", melt_lattice_dt$variable)


groupCols <- trellis.par.get("superpose.line")$col[1:length(unique(melt_lattice_dt$variable))]
names(groupCols) <- levels(as.factor(melt_lattice_dt$variable))

outFile <- file.path(outFolder, paste0("allDS_nSignifTADs_density_lattice.", plotType))
do.call(plotType, list(outFile, height=myHeight*3, width=myWidth*3))

densityplot(~nSignifTADs |thresholdComb_label, groups  = variable,data = melt_lattice_dt, #auto.key = TRUE, 
            # par.strip.text=list(cex=1), # width of the strip bar
            par.strip.text = list(cex = 1, font = 4, col = "brown"),
            layout = c(5, 5),
            scales=list(y=list(relation="free"),
                        x=list(relation="free")
                        ),
            panel = function(x,groups,subscripts, ...) {
              panel.densityplot(x,groups=groups, subscripts=subscripts)
              obs_nbr <-  as.numeric(table(groups[subscripts]))
              obs_cols <- groupCols[names(table(groups[subscripts]))]
              draw.key(key = list(text = list(
                paste0("n=", as.character(obs_nbr))), col= obs_cols), 
                draw=TRUE,  
                vp=viewport(x=0.9,y=0.9))
            },
            auto.key=list(title="", space = "bottom", cex=1.0, columns=length(unique(melt_lattice_dt$variable))),
            main = paste0("# signif. TADs"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))





                       


# which are the top ranking FCC AUC datasets
# load("../Yuanlong_Cancer_HiC_data_TAD_DA/EXPR_VARIANCE/LOG2FPKM/aucFCC.Rdata")
# head(sort(aucFCC, decreasing = TRUE),10)
# ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_wt_mutCTNNB1 ***
# 1.471168 
# ENCSR862OGI_RPMI-7951_40kb_TCGAskcm_wt_mutCTNNB1 
# 1.463047 
# ENCSR444WCZ_A549_40kb_TCGAluad_mutKRAS_mutEGFR ***
# 1.450056 
# ENCSR862OGI_RPMI-7951_40kb_TCGAskcm_wt_mutBRAF 
# 1.443286 
# ENCSR489OCU_NCI-H460_40kb_TCGAluad_mutKRAS_mutEGFR 
# 1.437669 
# ENCSR312KHQ_SK-MEL-5_40kb_TCGAskcm_wt_mutBRAF 
# 1.424054 
# ENCSR079VIJ_G401_40kb_TCGAkich_norm_kich ***
# 1.420483 
# ENCSR401TBQ_Caki2_40kb_TCGAkich_norm_kich 
# 1.417405 
# GSE75070_MCF-7_shNS_40kb_TCGAbrca_lum_bas 
# 1.403271 
# ENCSR444WCZ_A549_40kb_TCGAluad_wt_mutKRAS 
# 1.400265 


##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

