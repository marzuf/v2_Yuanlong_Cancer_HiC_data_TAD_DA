# GO ANALYSIS ON THE GENES FROM SIGNIF TADs


# Rscript revision_GO_byCptmt_binary.R

printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}

setDir <- ""

script_name <- "revision_GO_byCptmt.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(clusterProfiler)
require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

my_box_theme <- theme(
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


buildTable <- TRUE

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 5)
myWidth <- 8
axisCex <- 1.4
myHeightGG <- 6
myWidthGG <- 7

### HARD-CODED - MAIN SETTINGS
# _final:  discussion 04.08.2020 Giovanni - keep aran CPE data, first vial only

purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- paste0("Aran - ", pm)

corMet <- "pearson"
transfExpr <- "log10"
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)


source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
# source("../MANUSCRIPT_FIGURES/settings.R")

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 5)
myWidth <- 8

script0_name <- "0_prepGeneData"

####################################################################### prepare compartment assignment

# inFolder <- "REVISION_CHANGES_CPTMTLABELS_ALLDS" ### >> contain only tads in the pipeline and is by hicds and exprds
# outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
# tmp_dt <- get(load(outFile))
# tmp_dt$region_ID <- file.path(tmp_dt$hicds, tmp_dt$region)

# inFolder <- "REVISION_CHANGES_CPTMTLABELS_ALLDS_ALLTADS" ## => contain all tads and is by hicds and exprds
# outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
# tad2cptmt_dt <- get(load(outFile))

# tmp2_dt <- merge(tmp_dt[,c("region_ID","tad_eightCptmtLab" )], 
#                  tad2cptmt_dt[,c("region_ID","tad_eightCptmtLab" )], by="region_ID", all=FALSE)
# stopifnot(tmp2_dt$tad_eightCptmtLab.x == tmp2_dt$tad_eightCptmtLab.y)
# stopifnot(nrow(tmp_dt) == nrow(tmp2_dt))

inFolder <- "REVISION_CHANGES_CPTMTLABELS_ALLDS" ### >> contain only tads in the pipeline and is by hicds and exprds
outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
tad2cptmt_dt <- get(load(outFile))
tad2cptmt_dt$region_ID <- file.path(tad2cptmt_dt$hicds, tad2cptmt_dt$exprds, tad2cptmt_dt$region)
stopifnot(!duplicated(tad2cptmt_dt$region_ID))
tad2cptmt <- setNames(tad2cptmt_dt$tad_eightCptmtLab, tad2cptmt_dt$region_ID)
####################################################################### FOR THE GO
source("enricher_settings.R")

padjVarGO <- "p.adjust" # p.adjust or qvalue ???

barplot_vars <- c("foldEnrichment", "geneRatio", "log10_pval")
barplot_vars_tit <- setNames(c("Fold enrichment", "Gene ratio", paste0("-log10(", padjVarGO,  ")")), barplot_vars)
barplot_vars <- c( "log10_pval")
barplot_vars_tit <- setNames(c(paste0("-log10(", padjVarGO,  ")")), barplot_vars)

plotMaxBars <- 10

par_bot <- 8

padjVarGO_plotThresh <- 0.05

# GO for BP nad MF [do not take c5_CC]
if(enricher_ontologyType == "BP" | enricher_ontologyType == "MF" | enricher_ontologyType == "BP_MF"){
  gmtFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2", paste0("c5.", tolower(enricher_ontologyType), ".v6.1.entrez.gmt"))
} else {
  stop(paste0(enricher_ontologyType, " is not a valid ontologyType\n"))
}
stopifnot(file.exists(gmtFile))
c5_msigdb <- read.gmt(gmtFile)


outFolder <- file.path("REVISION_GO_BYCPTMT_BINARY")
dir.create(outFolder, recursive = TRUE)

runFolder <- "."

logFile <- file.path(outFolder, "enrichr_logfile.txt")

txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... gmtFile:\t", gmtFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher background genes:\t", "universe", "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pvalueCutoff:\t", enricher_pvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pAdjustMethod:\t", enricher_pAdjustMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher minGSSize:\t", enricher_minGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher maxGSSize:\t", enricher_maxGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher qvalueCutoff:\t", enricher_qvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher ontologyType:\t", enricher_ontologyType, "\n")
printAndLog(txt, logFile)
# txt <- paste0("... pval thresh select GOs:\t", pvalSelectGO, "\n")
# printAndLog(txt, logFile)
txt <- paste0("... enricher_results_sortGOby:\t", enricher_results_sortGOby, "\n")
printAndLog(txt, logFile)

txt <- paste0("... padjVarGO_plotThresh:\t", padjVarGO_plotThresh, "\n")
printAndLog(txt, logFile)
txt <- paste0("... plotMaxBars:\t", plotMaxBars, "\n")
printAndLog(txt, logFile)





setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


purity_file <- file.path(runFolder,"ALLTADS_AND_PURITY_FINAL", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt$region_id <- file.path(merge_dt$dataset, merge_dt$region)
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

purity_flagged_tads <- merge_dt$region_id[merge_dt$purityCorr <= purityCorrThresh]

# take those that have puritycorr annot + not flagged
# tads_signif_notPurityFlagged <- merge_dt$region_id[ (! merge_dt$region_id %in% purity_flagged_tads) &
#                                                       merge_dt$adjPvalComb <= signifThresh]
# FOR CPTMT -> DO NOT SELECT BY SIGNIF.
tads_signif_notPurityFlagged <- merge_dt$region_id[ (! merge_dt$region_id %in% purity_flagged_tads) ]
stopifnot(!purity_flagged_tads %in% tads_signif_notPurityFlagged)

# discard the flagged -> means that I keep those for which no puritycorr annot
merge_dt_all <- merge(agg_purity, resultData, by=c("dataset", "region"), all=TRUE)  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt_all$region_id <- file.path(merge_dt_all$dataset, merge_dt_all$region)
stopifnot(purity_flagged_tads %in%  merge_dt_all$region_id)
# tads_signif_notPurityFlaggedAll <- merge_dt_all$region_id[ (! merge_dt_all$region_id %in% purity_flagged_tads ) &
#                                                              merge_dt_all$adjPvalComb <= signifThresh]
# FOR CPTMT -> DO NOT SELECT BY SIGNIF.
tads_signif_notPurityFlaggedAll <- merge_dt_all$region_id[ (! merge_dt_all$region_id %in% purity_flagged_tads )]

stopifnot(tads_signif_notPurityFlagged %in% tads_signif_notPurityFlaggedAll)
stopifnot(! purity_flagged_tads  %in% tads_signif_notPurityFlaggedAll)
stopifnot(! purity_flagged_tads %in% tads_signif_notPurityFlagged)


all_data_list <- list(
  list(tad_list = tads_signif_notPurityFlagged, data_dt = merge_dt, name="keep not PF"),
  list(tad_list = tads_signif_notPurityFlaggedAll, data_dt = merge_dt_all, name = "discard PF")
)

save(all_data_list, file="all_data_list.Rdata", version=2)
all_data_list <- get(load("all_data_list.Rdata"))

tad2cptmt_dt$tad_eightCptmtLab <- substr(tad2cptmt_dt$tad_eightCptmtLab, 1, 1)
all_cptmts <- unique(tad2cptmt_dt$tad_eightCptmtLab)
cptmt = all_cptmts[1]

# rm("tad_list")
# rm("data_dt")

for(i in seq_along(all_data_list)) {
  
  init_tad_list <- all_data_list[[i]][["tad_list"]]
  init_data_dt <- all_data_list[[i]][["data_dt"]]
  myname <- all_data_list[[i]][["name"]]
  
  for(cptmt in all_cptmts) {
    
    data_dt <- init_data_dt
    tad_list <- init_tad_list
    
    cptmttads_to_keep <- tad2cptmt_dt$region_ID[tad2cptmt_dt$tad_eightCptmtLab == cptmt]
    nKeptCptmt <- length(cptmttads_to_keep )
    stopifnot(length(cptmttads_to_keep) > 0)
    stopifnot(cptmttads_to_keep %in% names(tad2cptmt))
    ######## ADDED HERE - ITERATE OVER COMPARTMENT
    stopifnot(data_dt$region_id %in% names(tad2cptmt))
    save(data_dt, file="data_dt.Rdata", version=2);
    # stop("--ok\n")
    data_dt <- data_dt[data_dt$region_id %in% cptmttads_to_keep,]
    
    tad_list <- tad_list[tad_list %in% cptmttads_to_keep]
    stopifnot(length(tad_list) > 0)
    
    stopifnot(tad_list %in% data_dt$region_id)
    signif_notPF_genes <- unique(unlist(strsplit(as.character(data_dt$region_genes[data_dt$region_id %in% tad_list]), split=",")))
    stopifnot(signif_notPF_genes %in% gff_dt$symbol)
    signif_notPF_entrez <- gff_dt$entrezID[gff_dt$symbol %in% signif_notPF_genes]
    
    # CHANGED HERE - WHEN DOING BY CPTMT -> UNIVERSE IS ALL GENES
    univers_genes <- unique(unlist(strsplit(as.character(init_data_dt$region_genes), split=",")))
    stopifnot(univers_genes %in% gff_dt$symbol)
    univers_entrez <- gff_dt$entrezID[gff_dt$symbol %in% univers_genes]
    
    
    go_univers_entrez <- as.character(univers_entrez)[as.character(univers_entrez) %in% as.character(c5_msigdb$gene)]
    stopifnot(length(go_univers_entrez)  > 0)
    go_signif_notPF_entrez <- as.character(signif_notPF_entrez)[as.character(signif_notPF_entrez) %in% as.character(c5_msigdb$gene)]
    stopifnot(length(go_signif_notPF_entrez)  > 0)
    
    nDS <- length(unique(data_dt$dataset))
    nTADs <- length(unlist(tad_list))
    nSignifGenes <- length(signif_notPF_entrez)
    nSignifGenesGO <- length(go_signif_notPF_entrez)
    nAllGenes <- length(univers_entrez)
    nAllGenesGO <- length(go_univers_entrez)
    
    
    cat("... available annot. for univers_entrez:\t", length(go_univers_entrez) , "/", length(univers_entrez), "\n")
    cat("... available annot. for univers_entrez:\t", length(go_signif_notPF_entrez) , "/", length(signif_notPF_entrez), "\n")
    
    
    stopifnot(go_signif_notPF_entrez %in% go_univers_entrez)
    
    cat(paste0(">  start enricher  \n"))
    
    if(length(go_signif_notPF_entrez) > 0) {
      
      entrez_signif_enrich <- enricher(gene = go_signif_notPF_entrez, 
                                       TERM2GENE=c5_msigdb,
                                       universe = go_univers_entrez,
                                       pvalueCutoff = enricher_pvalueCutoff, 
                                       pAdjustMethod = enricher_pAdjustMethod, 
                                       minGSSize = enricher_minGSSize, 
                                       maxGSSize = enricher_maxGSSize, 
                                       qvalueCutoff =enricher_qvalueCutoff)
      
      outfile <- file.path(outFolder, paste0(gsub("\\.", "", cptmt), "_", gsub(" ", "", myname), "_entrez_signif_enrich.Rdata"))
      save(entrez_signif_enrich, version=2, file=outfile)
      cat(paste0("... written: ", outfile, "\n"))
      # stop("-ok \n")
      
      entrez_signif_enrich_resultDT <- entrez_signif_enrich@result
      entrez_signif_enrich_resultDT <- entrez_signif_enrich_resultDT[order(entrez_signif_enrich_resultDT[,enricher_results_sortGOby], 
                                                                           decreasing=FALSE), ]
      entrez_signif_enrich_resultDT$log10_pval <- -log10(entrez_signif_enrich_resultDT[,paste0(padjVarGO)])
      entrez_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(entrez_signif_enrich_resultDT$GeneRatio, function(x) {
        gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
      entrez_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(entrez_signif_enrich_resultDT$BgRatio, function(x) {
        gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
      entrez_signif_enrich_resultDT$foldEnrichment <- entrez_signif_enrich_resultDT$geneRatio/entrez_signif_enrich_resultDT$bgRatio
      
      
      genes_signif_plotMax <- min(c(plotMaxBars, nrow(entrez_signif_enrich_resultDT)))
      
      
      if(genes_signif_plotMax > 0) {
        for(var_plot in barplot_vars) {
          
          entrez_signif_enrich_resultDT <- entrez_signif_enrich_resultDT[order(entrez_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
          
          my_x <- gsub("GO_", "", entrez_signif_enrich_resultDT$Description[1:genes_signif_plotMax])
          
          plot_dt <- entrez_signif_enrich_resultDT[1:genes_signif_plotMax,c(var_plot, "Description"), drop=FALSE]
          plot_dt$labs <-  gsub("GO_", "", plot_dt$Description)
          
          strwdth <- 25
          
          plot_dt$plot_labs <-  unlist(lapply(strwrap(gsub("_", " ", plot_dt$labs), 
                                                      width = strwdth, simplify=FALSE), 
                                              function(x) paste0(x, collapse="\n")))
          
          plot_dt <- plot_dt[order(plot_dt[,paste0(var_plot)], decreasing=TRUE),]
          plot_dt$plot_labs <- factor(plot_dt$plot_labs, levels=as.character(plot_dt$plot_labs))
          
          
          myTit <- paste0(cptmt, " - Signif. enriched GO - ", myname, " (top ",var_plot,  ")")
          
          
          subTit <- paste0("(# DS=", nDS, " - # TADs=", nTADs, " - # genes=",
                           nSignifGenes, "(", nSignifGenesGO, ")/", nAllGenes, "(", nAllGenesGO, "))")
          
          
          if(var_plot == "log10_pval") {
            my_ylab <- "adj. p-val [-log10]"
          } else {
            my_ylab <- var_plot
          }
          
          ggbar_p <-  ggbarplot(plot_dt, 
                                x="plot_labs", y=var_plot, fill="darkgrey") +
            ggtitle(myTit, subtitle=subTit)+
            my_box_theme +
            labs(x="" , y = my_ylab) + 
            theme(
              axis.text.x = element_text(hjust=1, vjust=0.5,size=10,angle=90)
            )
          
          outFile <- file.path(outFolder,paste0(gsub("\\.", "", cptmt), "_",gsub(" ", "", myname), "_genesFromDAsignif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplotGG", ".", plotType))
          ggsave(ggbar_p, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
          cat(paste0("... written: ", outFile, "\n"))
          
        } # iterating compartment
          
        } # iterating over plot_var
      } # if GO signif
    } # if not PF genes
  } # notPF and discard PF
  
    
    





