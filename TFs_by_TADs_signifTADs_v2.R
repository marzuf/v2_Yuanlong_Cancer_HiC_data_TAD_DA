startTime <- Sys.time()

cat(paste0("... start - ", startTime, "\n"))

require(foreach)
source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
require(doMC)
registerDoMC(40)

# Rscript TFs_by_TADs_signifTADs_v2.R crisp
# Rscript TFs_by_TADs_signifTADs_v2.R c3.mir
# Rscript TFs_by_TADs_signifTADs_v2.R c3.tft
# Rscript TFs_by_TADs_signifTADs_v2.R c3.all
# Rscript TFs_by_TADs_signifTADs_v2.R trrust
# Rscript TFs_by_TADs_signifTADs_v2.R tftg
# Rscript TFs_by_TADs_signifTADs_v2.R motifmap
# Rscript TFs_by_TADs_signifTADs_v2.R kegg
# Rscript TFs_by_TADs_signifTADs_v2.R chea3_all
# Rscript TFs_by_TADs_signifTADs_v2.R chea3_lung

# 


plotCex <- 1.4

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 400, 7)
myWidth <- ifelse(plotType == "png", 600, 9)
plotCex <- 1.4

nTop <- 10

fontFamily <- "Hershey"

require(ggsci)
top_col <- pal_d3()(2)[1]
last_col <- pal_d3()(2)[2]
# yarrr::transparent("grey", trans.val = .6)
mid_col <- "#BEBEBE66"

x_qt_val <- 0.2
y_qt_val <- 0.95


dsIn <- "crisp"
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1 | length(args) == 3)
dsIn <- args[1]
if(length(args) == 3) {
  all_hicds <- args[2]
  all_exprds <- args[3]

} else {
  all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
  all_hicds <- all_hicds[!grepl("_RANDOM", all_hicds)]
  all_hicds <- all_hicds[!grepl("_PERMUT", all_hicds)]

  all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))
}

stopifnot(dsIn %in% c("crisp", "c3.mir", "c3.all", "c3.tft", "trrust", "tftg", "motifmap", "kegg", "chea3_all", "chea3_lung"))

outFolder <- file.path(paste0("TFS_BY_TADS_SIGNIFTADS_V2_", toupper(dsIn)))
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

tad_signif_thresh <- 0.01


setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)

stopifnot(!duplicated(gff_dt$entrezID))
stopifnot(!duplicated(gff_dt$symbol))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
symb2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))


if(buildData){
  nRegFeat_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
  
    cat(paste0("> START - ", hicds,"\n"))
  
    if(dsIn == "crisp") {
      reg_file <- file.path("gene_set_library_crisp_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    } else if(dsIn == "chea3_all") {
      reg_file <- file.path("chea3_all_tissues_TFs_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    } else if(dsIn == "chea3_lung") {
      reg_file <- file.path("chea3_lung_TFs_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE)
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    } else if(dsIn == "trrust"){
      reg_file <- file.path("trrust_rawdata.human.tsv")
      reg_dt <- read.delim(reg_file, sep="\t", header=FALSE, stringsAsFactors = FALSE,
                           col.names = c("regSymbol", "targetSymbol", "direction", "ID"))
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt <- reg_dt[reg_dt$targetSymbol %in% names(symb2entrez),]
      cat(paste0("with Entrez: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
      reg_dt$targetEntrezID <- symb2entrez[reg_dt$targetSymbol]
      reg_dt$targetEntrezID <- as.character(reg_dt$targetEntrezID)
    } else if(dsIn == "tftg") {
      reg_file <- file.path("tftg_db_all_processed.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, 
                           col.names=c("regSymbol", "targetEntrezID"))
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    } else if(dsIn == "motifmap"){
      reg_file <- file.path("MOTIFMAP_ALLGENES/overlapDT_bp.Rdata")
      reg_dt <- get(load(reg_file))
      colnames(reg_dt)[colnames(reg_dt)=="entrezID"] <- "targetEntrezID"
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    } else if(dsIn == "kegg"){
      reg_file <- file.path("hsa_kegg_entrez.txt")
      reg_dt <- read.delim(reg_file, sep="\t", header=FALSE, stringsAsFactors = FALSE,
                           col.names = c("targetEntrezID", "regSymbol"))
      reg_dt$targetEntrezID <- gsub("hsa:", "",reg_dt$targetEntrezID )
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    }else {
      reg_file <- file.path(paste0(dsIn, ".v7.0.entrez_processed.txt"))
      reg_dt <- read.delim(reg_file, sep="\t", header=TRUE, stringsAsFactors = FALSE, 
                           col.names=c("regSymbol", "targetEntrezID"))
      cat(paste0("init nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    }
    
    g2t_file <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
    g2t_dt <- read.delim(g2t_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    
    g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
    stopifnot(!duplicated(g2t_dt$entrezID))
    g2t_vect <- setNames(g2t_dt$region, g2t_dt$entrezID)
    
    reg_dt <- reg_dt[reg_dt$targetEntrezID %in% g2t_dt$entrezID,]
    cat(paste0("with g2t assignment: nrow(reg_dt)", "\t=\t", nrow(reg_dt), "\n"))
    reg_dt$targetRegion <- g2t_vect[paste0(reg_dt$targetEntrezID)]
    stopifnot(!is.na(reg_dt))
      
    nbrReg_TADs_dt <- aggregate(regSymbol~targetRegion, data=reg_dt, function(x) length(unique(x)))
    
    hicds_reg_dt <- reg_dt  
    rm("reg_dt")

    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      
      if(dsIn == "chea3_lung") {
        if(! (grepl("lusc", exprds) | grepl("luad", exprds))) return(NULL)
      }
      
      plotTit <- paste0(hicds, "\n", exprds)
      
      result_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds, ]
      result_dt$tad_rank <- rank(result_dt$adjPvalComb, ties="min")
      result_dt$rev_tad_rank <- rank(-result_dt$adjPvalComb, ties="min")
      
      geneList <- get(load(file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata") ))
      stopifnot(geneList %in% g2t_dt$entrezID)
      gByTAD <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      
      
      nGbyT <- setNames(as.numeric(table(g2t_dt$region)), names(table(g2t_dt$region)))
      
      reg_dt <- hicds_reg_dt[hicds_reg_dt$targetEntrezID %in% geneList,]  # update 08.01.20 -> NEED ALSO TO SUBSET THE REGULATED FEATURES !
      
      # 1) # of genes in TAD
      tad_nGenes_dt <- aggregate(entrezID ~ region, data=gByTAD, FUN=function(x) length(x))
      colnames(tad_nGenes_dt)[colnames(tad_nGenes_dt) == "entrezID"] <- "nGenes"
      stopifnot(tad_nGenes_dt$nGenes >= 3)
      
      # 2) # of genes regulated within TAD
      tad_nRegGenes_dt <- aggregate(targetEntrezID~targetRegion, data=reg_dt, FUN=function(x)length(unique(x)) )
      colnames(tad_nRegGenes_dt)[colnames(tad_nRegGenes_dt) == "targetRegion"] <- "region"
      colnames(tad_nRegGenes_dt)[colnames(tad_nRegGenes_dt) == "targetEntrezID"] <- "nRegGenes"
      
      # 3) # of TFs within TAD
      tad_nTFs_dt <- aggregate(regSymbol~targetRegion, data=reg_dt, FUN=function(x)length(unique(x)) )
      colnames(tad_nTFs_dt)[colnames(tad_nTFs_dt) == "targetRegion"] <- "region"
      colnames(tad_nTFs_dt)[colnames(tad_nTFs_dt) == "regSymbol"] <- "nTFs"
      
      plot_dt <- merge(tad_nTFs_dt, merge(tad_nGenes_dt, tad_nRegGenes_dt,by="region"), by="region")
      
      stopifnot(plot_dt$nRegGenes <= plot_dt$nGenes)
      
      plot_dt$nTFs_byGenes <- plot_dt$nTFs/plot_dt$nGenes
      plot_dt$nRegGenes_byGenes <- plot_dt$nRegGenes/plot_dt$nGenes
      
      stopifnot(!duplicated(plot_dt$region))
      
      
      plot_dt$hicds <- hicds
      plot_dt$exprds <- exprds
      
      signif_plot_dt <- merge(plot_dt, result_dt[,c("hicds", "exprds", "region", "adjPvalComb")], by=c("hicds", "exprds", "region"))
      
      stopifnot(signif_plot_dt$region %in% names(nGbyT))
      signif_plot_dt$nGenes <- nGbyT[paste0(signif_plot_dt$region)]
      stopifnot(!is.na(signif_plot_dt$nGenes))
      
      
      signif_plot_dt$signif_lab <- ifelse(signif_plot_dt$adjPvalComb <= tad_signif_thresh, "signif.", "not signif.")
      
     
      nGenes_signif <- signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "signif."]
      nGenes_notSignif <- signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "not signif."]

      
      nTFs_signif <- signif_plot_dt$nTFs[signif_plot_dt$signif_lab == "signif."]
      nTFs_notSignif <- signif_plot_dt$nTFs[signif_plot_dt$signif_lab == "not signif."]
      
      nRegGenes_signif <- signif_plot_dt$nRegGenes[signif_plot_dt$signif_lab == "signif."]
      nRegGenes_notSignif <- signif_plot_dt$nRegGenes[signif_plot_dt$signif_lab == "not signif."]
      
      nTFsOVERnGenes_signif <- signif_plot_dt$nTFs[signif_plot_dt$signif_lab == "signif."]/signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "signif."]
      nTFsOVERnGenes_notSignif <- signif_plot_dt$nTFs[signif_plot_dt$signif_lab == "not signif."]/signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "not signif."]
      
      nRegGenesOVERnGenes_signif <- signif_plot_dt$nRegGenes[signif_plot_dt$signif_lab == "signif."]/signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "signif."]
      nRegGenesOVERnGenes_notSignif <- signif_plot_dt$nRegGenes[signif_plot_dt$signif_lab == "not signif."]/signif_plot_dt$nGenes[signif_plot_dt$signif_lab == "not signif."]
      
      if(hicds == "LG2_40kb" & exprds == "TCGAluad_nonsmoker_smoker") save(signif_plot_dt, file ="signif_plot_dt.Rdata", version=2)
      
      data.frame(
        hicds = hicds,
        exprds = exprds, 
        mean_nTFs_signif = mean(nTFs_signif),
        mean_nTFs_notSignif = mean(nTFs_notSignif),
        mean_nRegGenes_signif = mean(nRegGenes_signif),
        mean_nRegGenes_notSignif = mean(nRegGenes_notSignif),
        mean_nTFsOVERnGenes_signif = mean(nTFsOVERnGenes_signif),
        mean_nTFsOVERnGenes_notSignif = mean(nTFsOVERnGenes_notSignif),
        mean_nRegGenesOVERnGenes_signif = mean(nRegGenesOVERnGenes_signif),
        mean_nRegGenesOVERnGenes_notSignif = mean(nRegGenesOVERnGenes_notSignif),
        mean_nGenes_signif = mean(nGenes_signif),
        mean_nGenes_notSignif = mean(nGenes_notSignif),
        stringsAsFactors = FALSE
      )
    }# end-for iterating over exprds
    exprds_dt
  } # end-for iterating over hicds
  outFile <- file.path(outFolder, "nRegFeat_dt.Rdata")  
  save(nRegFeat_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
outFile <- file.path(outFolder, "nRegFeat_dt.Rdata")  
nRegFeat_dt <- get(load(outFile))
}  
# load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/nRegFeat_dt.Rdata")
outFile <- file.path(outFolder, paste0("nRegFeat_boxplot_allDS.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(mar=par()$mar+c(9,0,0,0))
boxplot(nRegFeat_dt[,!colnames(nRegFeat_dt) %in% c("hicds", "exprds")], las=2, main=paste0("all obs. ds (n=", length(unique(file.path(nRegFeat_dt$hicds, nRegFeat_dt$exprds))),")"),  cex.axis=0.8)
mtext(side=3, text = paste0(dsIn))
cat(paste0("... written: ", outFile, "\n"))

# load("TFS_BY_TADS_SIGNIFTADS_C3.TFT/nRegFeat_dt.Rdata")

keepCols <- c("mean_nTFs_signif", "mean_nTFs_notSignif", "mean_nGenes_signif", "mean_nGenes_notSignif", "mean_nTFsOVERnGenes_signif", "mean_nTFsOVERnGenes_notSignif")

outFile <- file.path(outFolder, paste0("nRegFeat_boxplot_allDS_keepCols.", plotType))  
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(mar=par()$mar+c(9,0,0,0))
boxplot(nRegFeat_dt[, keepCols], las=2, main=paste0("all obs. ds (n=", length(unique(file.path(nRegFeat_dt$hicds, nRegFeat_dt$exprds))),")"),  cex.axis=0.8)
mtext(side=3, text = paste0(dsIn))
cat(paste0("... written: ", outFile, "\n"))


# hicds                    exprds       region nTFs nGenes nRegGenes nTFs_byGenes nRegGenes_byGenes adjPvalComb
# 370 LG2_40kb TCGAluad_nonsmoker_smoker chr11_TAD468  409     23        13     29.21429         0.9285714   0.1671207
# 979 LG2_40kb TCGAluad_nonsmoker_smoker  chr19_TAD93  371     16         9     37.10000         0.9000000   0.1904504
# 807 LG2_40kb TCGAluad_nonsmoker_smoker  chr17_TAD38  366     19        13     26.14286         0.9285714   0.2993208
# yy=reg_dt[reg_dt$targetRegion == "chr11_TAD468,]
#  unique(yy$regSymbol)
#   [1] "TFAP2E"   "ASH1L"    "DOT1L"    "HMGA1"    "PRR12"    "SETBP1"   "SRCAP"    "ZBED1"    "ASCL2"    "BHLHA15" 
#  [11] "BHLHA9"   "HAND2"    "HES1"     "HES5"     "MLXIP"    "MNT"      "MYCN"     "NCOA2"    "SREBF1"   "SREBF2"  
#  [21] "TCF23"    "TCF4"     "TWIST2"   "POGK"     "ATF5"     "CEBPG"    "CREB1"    "CREB3"    "FOSL1"    "JDP2"    
#  [31] "MAF"      "MAFG"     "MAFK"     "NFE2L1"   "AEBP2"    "ANKZF1"   "CHAMP1"   "DPF1"     "E4F1"     "EEA1"    
#  [41] "FIZ1"     "GLIS2"    "GTF3A"    "GZF1"     "HIC2"     "IKZF2"    "KAT7"     "KIN"      "KLF16"    "MAZ"     
#  [51] "MYT1L"    "MZF1"     "OSR2"     "OVOL2"    "PEG3"     "PRDM15"   "PRDM2"    "REPIN1"   "SLC2A4RG" "SP1"     
#  [61] "SP3"      "SP5"      "VEZF1"    "WIZ"      "ZBTB21"   "ZBTB25"   "ZBTB3"    "ZBTB32"   "ZBTB40"   "ZBTB41"  
#  [71] "ZBTB43"   "ZBTB44"   "ZBTB45"   "ZBTB49"   "ZBTB6"    "ZBTB7A"   "ZBTB7B"   "ZBTB8A"   "ZBTB8B"   "ZFP14"   
#  [81] "ZFP28"    "ZFP3"     "ZFP41"    "ZFP82"    "ZFP90"    "ZFP92"    "ZFPM1"    "ZFX"      "ZFY"      "ZIK1"    
#  [91] "ZKSCAN3"  "ZKSCAN7"  "ZKSCAN8"  "ZMAT1"    "ZNF10"    "ZNF121"   "ZNF142"   "ZNF160"   "ZNF165"   "ZNF169"  
# [101] "ZNF174"   "ZNF180"   "ZNF189"   "ZNF197"   "ZNF200"   "ZNF208"   "ZNF212"   "ZNF213"   "ZNF22"    "ZNF225"  
# [111] "ZNF23"    "ZNF230"   "ZNF233"   "ZNF234"   "ZNF24"    "ZNF248"   "ZNF250"   "ZNF26"    "ZNF266"   "ZNF267"  
# [121] "ZNF273"   "ZNF28"    "ZNF282"   "ZNF283"   "ZNF32"    "ZNF320"   "ZNF322"   "ZNF326"   "ZNF329"   "ZNF33A"  
# [131] "ZNF33B"   "ZNF34"    "ZNF345"   "ZNF347"   "ZNF354C"  "ZNF37A"   "ZNF382"   "ZNF383"   "ZNF384"   "ZNF385C" 
# [141] "ZNF396"   "ZNF408"   "ZNF41"    "ZNF415"   "ZNF419"   "ZNF425"   "ZNF426"   "ZNF430"   "ZNF431"   "ZNF433"  
# [151] "ZNF436"   "ZNF439"   "ZNF445"   "ZNF446"   "ZNF45"    "ZNF451"   "ZNF48"    "ZNF483"   "ZNF486"   "ZNF488"  
# [161] "ZNF490"   "ZNF497"   "ZNF500"   "ZNF506"   "ZNF510"   "ZNF511"   "ZNF512B"  "ZNF517"   "ZNF521"   "ZNF527"  
# [171] "ZNF529"   "ZNF530"   "ZNF544"   "ZNF547"   "ZNF548"   "ZNF549"   "ZNF550"   "ZNF555"   "ZNF556"   "ZNF557"  
# [181] "ZNF561"   "ZNF562"   "ZNF565"   "ZNF566"   "ZNF568"   "ZNF57"    "ZNF576"   "ZNF579"   "ZNF580"   "ZNF583"  
# [191] "ZNF585A"  "ZNF585B"  "ZNF587"   "ZNF592"   "ZNF595"   "ZNF596"   "ZNF598"   "ZNF600"   "ZNF613"   "ZNF620"  
# [201] "ZNF628"   "ZNF639"   "ZNF641"   "ZNF644"   "ZNF655"   "ZNF665"   "ZNF667"   "ZNF668"   "ZNF674"   "ZNF675"  
# [211] "ZNF677"   "ZNF678"   "ZNF682"   "ZNF684"   "ZNF687"   "ZNF69"    "ZNF691"   "ZNF7"     "ZNF704"   "ZNF705A" 
# [221] "ZNF705E"  "ZNF706"   "ZNF713"   "ZNF714"   "ZNF716"   "ZNF718"   "ZNF726"   "ZNF727"   "ZNF729"   "ZNF735"  
# [231] "ZNF747"   "ZNF75D"   "ZNF761"   "ZNF765"   "ZNF766"   "ZNF768"   "ZNF771"   "ZNF772"   "ZNF777"   "ZNF783"  
# [241] "ZNF784"   "ZNF785"   "ZNF787"   "ZNF789"   "ZNF793"   "ZNF8"     "ZNF813"   "ZNF814"   "ZNF816"   "ZNF827"  
# [251] "ZNF829"   "ZNF83"    "ZNF835"   "ZNF84"    "ZNF843"   "ZNF844"   "ZNF846"   "ZNF85"    "ZNF850"   "ZNF865"  
# [261] "ZNF878"   "ZNF883"   "ZNF891"   "ZNF93"    "ZNF99"    "ZSCAN12"  "ZSCAN16"  "ZSCAN22"  "ZSCAN26"  "ZSCAN32" 
# [271] "ZSCAN5A"  "ZSCAN9"   "ZXDC"     "PATZ1"    "ZNF512"   "ZEB2"     "ZFHX4"    "MBNL2"    "CENPBD1"  "JRK"     
# [281] "TIGD3"    "LIN28A"   "YBX2"     "CUX1"     "ONECUT3"  "DNMT1"    "FBXL19"   "KMT2A"    "KMT2B"    "E2F1"    
# [291] "E2F4"     "E2F6"     "E2F7"     "TFDP1"    "TFDP2"    "ELF2"     "ELF4"     "ELK1"     "ELK4"     "ERF"     
# [301] "ETS1"     "ETV2"     "GABPA"    "SPDEF"    "SPIB"     "FLYWCH1"  "FOXK1"    "FOXK2"    "FOXM1"    "FOXN3"   
# [311] "FOXP4"    "FOXQ1"    "FOXS1"    "GATAD2A"  "TRPS1"    "TFCP2"    "GTF2IRD1" "CIC"      "HMGN3"    "SOX12"   
# [321] "SOX15"    "SOX6"     "TCF7L2"   "ALX3"     "CDX1"     "CRX"      "DLX4"     "HOXA10"   "HOXA3"    "HOXA6"   
# [331] "HOXA7"    "HOXA9"    "HOXB9"    "HOXC5"    "IRX1"     "LBX1"     "NANOG"    "PKNOX1"   "PKNOX2"   "SIX5"    
# [341] "VSX1"     "ZHX1"     "POU5F1B"  "POU5F2"   "POU6F1"   "HSF1"     "HSFY2"    "MSANTD1"  "MBD4"     "PIN1"    
# [351] "BAZ2A"    "MTERF2"   "MTERF3"   "DMTF1"    "MYB"      "MYBL2"    "MYPOP"    "MYSM1"    "SNAPC4"   "TERF1"   
# [361] "ESRRA"    "NR0B1"    "NR1D1"    "NR2C1"    "NR2E3"    "NR2F2"    "NR2F6"    "RARG"     "RORB"     "RXRA"    
# [371] "ARHGAP35" "CC2D1A"   "DRAP1"    "GPBP1L1"  "NME2"     "PA2G4"    "PCGF2"    "PREB"     "PURB"     "RBCK1"   
# [381] "SAFB2"    "SKIL"     "SKOR1"    "SNAPC5"   "SON"      "SPEN"     "THYN1"    "TMF1"     "TSC22D1"  "XPA"     
# [391] "TP53"     "PAX8"     "PAX9"     "NFATC3"   "RFX1"     "RFX5"     "DEAF1"    "GMEB2"    "NFIB"     "SMAD1"   
# [401] "SMAD3"    "SMAD4"    "STAT1"    "STAT2"    "TBX6"     "TBP"      "THAP4"    "THAP6"    "THAP7"   
# 
# > unique(yy$targetSymbol)
#  [1] "VPS11"   "DDX6"    "BCL9L"   "H2AFX"   "SLC37A4" "HYOU1"   "HINFP"   "DPAGT1"  "TRAPPC4" "HMBS"    "CCDC84"  "RPS25"  
# [13] "CXCR5" 

#####################################################################
cat("*** DONE\n")
cat(paste0("... end - ", Sys.time(), "\n"))
