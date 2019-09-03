options(scipen=100)

setDir=""

# Rscript intersect_FDR_combPval_lolliPlot_conserv_annot.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript intersect_FDR_combPval_lolliPlot_conserv_annot.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich
# Rscript intersect_FDR_combPval_lolliPlot_conserv_annot.R   # to run all datasets in one shot

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"
hicds="ENCSR079VIJ_G401_40kb"
exprds="TCGAkich_norm_kich"


script_name <- "intersect_FDR_combPval_lolliPlot_conserv.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE

buildData <- TRUE
separateHeatmap <- TRUE

require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
require(reshape2)

# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# # source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
# source("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2/TAD_DE_utils.R")
# source("../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("plot_lolliTAD_funct.R")

script0_name <- "0_prepGeneData"
# script3_name <- "3_runMeanTADLogFC"
# script4_name <- "4_runMeanTADCorr"
# script9_name <- "9_runEmpPvalMeanTADLogFC"
# script19_name <- "19onlyFC_SAM_emp_measurement"
# script19sameNbr_name <- "19sameNbr_SAM_emp_measurement"
# script11same_name <- "11sameNbr_runEmpPvalCombined"

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

pipFolder <- file.path(".")
stopifnot(dir.exists(pipFolder))

pipOutFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "INTERSECT_FDR_COMBPVAL_LOLLIPLOT_CONSERV_ANNOT"
dir.create(outFolder, recursive=TRUE)


pvalThresh <- 0.01
matching_ratio_thresh <- 0.8

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

final_table_file <- file.path("CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

final_table_DT[, paste0("signifEmpPval_", pvalThresh)] <- final_table_DT$adjPvalComb <= pvalThresh

annot1_name <- "cancermine"
annot1_file <- file.path("db_gene_annot/cancermine_collated.tsv")
annot1_dt <- read.delim(annot1_file, header=TRUE, stringsAsFactors = FALSE)
annot1_dt$gene_entrez_id <- as.character(annot1_dt$gene_entrez_id)
ambiguous_entrez <- annot1_dt$gene_entrez_id[duplicated(annot1_dt$gene_entrez_id)]
annot1_dt <- annot1_dt[!annot1_dt$gene_entrez_id %in% ambiguous_entrez,]
stopifnot(!duplicated(annot1_dt$gene_entrez_id))
annot1_entrez <- setNames(annot1_dt$role, annot1_dt$gene_entrez_id)

annot2_name <- "cigene"
annot2_file <- file.path("db_gene_annot/cigene_human.txt")
annot2_dt <- read.delim(annot2_file, header=TRUE, stringsAsFactors = FALSE)
all_genes1 <- annot2_dt$GeneSymbol
all_genes2 <- unlist(strsplit(annot2_dt$Alias, split="\\|"))
annot2_symbols <- setNames(rep("cancer_initiation", length(c(all_genes1, all_genes2))), c(all_genes1, all_genes2))

annot3_name <- "cgi"
annot3_file <- file.path("db_gene_annot/cgi-bin/cgi_bin_gene_annot.txt")
annot3_dt <- read.delim(annot3_file, header=FALSE, stringsAsFactors = FALSE, col.names = c("symbol", "role", "description") )
stopifnot( !duplicated(annot3_dt$symbol))
annot3_symbols <- setNames(annot3_dt$role, annot3_dt$symbol)


entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
stopifnot(!duplicated(entrez2symb_dt$entrezID))
entrez2symb <- setNames(entrez2symb_dt$symbol, entrez2symb_dt$entrezID)


inFolder <- "TAD_MATCHING_ACROSS_HICDS"
inFile <- file.path(inFolder, "all_matching_dt.Rdata")
stopifnot(file.exists(inFile))

all_matching_dt <- get(load(inFile))
all_matching_dt$overlapBp_ratio <- all_matching_dt$overlapBP/all_matching_dt$ref_totBp
all_matching_dt$overlapGenes_ratio <- all_matching_dt$nOverlapGenes/all_matching_dt$ref_nGenes
stopifnot(na.omit(all_matching_dt$overlapBp_ratio) >= 0 & na.omit(all_matching_dt$overlapBp_ratio) <= 1)
stopifnot(na.omit(all_matching_dt$overlapGenes_ratio) >= 0 & na.omit(all_matching_dt$overlapGenes_ratio <= 1))


all_matching_dt$overlapBp_overMatchRatio <- all_matching_dt$overlapBp_ratio >= matching_ratio_thresh
all_matching_dt$overlapGenes_overMatchRatio <- all_matching_dt$overlapGenes_ratio >= matching_ratio_thresh

nDS <- length(unique(all_matching_dt$ref_dataset))

nOverThreshBp_dt <- aggregate(overlapBp_overMatchRatio ~ ref_dataset + refID, data = all_matching_dt, FUN=sum)
stopifnot(range(nOverThreshBp_dt$overlapBp_overMatchRatio) <= nDS)

nOverThreshGenes_dt <- aggregate(overlapGenes_overMatchRatio ~ ref_dataset + refID, data = all_matching_dt, FUN=sum)
stopifnot(range(nOverThreshGenes_dt$overlapGenes_overMatchRatio) <= nDS)

nOverThresh_dt <- merge(nOverThreshGenes_dt, nOverThreshBp_dt,by=c("ref_dataset", "refID"), all=TRUE)
nOverThresh_dt$hicds <- dirname(nOverThresh_dt$ref_dataset)
nOverThresh_dt$exprds <- basename(nOverThresh_dt$ref_dataset)
nOverThresh_dt$ref_dataset <- NULL

colnames(nOverThresh_dt)[colnames(nOverThresh_dt) == "refID"] <- "region"

sub_final_dt <- final_table_DT[,c("hicds", "exprds","region", "signifFDR_0.2", "signifEmpPval_0.01")]

plot_dt <- merge(sub_final_dt, nOverThresh_dt, by=c("hicds", "exprds", "region"))




############# > PLOT FOR EACH DATASET SEPARATELY


allDS_intersect_DT <- foreach(hicds = all_hicds) %dopar% {
  exprds_intersect_DT <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
    cat("... start ", hicds, " - " , exprds, "\n")
    ds_dt <- plot_dt[plot_dt$hicds == hicds & plot_dt$exprds == exprds,]
    stopifnot(nrow(ds_dt) > 0)
    signif_DT <- ds_dt[ds_dt$signifFDR_0.2 & ds_dt$signifEmpPval_0.01,]
    signif_DT <- signif_DT[order(signif_DT$overlapGenes_overMatchRatio, signif_DT$overlapBp_overMatchRatio, decreasing=TRUE),]
    
    if(nrow(signif_DT) == 0) return(NULL)
    
    
    stopifnot(dir.exists(file.path(pipOutFolder, hicds, exprds, script0_name)))
    geneListFile <- file.path(pipOutFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
    stopifnot(file.exists(geneListFile))
    pipeline_geneList <- eval(parse(text = load(geneListFile))) 
    
    g2tFile <- file.path(pipFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(g2tFile))
    g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
    g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
    stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
    g2t_DT <- g2t_DT[ g2t_DT$entrezID %in% pipeline_geneList,]
    
    # g2t_DT$annot <- annot_entrez[g2t_DT$entrezID]
    
    outfile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", "cancer_annot_signifTADs.txt"))
    file.remove(outfile)
    i=1
    
    for(i in 1:nrow(signif_DT)) {
      cat(i,"\n")
      curr_tad <- signif_DT$region[i]
      curr_entrez <- g2t_DT$entrezID[g2t_DT$region == curr_tad]
      stopifnot(!duplicated(curr_entrez))
      stopifnot(paste0(curr_entrez) %in% names(entrez2symb))
      curr_genes <- setNames(entrez2symb[paste0(curr_entrez)], paste0(curr_entrez))
      
      if(any(names(curr_genes) %in% names(annot1_entrez))) {
        
        a1_genes <- curr_genes[names(curr_genes) %in% names(annot1_entrez)]
        
        dt <- data.frame(entrezID = names(a1_genes),
                         symbol = as.character(a1_genes),
                         annot = as.character(annot1_entrez[names(a1_genes)]),
                        stringsAsFactors = FALSE)
        stopifnot(!is.na(dt))
        mycon <- file(outfile,"a")
        writeLines(text=paste0("\n> ", curr_tad, " - ", annot1_name, " annotation:"), con=mycon)
        close(mycon)
        write.table(dt, file = outfile, append=T, quote=F,col.names = FALSE,row.names=FALSE,sep="\t")
        
      }
      
      
      
      if(any(paste0(curr_genes) %in% names(annot2_symbols))) {
        
        a2_genes <- curr_genes[paste0(curr_genes) %in% names(annot2_symbols)]
        
        dt <- data.frame(entrezID = names(a2_genes),
                         symbol = as.character(a2_genes),
                         annot = as.character(annot2_symbols[paste0(a2_genes)]),
                         stringsAsFactors = FALSE)
        stopifnot(!is.na(dt))
        mycon <- file(outfile,"a")
        writeLines(text=paste0("\n> ", curr_tad, " - ", annot2_name, " annotation:"), con=mycon)
        close(mycon)
        write.table(dt, file = outfile, append=T, quote=F,col.names = FALSE,row.names=FALSE,sep="\t")
        
      }
      
      if(any(paste0(curr_genes) %in% names(annot3_symbols))) {
        
        a3_genes <- curr_genes[paste0(curr_genes) %in% names(annot3_symbols)]
        
        dt <- data.frame(entrezID = names(a3_genes),
                         symbol = as.character(a3_genes),
                         annot = as.character(annot3_symbols[paste0(a3_genes)]),
                         stringsAsFactors = FALSE)
        stopifnot(!is.na(dt))
        
        mycon <- file(outfile,"a")
        writeLines(text=paste0("\n> ", curr_tad, " - ",  annot3_name, " annotation:"), con=mycon)
        close(mycon)
        write.table(dt, file = outfile, append=T, quote=F,col.names = FALSE,row.names=FALSE,sep="\t")
        
      }
      
      
      
    }
    cat(paste0("... written:\t", outfile, "\n"))

      

    
    
    
  } # end-foreach iterating over exprds
}# end-foreach iterating over hicds => allDS_intersect_DT







##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

