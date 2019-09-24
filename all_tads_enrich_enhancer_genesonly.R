startTime <- Sys.time()
cat(paste0("> Rscript all_tads_enrich_enhancer_genesonly.R\n"))

# Rscript all_tads_enrich_enhancer_genesonly.R

options(scipen=100)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(40)
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

buildTable <- FALSE



plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightGG <- 7
myWidthGG <- myHeightGG*1.2


col1 <- "#00AFBB"
col2 <- "#FC4E07"

outFolder <- file.path("ALL_TADS_ENRICH_ENHANCER_ANNOT_GENESONLY")
dir.create(outFolder, recursive = TRUE)
outfile <- file.path(outFolder, "all_tads_enrich_enhancer_annot_genesonly.txt")
file.remove(outfile)

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))

inFolder <- "ENRICH_SAME_DIFF_TAD"
inFile <- file.path(inFolder, "full_enhancerDT.Rdata")
stopifnot(file.exists(inFile))
full_enhancer_DT <- get(load(inFile))
full_enhancer_DT$entrezID <- as.character(full_enhancer_DT$entrezID)

pipFolder <- file.path("PIPELINE","OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))

setDir="/media/electron"
setDir=""

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
stopifnot(!duplicated(entrez2symb_dt$entrezID))
stopifnot(!duplicated(entrez2symb_dt$symbol))
entrez2symb <- setNames(entrez2symb_dt$symbol, entrez2symb_dt$entrezID)


final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
stopifnot(signif_column %in% colnames(final_table_DT))
final_table_DT[,paste0(signifcol)] <- final_table_DT[,paste0(signif_column)] <= signifThresh

signif_final_table_DT <- final_table_DT[,c("hicds", "exprds", "region", "signifFDR_0.2", signifcol)]

all_hicds <- list.files(pipFolder)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

# all_hicds=all_hicds[1]

if(buildTable) {
  
  hicds="GSE105381_HepG2_40kb"
  all_enh_tad_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds="TCGAlihc_norm_lihc"  
    ds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat("... start ", hicds, " - ", exprds, "\n")
      
      
      g2tFile <- file.path(hicds, "genes2tad", "all_genes_positions.txt")
      stopifnot(file.exists(g2tFile))
      g2t_DT <- read.delim(g2tFile, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_DT$entrezID <- as.character(g2t_DT$entrezID)
      
      
      tadposFile <- file.path(hicds, "genes2tad", "all_assigned_regions.txt")
      stopifnot(file.exists(tadposFile))
      tadpos_DT <- read.delim(tadposFile, header=F, col.names = c( "chromo", "region", "start", "end"), stringsAsFactors = FALSE)
      tadpos_DT <- tadpos_DT[grepl("_TAD", tadpos_DT$region),]
      
      geneList_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      geneList <- eval(parse(text = load(geneList_file)))
      
      stopifnot(geneList %in% g2t_DT$entrezID)
      
      g2t_DT <- g2t_DT[g2t_DT$entrezID %in% geneList,]
      stopifnot(g2t_DT$entrezID %in% entrez2symb_dt$entrezID)
      g2t_DT$symbol <- sapply(g2t_DT$entrezID, function(x) entrez2symb_dt$symbol[entrez2symb_dt$entrezID == x])
      
      
      all_regs <- unique(as.character(g2t_DT$region))
      

            reg = all_regs[1]
      all_regs_dt <- foreach(reg = all_regs, .combine='rbind') %dopar% {
        
        
        reg_start <- tadpos_DT$start[tadpos_DT$region == reg]
        stopifnot(length(reg_start) == 1)
        stopifnot(is.numeric(reg_start) )
        reg_end <- tadpos_DT$end[tadpos_DT$region == reg]
        stopifnot(length(reg_end) == 1)
        stopifnot(is.numeric(reg_end) )
        stopifnot(reg_end > reg_start)
        
        
        curr_entrez <- g2t_DT$entrezID[g2t_DT$region == reg]
        stopifnot(length(curr_entrez) > 0)
        stopifnot(curr_entrez %in% names(entrez2symb))
        
        curr_symbols <- entrez2symb[paste0(curr_entrez)]
        stopifnot(length(curr_entrez) == length(curr_symbols))
        
        
        
        matching_gene_enhancer_dt <- full_enhancer_DT[full_enhancer_DT$entrezID %in% curr_entrez | 
                                                        full_enhancer_DT$geneSymbol %in% curr_symbols, ]
        
        
        
        matching_gene_enhancer_dt$enhancer_in_tad <- matching_gene_enhancer_dt$enhancer_start >= reg_start &  matching_gene_enhancer_dt$enhancer_start <= reg_end
        
        if(nrow(matching_gene_enhancer_dt) > 0) {
          
          inTAD_dt <- aggregate(enhancer_in_tad ~ geneSymbol, data=matching_gene_enhancer_dt, FUN=sum)
          colnames(inTAD_dt)[colnames(inTAD_dt) == "enhancer_in_tad"] <- "nEnhancersSameTAD"
          
          outTAD_dt <- aggregate(enhancer_in_tad ~ geneSymbol, data=matching_gene_enhancer_dt, FUN=function(x) sum(!x))
          colnames(outTAD_dt)[colnames(outTAD_dt) == "enhancer_in_tad"] <- "nEnhancersDiffTAD"
          
          stopifnot(sum(inTAD_dt$nEnhancersSameTAD) + sum(outTAD_dt$nEnhancersDiffTAD) == nrow(matching_gene_enhancer_dt))
          
          
          out_dt_tmp <- merge(inTAD_dt, outTAD_dt, by = "geneSymbol")  
        } else {
          out_dt_tmp <-  data.frame(geneSymbol = c(),
                                    nEnhancersDiffTAD = c(),
                                    nEnhancersSameTAD = c(),
                                    stringsAsFactors = FALSE)
        }
        
        
        
        missingSymb <- curr_symbols[ ! (curr_symbols %in% out_dt_tmp$geneSymbol)]
        
        if(length(missingSymb) > 0 ) {
          missingDT <- data.frame(geneSymbol = missingSymb,
                                  nEnhancersDiffTAD = NA,
                                  nEnhancersSameTAD = NA,
                                  stringsAsFactors = FALSE
          )
          out_dt <- rbind(out_dt_tmp, missingDT)
          
        } else {
          out_dt <- out_dt_tmp
        }
        
        out_dt <- out_dt[order(out_dt$nEnhancersSameTAD, out_dt$nEnhancersDiffTAD, decreasing = TRUE),]
        
        all_cols <- colnames(out_dt)
        
        out_dt$region <- reg
        out_dt$hicds <- hicds
        out_dt$exprds <- exprds
        
        out_dt[,c("hicds", "exprds", "region", all_cols)]

        
        
      } # end-foreach iterating over regions
      
      all_regs_dt
      
    } # end-foreach exprds
    ds_dt
  } # end-foreach hicds
  
  outFile <- file.path(outFolder, "all_enh_tad_dt.Rdata")
  save(all_enh_tad_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_enh_tad_dt.Rdata")
  all_enh_tad_dt <- eval(parse(text = load(outFile)))
  
}


idcols <- c("hicds", "exprds", "region")

all_enh_signif_dt <- merge(all_enh_tad_dt, signif_final_table_DT, by=idcols, all.x=TRUE)

signif_cols <- c( "signifFDR_0.2", signifcol)
signcol=signif_cols[1]

for(signcol in signif_cols) {
  
  sub_dt <- all_enh_signif_dt[,c(idcols, signcol)]
  
  all_p <- c("nEnhancersSameTAD", "nEnhancersDiffTAD")
  
  for(to_p in all_p) {
    p <- ggdensity(all_enh_signif_dt, 
                    title = paste0("", signcol),
                    subtitle=paste0(""),
                    x = to_p,
                    xlab=paste0(to_p),
                    color = signcol, fill = signcol,
                    #add = "mean", rug = TRUE,
                    palette = c(col1, col2))
    
    outFile <- file.path(outFolder, paste0("all_ds_signif_",to_p, "_",  signcol, "_density.", plotType))
    ggsave(p, file = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile,"\n"))
    
  }
  
  
}





######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

