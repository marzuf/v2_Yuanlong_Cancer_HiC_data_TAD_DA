# Rscript tad_genes_coord.R

require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)

# chain <- import.chain("hg19ToHg18.over.chain")
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# tx_hg19 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
# tx_hg18 <- liftOver(tx_hg19, chain)

outFolder <- file.path("TAD_GENES_COORD")
dir.create(outFolder, recursive = TRUE)

buildTable <- TRUE

all_hicds <- c("ENCSR444WCZ_A549_40kb", "LG1_40kb", "LG2_40kb", "ENCSR489OCU_NCI-H460_40kb")
hicds="LG1_40kb"
gene_oi <- "AKR1C1"

runFolder <- file.path("..")

normQt <- 0.95
qvalThresh <- 0.01
pvalThresh <- 0.01

plotType <-  "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4

# my_cols <- setNames(pal_jama()(5)[c(3, 2)], c(check_exprds, "other"))
# my_cols <- setNames(pal_jama()(5)[c(3, 2)], c(check_exprds, "other"))

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

stopifnot(gene_oi %in% entrez2symb)

entrez_oi <- names(entrez2symb)[entrez2symb == gene_oi]

if(buildTable) {
  
    
    hicds_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
      
      # cat(paste0("... start\t", ds_name, "\t", hicds, "\n"))
      
      g2t_dt_file <- file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt")
      g2t_dt <- read.delim(g2t_dt_file, stringsAsFactors = FALSE, header=FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
      g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
      tad_g2t_dt <- g2t_dt[grepl("_TAD", g2t_dt$region),]
      stopifnot(!duplicated(g2t_dt$entrezID))
      entrez2tad <- setNames(g2t_dt$region, g2t_dt$entrezID)
      
      stopifnot(entrez_oi %in% names(entrez2tad))
      
      tad_oi <- entrez2tad[names(entrez2tad) == entrez_oi]
      
      
      tad_dt <- read.delim(file.path(runFolder, hicds, "genes2tad/all_assigned_regions.txt"), 
                           col.names=c("chromo", "region", "start", "end"), header=F, stringsAsFactors=FALSE)
      
      stopifnot(tad_oi %in% tad_dt$region)
      tadChr_oi <- tad_dt$chromo[tad_dt$region == tad_oi]
      hg19_tadStart_oi <- tad_dt$start[tad_dt$region == tad_oi]
      hg19_tadEnd_oi <- tad_dt$end[tad_dt$region == tad_oi]

      oi_g2t_dt <- g2t_dt[g2t_dt$region == tad_oi,]
      stopifnot(oi_g2t_dt$entrezID %in% gff_dt$entrezID)
      oi_gff_dt <- gff_dt[gff_dt$entrezID %in% oi_g2t_dt$entrezID,]
      
      oi_gff_dt$region <- tad_oi
      oi_gff_dt$regionStart <- hg19_tadStart_oi
      oi_gff_dt$regionEnd <- hg19_tadEnd_oi
      oi_gff_dt$assembly <- NULL
      oi_gff_dt$entrezID <- NULL
            
      
      outFile <- file.path(outFolder, paste0(gene_oi, "_",tad_oi, "_genes2tad_", hicds, "_dt.txt"))
      write.table(oi_gff_dt, file=outFile, sep="\t", quote=F)
      cat(paste0("... written: ", outFile, "\n"))
      
      
    }
  #   hicds_dt
  # }
} 
