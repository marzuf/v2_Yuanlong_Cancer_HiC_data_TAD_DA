# Rscript liftover_coord.R


require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)

# chain <- import.chain("hg19ToHg18.over.chain")
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# tx_hg19 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
# tx_hg18 <- liftOver(tx_hg19, chain)
require(liftOver)
require(GenomicRanges)
chain <- import.chain("hg19ToHg38.over.chain")

outFolder <- file.path("LIFTOVER_COORD")
dir.create(outFolder, recursive = TRUE)

buildTable <- TRUE

all_hicds <- c("ENCSR444WCZ_A549_40kb", "LG1_40kb", "LG2_40kb", "ENCSR489OCU_NCI-H460_40kb")
hicds="LG1_40kb"
gene_oi <- "AKR1C1"

all_files <- list.files("encode_data_38", pattern="\\.bed$", full.names=TRUE)

runFolder <- file.path("..", "..")

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

encode_annot_dt <- read.delim("metadata_filter.csv", header=FALSE, sep="\t", stringsAsFactors = FALSE)
bed2annot <- setNames(as.character(encode_annot_dt$V13), as.character(encode_annot_dt$V1))
bed2annot <- gsub("-human", "",  bed2annot)
bed2cl <- setNames(as.character(encode_annot_dt$V10), as.character(encode_annot_dt$V1))

cl_levels <- c("lung",
               "left lung",
               "upper lobe of left lung",
               "lower lobe of left lung",
               "AG04450", # fetal
               "IMR-90" , # normal
               "A549", # kras mut
               "PC-9") # egfr mut

stopifnot(gsub("\\.bed$", "", basename(basename(all_files))) %in% names(bed2annot))

filename = all_files[1]

if(buildTable) {
  
  # all_histOverlap_dt <- foreach(filename=all_files, .combine='rbind') %dopar% {
  #   
  #   ds_name <- gsub("\\.bed$", "", basename(basename(filename)))
  #   stopifnot(ds_name %in% names(bed2annot))
    # 
    # hist_dt <- read.delim(filename, header=F, stringsAsFactors = FALSE,
    #                       col.names=c("chromo", "chromStart", "chromEnd", "name",
    #                                   "score","strand","signalValue", "pValue_log10", "qValue_log10", "peak"))
    # #chr1    564540  564690  .       0       .       7       5.47141 -1      -1
    # #chr1    569820  569970  .       0       .       10      5.00553 -1      -1
    # hist_dt$signalValue_qtNorm <- hist_dt$signalValue/quantile(hist_dt$signalValue, probs = normQt)
    # hist_dt$signalValue_zNorm <- as.numeric(scale(hist_dt$signalValue, center = TRUE, scale = TRUE))
    # 
    
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
      
      tad_gr <- GRanges(seqnames=tadChr_oi, ranges=IRanges(start=hg19_tadStart_oi, end=hg19_tadEnd_oi))
      tad_hg38 <- liftOver(tad_gr, chain)
      tadStart_oi <- as.numeric(start(tad_hg38))
      tadEnd_oi <- as.numeric(end(tad_hg38))
      
      data.frame(
        hicds = hicds,
        gene_oi=gene_oi,
        tad_oi = tad_oi,
        hg19_tadStart_oi=hg19_tadStart_oi,
        hg19_tadEnd_oi=hg19_tadEnd_oi,
        tadStart_oi=tadStart_oi,
        tadEnd_oi=tadEnd_oi,
        stringsAsFactors = FALSE
      )
      
    }
  #   hicds_dt
  # }
  outFile <- file.path(outFolder, paste0(gene_oi, "_hicds_dt.txt"))
  write.table(hicds_dt, file=outFile, sep="\t", quote=F)
  cat(paste0("... written: ", outFile, "\n"))
} 
