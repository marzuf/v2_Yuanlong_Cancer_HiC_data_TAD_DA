# Rscript prep_targetGenes.R

outFolder <- "PREP_TARGETGENES"
dir.create(outFolder, recursive = TRUE)

# ensg_entrez_dt <- read.delim("../ensembl_entrez_biomart.txt", stringsAsFactors = FALSE, sep=",")
# ensg_entrez_dt <- unique(ensg_entrez_dt[,c("Gene.stable.ID", "EntrezGene.ID")])
# ensg_entrez_dt <- na.omit(ensg_entrez_dt)
# mean(ep_dt$target_ensembl %in% ensg_entrez_dt$Gene.stable.ID & ep_dt$target_symbol %in% gff_dt$symbol)
# 0.7500455

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$symbol))
symbol2entrez <- setNames(gff_dt$entrezID, gff_dt$symbol)

ep_dt <- read.delim("enhanceratlas_AllEPs_Lung_EP.txt", sep="$", stringsAsFactors = FALSE)
ep_dt <- ep_dt[,1:2]
colnames(ep_dt) <- c("enhancer", "target_symbol")
ep_dt$enhancer_chromo <- gsub("(.+):.+-.+_ENSG.+","\\1", ep_dt$enhancer )
ep_dt$enhancer_start <- as.numeric(gsub(".+:(.+)-.+_ENSG.+","\\1", ep_dt$enhancer ))
stopifnot(!is.na(ep_dt$enhancer_start))
ep_dt$enhancer_end <- gsub(".+:.+-(.+)_ENSG.+","\\1", ep_dt$enhancer )
stopifnot(!is.na(ep_dt$enhancer_end))
ep_dt$target_ensembl <- gsub(".+:.+-.+_(ENSG.+)","\\1", ep_dt$enhancer )

cat(paste0(nrow(ep_dt), " \n"))

ep_dt <- ep_dt[ep_dt$target_symbol %in% gff_dt$symbol,]

cat(paste0(nrow(ep_dt), " \n"))

# mean(ep_dt$target_symbol %in% gff_dt$symbol)
# 0.7783109

ep_dt$target_entrezID <- symbol2entrez[paste0(ep_dt$target_symbol)]
stopifnot(!is.na(ep_dt$target_entrezID))
ep_dt$target_entrezID <- as.character(ep_dt$target_entrezID)
ep_dt <- ep_dt[, c("enhancer_chromo", "enhancer_start", "enhancer_end", "target_entrezID")]

enhancerEntrezTarget_dt <- ep_dt

outFile <- file.path(outFolder, "enhancerEntrezTarget_dt.Rdata")
save(enhancerEntrezTarget_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

