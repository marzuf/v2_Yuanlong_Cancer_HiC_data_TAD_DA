# whant ot have data in the format
# enhancer_chromo enhancer_start enhancer_end target_entrezID
# 1            chr1         949920       950500           84069

# Rscript process_jeme.R

outFolder <- "JEME_DATA"
dir.create(outFolder, recursive = TRUE)

setDir <- "/media/electron"
setDir <- ""
entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb_dt$symbol <- as.character(entrez2symb_dt$symbol)
stopifnot(!duplicated(entrez2symb_dt$entrezID))
entrez2symb <- setNames(entrez2symb_dt$entrezID, entrez2symb_dt$symbol)

dt2 <- read.delim("fantom5_elasticnet.646_NCI-H460.csv", sep=",", header=F, col.names=c("enhancerID", "targetID", "score"), stringsAsFactors = FALSE)
nrow(dt2)
# 1448
dt2$targetSymbol <- sapply(dt2$targetID, function(x) unlist(strsplit(x, split="\\$"))[2])
dt2 <- dt2[dt2$targetSymbol %in% names(entrez2symb),]
nrow(dt2)
# 1427
dt2$target_entrezID <- entrez2symb[paste0(dt2$targetSymbol)]
stopifnot(!is.na(dt2$target_entrezID))
dt2$enhancer_start <- as.numeric(gsub(".+:(.+)-.+", "\\1", dt2$enhancerID))
stopifnot(!is.na(dt2$enhancer_start))
dt2$enhancer_end <- as.numeric(gsub(".+:.+-(.+)", "\\1", dt2$enhancerID))
stopifnot(!is.na(dt2$enhancer_end))
dt2$enhancer_chromo <- gsub("(.+):.+", "\\1", dt2$enhancerID)
stopifnot(grepl("chr", dt2$enhancer_chromo))

jeme_elasticnet_dt <- dt2[,c("enhancer_chromo", "enhancer_start", "enhancer_end", "target_entrezID")]
outFile <- file.path(outFolder, "jeme_elasticnet_dt.Rdata")
save(jeme_elasticnet_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


dt1 <- read.delim("fantom5_lasso.646_NCI-H460.csv", sep=",", header=F, col.names=c("enhancerID", "targetID", "score"), stringsAsFactors = FALSE)
nrow(dt1)
# 1448
dt1$targetSymbol <- sapply(dt1$targetID, function(x) unlist(strsplit(x, split="\\$"))[2])
dt1 <- dt1[dt1$targetSymbol %in% names(entrez2symb),]
nrow(dt1)
# 1427
dt1$target_entrezID <- entrez2symb[paste0(dt1$targetSymbol)]
stopifnot(!is.na(dt1$target_entrezID))
dt1$enhancer_start <- as.numeric(gsub(".+:(.+)-.+", "\\1", dt1$enhancerID))
stopifnot(!is.na(dt1$enhancer_start))
dt1$enhancer_end <- as.numeric(gsub(".+:.+-(.+)", "\\1", dt1$enhancerID))
stopifnot(!is.na(dt1$enhancer_end))
dt1$enhancer_chromo <- gsub("(.+):.+", "\\1", dt1$enhancerID)
stopifnot(grepl("chr", dt1$enhancer_chromo))

jeme_lasso_dt <- dt1[,c("enhancer_chromo", "enhancer_start", "enhancer_end", "target_entrezID")]
outFile <- file.path(outFolder, "jeme_lasso_dt.Rdata")
save(jeme_lasso_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))



