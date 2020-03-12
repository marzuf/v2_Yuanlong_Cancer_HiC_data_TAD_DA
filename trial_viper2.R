library(aracne.networks)
data(regulonluad)
library(bcellViper)
library(viper)

# Rscript trial_viper2.R

startTime <- Sys.time()

outFolder <- file.path("TRIAL_VIPER2")
dir.create(outFolder, recursive = TRUE)

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

pipFolder <- "PIPELINE/OUTPUT_FOLDER"
settingFolder<- "PIPELINE/INPUT_FILES"

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"

settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)
samp1 <- get(load(file.path(setDir, sample1_file)))
samp2 <- get(load(file.path(setDir, sample2_file)))

exprData <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")))
exprData <- exprData[,c(samp1, samp2)]
exprData_norm <- apply(exprData, 2, function(x)x/sum(x))
exprData_norm <- exprData_norm[ ! rowSums(exprData_norm) == 0,]

log_offset <- 0.01

stopifnot( abs(colSums(exprData_norm)-1) <= 10^-4)
log2_fpkmDT <- log2(exprData_norm*10^6+log_offset)


stopifnot(samp1 %in% colnames(log2_fpkmDT))
stopifnot(samp2 %in% colnames(log2_fpkmDT))

cond1_exprDT <- log2_fpkmDT[,samp1]
cond2_exprDT <- log2_fpkmDT[,samp2]

cond1_exprDT <- as.matrix(cond1_exprDT)
cond2_exprDT <- as.matrix(cond2_exprDT)

cond1_exprDT[1:5,1:5]
cond2_exprDT[1:5,1:5]

luad_signature <- rowTtest(x=cond2_exprDT, y=cond1_exprDT, alternative="two.sided")
# x: ExpressionSet object or Numerical matrix containing the test samples
# y: Optional numerical matrix containing the reference samples.
luad_signature <- (qnorm(luad_signature$p.value/2, lower.tail = FALSE) * sign(luad_signature$statistic))[, 1]

luad_nullmodel <- ttestNull(x = cond1_exprDT, y = cond2_exprDT, per = 1000, repos = T, verbose = F)
# x: ExpressionSet object or Matrix containing the test dataset; 
# y=reference dataset
save(luad_nullmodel, file = file.path(outFolder, "luad_nullmodel.Rdata"), version=2)

luad_mrs <- msviper(luad_signature, regulonluad, luad_nullmodel, verbose = F)
save(luad_mrs, file = file.path(outFolder, "luad_mrs.Rdata"), version=2)

summary(luad_mrs)
# Normalized Enrichment Score (NES)

mrs_res <- data.frame(summary(luad_mrs,length(luad_mrs$es$nes)))
save(mrs_res, file = file.path(outFolder, "mrs_res.Rdata"), version=2)

mrs_res$symbol <- entrez2symb[rownames(mrs_res)]

write.table(mrs_res, file = file.path(outFolder, "mrs_res_table.txt"), sep="\t", quote=F, col.names=TRUE, row.names=TRUE)

topp = data.frame(summary(luad_mrs,length(luad_mrs$es$nes)))
topp$gene_symbol <- entrez2symb[rownames(mrs_res)]
### Plotting
maxx = 20
topp_plot = head(topp[order(topp$p.value),],maxx)
topp_plot = topp_plot[order(-topp_plot$NES),]
# topp_plot = topp_plot[order(-sign(topp_plot$NES),topp_plot$p.value),]
topMrs = rownames(topp_plot)
pdf(file = file.path(outFolder, paste0("viper_topPlot_top",maxx,".pdf")))
plot(mrs_res, mrs = topMrs, include = "")
axis(4, at = c(1:length(topMrs)), tick=F,labels = rev(topp[topMrs,"gene_symbol"]), cex.axis = 0.5, las = 1)
dev.off()

### Leading-edge analysis: which are the target genes of the master regulators
target_luad_mrs = ledge(luad_mrs)
summary(target_luad_mrs)


cat(paste0(startTime, "\n", Sys.time(), "\n"))

# 
# library(bcellViper)
# library(viper)
# data(bcellViper, package="bcellViper")
# adjfile <- system.file("aracne", "bcellaracne.adj", package = "bcellViper")
# signature <- rowTtest(dset, "description", c("CB", "CC"), "N")
# # It can also take two matrixes as arguments, the first one containing the ‘test’ samples and the second the
# # ‘reference’ samples.
# # While we could define the Gene Expression Signature (GES) by using the t-statistic, to be consistent
# # with the z-score based null model for msVIPER (see section 6.2), we will estimate z-score values for the
# # GES:
#  # signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
# 
# 
#  #
#  nullmodel <- ttestNull(dset, "description", c("CB", "CC"), "N", per = 1000,
#                           repos = TRUE, verbose = FALSE)
# 
#  vpres <- viper(dset, regulon, verbose = FALSE)
#  # The viper function generates a matrix – or ‘ExpressionSet’ object in case an ‘ExpressionSet’ object is
#  # given as input – of regulator’s activity, containing 621 regulators x 211 samples in our example.
#  dim(vpres)
#  # Features
#  # 621
#  # Samples
#  # 211
#  # The differential activity of regulatory proteins between groups of samples, for example between germinal
#  # center B-cell and Naı̈ve B-cells, can be obtained by any hypothesis testing statistical method, like for example
#  # the Student’s t-test:
# tmp <- rowTtest(vpres, "description", c("CB", "CC"), "N")
# 
# data.frame(Gene = rownames(tmp$p.value), t = round(tmp$statistic, 2),
#               "p-value" = signif(tmp$p.value, 3))[order(tmp$p.value)[1:10], ]
# 
# 
# mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE)
#  # The reults can be summarized by the generic function summary, which takes the msviper object and
#  # either the number of top regulators to report or a specific set of regulators to list. The default for this
#  # parameter is the top 10 master regulators (MRs).
# summary(mrs)
