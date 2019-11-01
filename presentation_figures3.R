

options(scipen=100)

# Rscript presentation_figures3.R

script_name <- "presentation_figures3.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
require(reshape2)

registerDoMC(4)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")

all_cols[all_cols == "red"] <- "brown3"
all_cols[all_cols == "blue"] <- "darkblue"
all_cols[all_cols == "green"] <- "forestgreen"

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

setDir <- "/media/electron"
setDir <- ""

mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- "PRESENTATION_FIGURES3"
dir.create(outFolder, recursive = TRUE)


ex_hicds <- "GSE99051_786_O_40kb"
ex_exprds <- "TCGAkich_norm_kich"
hicds_tit <- "786-O"
exprds_tit <- "norm_vs_kich"


plotCex <- 1.2

signifThresh <- 0.01


final_DT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
final_DT_signif <- final_DT[final_DT$adjPvalComb <= signifThresh,]


ex_DT <- final_DT[final_DT$hicds == ex_hicds &
                    final_DT$exprds == ex_exprds,
                  ]


ex_DT <- ex_DT[order(ex_DT$adjPvalComb),]

ex_DT[order(ex_DT$adjPvalComb),][1:10,]

ex_DT_signif <- ex_DT[ex_DT$adjPvalComb <= signifThresh,]

mean(ex_DT$meanLogFC)
mean(ex_DT$meanCorr)
mean(ex_DT_signif$meanLogFC)
mean(ex_DT_signif$meanCorr)

mean(final_DT$meanLogFC)
mean(final_DT$meanCorr)
mean(final_DT_signif$meanLogFC)
mean(final_DT_signif$meanCorr)

mean(abs(ex_DT$meanLogFC))
mean(abs(ex_DT_signif$meanLogFC))
mean(abs(final_DT$meanLogFC))
mean(abs(final_DT_signif$meanLogFC))

nrow(ex_DT_signif)/nrow(ex_DT) * 100
nrow(final_DT_signif)/nrow(final_DT) * 100



min(abs(ex_DT$meanLogFC))
min(ex_DT$meanCorr)
min(abs(ex_DT_signif$meanLogFC))
min(ex_DT_signif$meanCorr)

min(abs(final_DT$meanLogFC))
min(final_DT$meanCorr)
min(abs(final_DT_signif$meanLogFC))
min(final_DT_signif$meanCorr)


outFile <- file.path(outFolder, paste0(ex_hicds, "_", ex_exprds, "_meanIntraCorr_meanLogFC_dotplot_with_signif.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

ex_DT$dotcol <- ifelse(ex_DT$adjPvalComb <= signifThresh,"red", "darkgrey")
par(bty="l")
plot(
  main = paste0(hicds_tit, " - ", exprds_tit),
  x = ex_DT$meanLogFC,
  y = ex_DT$meanCorr,
  xlab = "meanLogFC",
  ylab = "meanIntraCorr",
  pch=16,
  col = ex_DT$dotcol,
  cex=0.7,
  cex.axis=plotCex,
  cex.lab=plotCex
)
# text(x = ex_DT$meanLogFC[1:10],
#      y = ex_DT$meanCorr[1:10],
#      labels = ex_DT$region[1:10], pos = 3)
points(x = ex_DT$meanLogFC[1:10],
     y = ex_DT$meanCorr[1:10],
     col="red", cex=1.4, lwd=2)
mtext(side=3, text = paste0("# TADs = ", nrow(ex_DT), "; # signif. TADs = ", nrow(ex_DT_signif)))
# legend("topleft", legend=paste0("adj. p-val <= ", signifThresh), pch=16, col="red", bty="n")
# legend("bottomright", legend=paste0("adj. p-val <= ", signifThresh), pch=16, col="red", bty="n")
# legend("bottomright", legend=paste0("top 10"), pch=1, col="red", bty="n")
legend("bottomright", legend=c(paste0("adj. p-val <= ", signifThresh), "top 10"), pch=c(16, 1), col="red", bty="n")
# legend("bottomright", legend=paste0("top 10"), pch=1, col="red", bty="n")



foo <- dev.off()

stop("--ok\n")


final_DT_signif <- final_DT[final_DT$adjPvalComb <=signifThresh,]

nrow(final_DT)# 100884
nrow(final_DT_signif) # 1453
nrow(final_DT_signif)/58
# [1] 25.05172
nrow(final_DT)/58
# [1] 1739.379


Rscript plot_lolli_by_tad_dataset.R chr11_TAD482 GSE99051_786_O_40kb TCGAkich_norm_kich
Rscript plot_lolli_by_tad_dataset.R chr14_TAD162 GSE99051_786_O_40kb TCGAkich_norm_kich no annot
Rscript plot_lolli_by_tad_dataset.R chr16_TAD160 GSE99051_786_O_40kb TCGAkich_norm_kich MT1 -> GOOD
Rscript plot_lolli_by_tad_dataset.R chr16_TAD4 GSE99051_786_O_40kb TCGAkich_norm_kich no annot
Rscript plot_lolli_by_tad_dataset.R chr18_TAD156 GSE99051_786_O_40kb TCGAkich_norm_kich no annot
Rscript plot_lolli_by_tad_dataset.R chr19_TAD207 GSE99051_786_O_40kb TCGAkich_norm_kich no annot
Rscript plot_lolli_by_tad_dataset.R chr7_TAD116 GSE99051_786_O_40kb TCGAkich_norm_kich  no annot
Rscript plot_lolli_by_tad_dataset.R chr7_TAD438 GSE99051_786_O_40kb TCGAkich_norm_kich  MET, CAV1 -> GOOD
Rscript plot_lolli_by_tad_dataset.R chr10_TAD23 GSE99051_786_O_40kb TCGAkich_norm_kich akr1c
Rscript plot_lolli_by_tad_dataset.R chr6_TAD124 GSE99051_786_O_40kb TCGAkich_norm_kich HLA


# MET, CAv1, cav2; akr1c1, akr1c2, akr1c3
tolookgenes <- c("MET", "CAV1", "CAV2", "AKR1C1", "AKR1C2", "AKR1C3")


setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

tolookentrez <- names(entrez2symb)[entrez2symb %in% tolookgenes]


load("GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")
gene_tad_dt <- all_gene_tad_signif_dt[all_gene_tad_signif_dt$hicds == ex_hicds & all_gene_tad_signif_dt$exprds == ex_exprds,]

sub_gene_tad_dt <- gene_tad_dt[as.character(gene_tad_dt$entrezID) %in% as.character(tolookentrez),]
sub_gene_tad_dt$geneSymb <- entrez2symb[sub_gene_tad_dt$entrezID]

# gene_rank geneSymb
# 327646      5560   AKR1C1
# 327647       706   AKR1C2
# 329760       829      MET
# 333383      3503     CAV1
# 333389       465     CAV2
# 333408        26   AKR1C3


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))




