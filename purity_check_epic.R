########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript purity_check_epic.R\n"))

script_name <- "purity_check_epic.R"

### !!! NEED TO HAVE RUN DELTA_SAMPLE_PURITY BEFORE !!!

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# Rscript purity_check_epic.R <gene_symbol>
# Rscript purity_check_epic.R GIMAP4
all_gene_symbols=c("GIMAP2", "GIMAP1")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0)
all_gene_symbols <- args[1:length(args)]

# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotType <- "png"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4
myHeightGG <- 7
myWidthGG <- 7

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

mainFolder <- file.path(".")
pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path("PIPELINE", "INPUT_FILES")

outFolder <- file.path("PURITY_CHECK_EPIC")
dir.create(outFolder, recursive = TRUE)

final_dt_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
final_dt <- get(load(final_dt_file))

purity_file <- file.path("DELTA_SAMPLE_PURITY_EPIC", "purity_delta_dt.Rdata")
purity_dt <- get(load(purity_file))
purity_dt <- data.frame(purity_dt)
all_pm_metrics <- colnames(purity_dt) #"Bcells"      "CAFs"        "CD4_Tcells"  "CD8_Tcells"  "Endothelial" "Macrophages" "NKcells"     "otherCells" 
pm <- "otherCells"


matching_symbols <- sapply(final_dt$region_genes, function(x) any(all_gene_symbols %in% unlist(strsplit(x=x, split=","))))
matching_dt <- final_dt[matching_symbols,]
  
if(nrow(matching_dt) == 0) stop(paste0("... could not find any TAD for:\t", paste0(all_gene_symbols, collapse=","), "\n"))

purity_dt$hicds <- dirname(rownames(purity_dt))
purity_dt$exprds <- basename(rownames(purity_dt))

stopifnot(matching_dt$hicds %in% purity_dt$hicds)
stopifnot(matching_dt$exprds %in% purity_dt$exprds)

matching_purity_dt <- merge(matching_dt, purity_dt, by=c("hicds", "exprds"), all.x=TRUE, all.y=FALSE)
matching_purity_dt$labels <- paste0(matching_purity_dt$hicds, "\n", matching_purity_dt$exprds)

matching_purity_dt$adjPvalComb_log10 <- -log10(matching_purity_dt$adjPvalComb)

print_dt <- matching_purity_dt[,c("hicds", "exprds", "region", "region_genes", "meanLogFC", "adjPvalComb", all_pm_metrics)]

plotSubTit <- ""

tad_val="meanLogFC"
pm="ABSOLUTE"

for(tad_val in c("meanLogFC", "adjPvalComb_log10")) {
  stopifnot(tad_val %in% colnames(matching_purity_dt))
  for(pm in all_pm_metrics) {
    stopifnot(pm %in% colnames(matching_purity_dt))
    myylab <- paste0("delta purity (", pm, ")")
    myxlab <- paste0(tad_val)
    
    outFile <- file.path(outFolder, paste0(paste0(all_gene_symbols, collapse="_"), "_", pm, "_vs_", tad_val, ".", plotType))
    do.call(plotType, list(file=outFile, height=myHeight, width=myWidth))
    par(bty="l")
    plot(
      as.formula(paste0(pm, "~", tad_val)), data = matching_purity_dt,
      xlab = myxlab,
      ylab = myylab,
      cex=0.7,
      pch=16,
      cex.lab = plotCex,
      cex.axis = plotCex
      )
    text(x=matching_purity_dt[,paste0(tad_val)],
         y=matching_purity_dt[,paste0(pm)],
         labels = matching_purity_dt[,paste0("labels")],
         cex=0.7
         )
    mtext(side=3, text = plotSubTit, cex=1.4)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))    

  } # end-for iterating over purity metrics
} # end-for iterating over TAD values

write.table(print_dt[,c("hicds", "exprds", "region", "region_genes")], sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

outFile <- file.path(outFolder, paste0(paste0(all_gene_symbols, collapse="_"), "_matching_purity_dt.txt"))
write.table(print_dt, file=outFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
cat(paste0("... written: ", outFile, "\n"))    


##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

