
require(foreach)
require(doMC)
registerDoMC(40)
require(ggpubr)
  # Rscript cmp_pval_permut_allDS.R

plotType <- "svg"
myHeight <- myWidth <- 7
myHeightGG <- 7
myWidthGG <- 9

outFolder <- "CMP_PVAL_PERMUT_ALLDS"
dir.create(outFolder, recursive = TRUE)
  
all_hicds <- list.files("PIPELINE/OUTPUT_FOLDER")
all_hicds1 <- all_hicds[!(grepl("PERMUT", all_hicds) | grepl("RANDOM", all_hicds))]
all_hicds2 <- all_hicds[grepl("RANDOMMIDPOS", all_hicds)]
all_hicds <- c(all_hicds1, all_hicds2)
all_exprds <- sapply(all_hicds, function(x) list.files(file.path("PIPELINE/OUTPUT_FOLDER", x)))


# 
rd_patterns <- c("RANDOMMIDPOS", "RANDOMNBRGENES", "RANDOMSHIFT", "PERMUTG2T", "RANDOMMIDPOSDISC" , "RANDOMMIDPOSSTRICT")



source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

all_ds_adjPval_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
  exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {

	pval_file <- file.path("PIPELINE", "OUTPUT_FOLDER", hicds, exprds, "11sameNbr_runEmpPvalCombined", "emp_pval_combined.Rdata")
	stopifnot(file.exists(pval_file))
	pvals <- get(load(pval_file))
	adj_pvals <- p.adjust(pvals, method="BH")

     data.frame(
   hicds = hicds,
   exprds = exprds,
   adjPval = adj_pvals,
   stringsAsFactors = FALSE)
  }
  exprds_dt
}

all_ds_adjPval_dt$hicds_lab <- gsub(".+RANDOMSHIFT_40kb", "RANDOMSHIFT", 
                               gsub(".+RANDOMNBRGENES_40kb", "RANDOMNBRGENES",
                                    gsub(".+PERMUTG2T_40kb", "PERMUTG2T", 
                                         gsub(".+RANDOMMIDPOSDISC_40kb", "RANDOMMIDPOSDISC",
                                         gsub(".+RANDOMMIDPOSSTRICT_40kb", "RANDOMMIDPOSSTRICT",
                               gsub(".+RANDOMMIDPOS_40kb", "RANDOMMIDPOS", all_ds_adjPval_dt$hicds))))))
all_ds_adjPval_dt$hicds_lab[! all_ds_adjPval_dt$hicds_lab %in% rd_patterns] <- "OBSERVED"



outFile <- file.path(outFolder, paste0("allDS_adjPval_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_ds_adjPval_dt$adjPval, all_ds_adjPval_dt$hicds_lab),
  plotTit = paste0("all DS - n=", length(unique(file.path(all_ds_adjPval_dt$hicds, all_ds_adjPval_dt$exprds))), " - TAD adjPval"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

all_ds_adjPval_dt$adjPval_log10 <- -log10(all_ds_adjPval_dt$adjPval)

outFile <- file.path(outFolder, paste0("allDS_adjPval_density_log10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.4))
par(bty="L")
plot_multiDens(
  split(all_ds_adjPval_dt$adjPval_log10, all_ds_adjPval_dt$hicds_lab),
  plotTit = paste0("all DS - n=", length(unique(file.path(all_ds_adjPval_dt$hicds, all_ds_adjPval_dt$exprds))), " - TAD adjPval [-log10]"), legPos = "topright")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


