options(scipen=100)

setDir=""

# Rscript homer_topTADs_ds.R <hicds> <exprds>
# Rscript homer_topTADs_ds.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF
# Rscript homer_topTADs_ds.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript homer_topTADs_ds.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich

script_name <- "homer_topTADs_ds.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE


require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
source("subtype_cols.R")

nToPlot <- 10

outFolder <- "HOMER_TOPTADS_DS"
dir.create(outFolder, recursive=TRUE)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

hicds <- args[1]
exprds <- args[2]

final_DT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

final_DT <- final_DT[final_DT$hicds == hicds & final_DT$exprds == exprds,]
stopifnot(nrow(final_DT) > nToPlot)

final_DT <- final_DT[order(final_DT$adjPvalComb),]

extendBp <- 1000

### BUILD SIGNIF ALONG FDR THRESH
cat("... start retrieving FDR signif. TADs\n")

plotList <- list()

# bin/findMotifsGenome.pl Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/Barutcu_MCF-10A_40kb_TCGAbrca_lum_bas_adjPvalComb_0.01_plus.txt hg19r Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/adjPvalComb_0.01/MotifOutput_plus -size 200 -p 80
binExec <- file.path("HOMER", "bin","findMotifsGenome.pl" )
homerSize <- 200
homerP <- 80
nToPlot=10
for(i_tad in 1:nToPlot) {
  
  cat("... start top # ", i_tad, "\n")
  
  curr_start <- final_DT$start[i_tad]
  curr_end <- final_DT$end[i_tad]
  curr_exprds <- final_DT$exprds[i_tad]
  curr_hicds <- final_DT$hicds[i_tad]
  curr_TAD <- final_DT$region[i_tad]
  curr_chromo <- gsub("(chr.+)_.+", "\\1", curr_TAD)
  
  # chr6_TAD103	chr6	26200001	26360000	-
  
  outdt <- data.frame(region = c(curr_TAD, curr_TAD),
                        chr = c(curr_chromo, curr_chromo),
                        start = c(curr_start-extendBp, curr_start-extendBp),
                        end = c(curr_end+extendBp, curr_end+extendBp),
                        strand=c("+", "-"),
                        stringsAsFactors = FALSE)
  
  outFile <- file.path(outFolder, 
                       paste0(i_tad, "_", curr_hicds, "_", curr_exprds, "_", curr_TAD, "_plusminus.txt"))
  write.table(outdt, file=outFile, sep="\t", append=FALSE, quote=FALSE, col.names=FALSE, row.names=FALSE)
  cat(paste0("... written: ", outFile, "\n"))
  
  outHomer <- file.path(outFolder, 
                                   paste0(i_tad, "_", curr_hicds, "_", curr_exprds, "_", curr_TAD, "_MotifOutplut_plusminus"))
  
  cmd <- paste(binExec, outFile, "hg19r", outHomer, "-size", homerSize, "-p", homerP)
  
  cat(paste0("> ", cmd, "\n"))
  system(cmd)
  
} # end-for iterating over TADs to plot

##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

