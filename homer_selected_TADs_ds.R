options(scipen=100)

setDir=""

# Rscript homer_selected_TADs_ds.R <hicds> <exprds> <TAD>
# Rscript homer_selected_TADs_ds.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF chr1_TAD12
# Rscript homer_selected_TADs_ds.R LG1_40kb TCGAluad_norm_luad chr12_TAD32

# Rscript homer_selected_TADs_ds.R LG1_40kb TCGAluad_norm_luad chr12_TAD32
# Rscript homer_selected_TADs_ds.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr13_TAD244
# Rscript homer_selected_TADs_ds.R LG2_40kb TCGAlusc_norm_lusc chr14_TAD144
# Rscript homer_selected_TADs_ds.R ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad chr19_TAD192
# Rscript homer_selected_TADs_ds.R GSE118514_RWPE1_40kb TCGAprad_norm_prad chr17_TAD174
# Rscript homer_selected_TADs_ds.R GSE118514_RWPE1_40kb TCGAprad_norm_prad chr1_TAD776
# Rscript homer_selected_TADs_ds.R GSE118514_22Rv1_40kb TCGAprad_norm_prad chr17_TAD269
# Rscript homer_selected_TADs_ds.R GSE118514_22Rv1_40kb TCGAprad_norm_prad chr12_TAD196
# Rscript homer_selected_TADs_ds.R GSE118514_RWPE1_40kb TCGAprad_norm_prad chr12_TAD194
# Rscript homer_selected_TADs_ds.R PA2_40kb TCGApaad_wt_mutKRAS chr21_TAD121
# Rscript homer_selected_TADs_ds.R LG2_40kb TCGAlusc_norm_lusc chr1_TAD520
# Rscript homer_selected_TADs_ds.R GSE105381_HepG2_40kb TCGAlihc_norm_lihc chr16_TAD164
# Rscript homer_selected_TADs_ds.R LI_40kb TCGAlihc_norm_lihc chr4_TAD723
# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc chr18_TAD245
# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAluad_norm_luad chr8_TAD83
# Rscript homer_selected_TADs_ds.R LG2_40kb TCGAluad_norm_luad chr16_TAD122
# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc chr11_TAD398
# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc chr1_TAD721
# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc chr3_TAD717
# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc chr14_TAD157
# Rscript homer_selected_TADs_ds.R LG1_40kb TCGAlusc_norm_lusc chr11_TAD390
# Rscript homer_selected_TADs_ds.R LG1_40kb TCGAlusc_norm_lusc chr3_TAD718

# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAluad_norm_luad chr5_TAD170
# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAluad_norm_luad chr11_TAD398
# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAluad_norm_luad chr1_TAD721
# Rscript homer_selected_TADs_ds.R ENCSR444WCZ_A549_40kb TCGAluad_norm_luad chr15_TAD40

# Rscript homer_selected_TADs_ds.R LG1_40kb TCGAluad_norm_luad chr15_TAD44
# Rscript homer_selected_TADs_ds.R LG1_40kb TCGAluad_norm_luad chr12_TAD174

# Rscript homer_selected_TADs_ds.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich chr16_TAD153

# Rscript homer_selected_TADs_ds.R ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad chr12_TAD220

# Rscript homer_selected_TADs_ds.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf chr1_TAD84
# Rscript homer_selected_TADs_ds.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf chr7_TAD568




script_name <- "homer_selected_TADs_ds.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE


require(foreach)
require(doMC)
registerDoMC(ifelse(SSHFS, 2, 40))
source("subtype_cols.R")

nToPlot <- 10

outFolder <- "HOMER_SELECTED_TADS_DS"
dir.create(outFolder, recursive=TRUE)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 3)

hicds <- args[1]
exprds <- args[2]
curr_TAD <- args[3]

final_DT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

final_DT <- final_DT[final_DT$hicds == hicds & final_DT$exprds == exprds,]
stopifnot(nrow(final_DT) > nToPlot)

final_DT <- final_DT[order(final_DT$adjPvalComb),]

stopifnot(curr_TAD %in% final_DT$region)

final_DT <- final_DT[final_DT$region == curr_TAD,]
stopifnot(nrow(final_DT) == 1)

extendBp <- 1000

### BUILD SIGNIF ALONG FDR THRESH
cat("... start retrieving FDR signif. TADs\n")

plotList <- list()

# bin/findMotifsGenome.pl Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/Barutcu_MCF-10A_40kb_TCGAbrca_lum_bas_adjPvalComb_0.01_plus.txt hg19r Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas/adjPvalComb_0.01/MotifOutput_plus -size 200 -p 80
binExec <- file.path("HOMER", "bin","findMotifsGenome.pl" )
homerSize <- 200
homerP <- 80

# for(i_tad in 1:nToPlot) {
i_tad = 1
  
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
  
#} # end-for iterating over TADs to plot

##############################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

