require(foreach)
require(doMC)

# Rscript top_result_table_rankDiff.R

registerDoMC(40)

setDir <- ""

outFolder <- "TOP_RESULT_TABLE"
dir.create(outFolder, recursive = TRUE)


source("subtype_cols.R")


options(scipen=100)

SSHFS=F

buildData <- TRUE
# FOCUS OSN TAD !
# Rscript top_result_table_rankDiff.R


all_hicds_norm <- c(
  "LI_40kb",
  "LG1_40kb",
  "LG2_40kb",
  "LG1_40kb",
  "LG2_40kb",
  "LG1_40kb",
  "LG2_40kb",
  "LG1_40kb",
  "LG2_40kb",
  "GSE118514_RWPE1_40kb",
  "GSE118514_RWPE1_40kb"
)
all_hicds_tumor <- c(
  "GSE105381_HepG2_40kb",
  "ENCSR444WCZ_A549_40kb",
  "ENCSR444WCZ_A549_40kb",
  "ENCSR489OCU_NCI-H460_40kb",
  "ENCSR489OCU_NCI-H460_40kb",
  "ENCSR444WCZ_A549_40kb",
  "ENCSR444WCZ_A549_40kb",
  "ENCSR489OCU_NCI-H460_40kb",
  "ENCSR489OCU_NCI-H460_40kb",
  "ENCSR346DCU_LNCaP_40kb",
  "GSE118514_22Rv1_40kb"
  
)
all_exprds <- c(
  "TCGAlihc_norm_lihc",
  "TCGAluad_norm_luad",
  "TCGAluad_norm_luad",
  "TCGAluad_norm_luad",
  "TCGAluad_norm_luad",
  "TCGAlusc_norm_lusc",
  "TCGAlusc_norm_lusc",
  "TCGAlusc_norm_lusc",
  "TCGAlusc_norm_lusc",
  "TCGAprad_norm_prad",
  "TCGAprad_norm_prad"
)
stopifnot(length(all_hicds_norm) == length(all_hicds_tumor))
stopifnot(length(all_hicds_norm) == length(all_exprds))

script_name <- "top_result_table_rankDiff.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(foreach)
require(doMC)
registerDoMC(40)

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")

source("plot_lolliTAD_funct.R")
source("my_heatmap.2.R")
source("annotated_TADs.R"); stopifnot(exists("annotated_tads"))

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.2

outFolder <- file.path("TOP_RESULT_TABLE_RANKDIFF")
dir.create(outFolder, recursive = TRUE)

pipFolder <- file.path("PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))


final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))



all_refs <- c("norm", "tumor")
require(dplyr)



inFolder <- "RANKDIFF_ACTIVDIFF"
stopifnot(dir.exists(inFolder))

nTop <- 10

all_out_dt <- foreach(i = 1:length(all_hicds_norm), .combine='rbind') %dopar% {
  
  hicds_norm <- all_hicds_norm[i]
  hicds_tumor <- all_hicds_tumor[i]
  exprds <- all_exprds[i]
  
  out_dt <- foreach(curr_ref = all_refs, .combine='rbind') %do% {
    
    inFile <- file.path(inFolder, paste0( hicds_norm, "_", hicds_tumor, "_", exprds, "_ref_", curr_ref, ".Rdata"))
    cat(paste0(inFile))
    if(!file.exists(inFile))cat(paste0(inFile))
    
    stopifnot(file.exists(inFile))
    signif_dt <- get(load(inFile))
    top_signif_dt <- signif_dt[1:min(nTop, nrow(signif_dt)),c("ref_hicds", "ref_exprds", "refID", "matching_hicds", "rankDiff")]
    
    final_top_signif_dt <- left_join(top_signif_dt, final_dt, by=c("ref_hicds"="hicds", "ref_exprds" = "exprds", "refID"="region"))
    final_top_signif_dt <- final_top_signif_dt[order(abs(final_top_signif_dt$rankDiff), decreasing = TRUE),]
    
    curr_hicds <- unique(final_top_signif_dt$ref_hicds)
    stopifnot(length(curr_hicds) == 1)
    
    curr_exprds <- unique(final_top_signif_dt$ref_exprds)
    stopifnot(length(curr_exprds) == 1)
    
    curr_matchHicds <- unique(final_top_signif_dt$matching_hicds)
    stopifnot(length(curr_matchHicds) == 1)
    
    stopifnot(!is.na(final_top_signif_dt))
    
    
    data.frame(
      hicds = curr_hicds,
      match_hicds = curr_matchHicds,
      exprds = curr_exprds,
      ds_type = all_cmps[curr_exprds],
      topTADs = paste0(final_top_signif_dt$refID, collapse="/"),
      topGenes = paste0(final_top_signif_dt$region_genes, collapse="/"),
      topAdjPvalComb = paste0(round(final_top_signif_dt$adjPvalComb, 4), collapse="/"),
      topSignifFDR02 = paste0(as.character(final_top_signif_dt$signifFDR_0.2), collapse="/"),
      stringsAsFactors = FALSE    
    )
  }
  
  out_dt
  
}


outFile <- file.path(outFolder, paste0("topTADs_rankDiff_top_", nTop, "_resultDT.txt"))
write.table(all_out_dt, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
cat(paste0("... written: ", outFile, "\n"))