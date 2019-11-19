require(foreach)
require(doMC)

# Rscript top_result_table.R

registerDoMC(40)

setDir <- ""

outFolder <- "TOP_RESULT_TABLE"
dir.create(outFolder, recursive = TRUE)

final_DT <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))

source("subtype_cols.R")

final_DT$dataset <- file.path(final_DT$hicds, final_DT$exprds)
final_DT <- final_DT[order(final_DT$dataset),]
all_ds <- unique(final_DT$dataset)

final_DT <- final_DT[order(final_DT$adjPvalComb),]

ds=all_ds[1]

nTop <- 10

out_dt <- foreach(ds = all_ds, .combine='rbind') %dopar% {
  hicds <- dirname(ds)
  exprds <- basename(ds)
  settingF <- file.path("PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
  stopifnot(file.exists(settingF))
  source(settingF)
  samp1 <- get(load(sample1_file))
  samp2 <- get(load(sample2_file))
  nSamp1 <- length(samp1)
  nSamp2 <- length(samp2)
  curr_dt <- final_DT[final_DT$dataset == ds,][1:nTop,]
  
  data.frame(
    hicds = hicds,
    exprds = exprds,
    ds_type = all_cmps[exprds],
    cond1 = cond1,
    cond2=cond2,
    nSamp1 = nSamp1,
    nSamp2 = nSamp2,
    topTADs = paste0(curr_dt$region, collapse="/"),
    topGenes = paste0(curr_dt$region_genes, collapse="/"),
    topAdjPvalComb = paste0(round(curr_dt$adjPvalComb, 4), collapse="/"),
    topSignifFDR02 = paste0(as.character(curr_dt$signifFDR_0.2), collapse="/"),
    stringsAsFactors = FALSE    
  )
}

outFile <- file.path(outFolder, paste0("topTADs_top_", nTop, "_resultDT.txt"))
write.table(out_dt, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
cat(paste0("... written: ", outFile, "\n"))