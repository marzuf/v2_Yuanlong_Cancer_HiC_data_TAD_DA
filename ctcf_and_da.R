
# Rscript ctcf_and_da.R

library("readxl")
library(doMC)
library(foreach)

registerDoMC(40)
# runFolder <- "../v2_Yuanlong_Cancer_HiC_data_TAD_DA_GM12878" #PIPELINE/OUTPUT_FOLDER/GM12878_40kb/TCGAluad_norm_luad/11sameNbr_runEmpPvalCombined/"
# hicds <- "GM12878_40kb"

runFolder <- "." 
hicds <- "ENCSR444WCZ_A549_40kb"
exprds <- "TCGAluad_mutKRAS_mutEGFR"

final_dt <- get(load("CREATE_FINAL_TABLE/all_result_dt.Rdata"))
ds_final_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds, ]
stopifnot(nrow(ds_final_dt) > 0)

buildTable <- TRUE

outFolder <- file.path("CTCF_AND_DA", hicds, exprds)
dir.create(outFolder, recursive = TRUE)

tad_dt <- read.delim(file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt"), stringsAsFactors = FALSE, 
                     header=F, col.names = c("chromo", "region", "start", "end"))

### KEEP ONLY TAD REGIONS
tad_dt <- tad_dt[grepl("_TAD", tad_dt$region),]
stopifnot(!grepl("BOUND", tad_dt$region))

ctcf_dt <- read_excel("13059_2020_2108_MOESM2_ESM.xlsx", sheet="CTCFs")
ctcf_dt <- as.data.frame(ctcf_dt)
ctcf_dt <- ctcf_dt[, 1:7]

# assign ctcf BS to tads
ctcf_dt$chr <- as.character(ctcf_dt$chr)
ctcf_dt <- ctcf_dt[ctcf_dt$chr %in% tad_dt$chromo,]
stopifnot(nrow(ctcf_dt) > 0)
stopifnot(is.numeric(ctcf_dt$start))
stopifnot(is.numeric(ctcf_dt$end))
stopifnot(ctcf_dt$start <= ctcf_dt$end)
ctcf_dt$region <- NA

i=1
if(buildTable){
  ctcf2tad_dt <- foreach(i = 1:nrow(ctcf_dt), .combine='rbind') %dopar% {
    
    chr <- ctcf_dt$chr[i]
    ctcf_start <- ctcf_dt$start[i]
    ctcf_end <- ctcf_dt$end[i]
    
    subtad_dt <- tad_dt[tad_dt$chromo == chr,]
    stopifnot(nrow(subtad_dt) > 0)
    
    # assign if start after tad start and end before tad end
    test1 <- which(ctcf_start >= subtad_dt$start & ctcf_end <= subtad_dt$end)
    test2 <- which(subtad_dt$start <= ctcf_start & subtad_dt$end >= ctcf_end)
    stopifnot(test1==test2)
    stopifnot(length(test1) == 0 | length(test1) == 1)
    if(length(test1) == 1) {
      ctcf_dt$region[i] <- subtad_dt$region[test1]
    } else {
      ctcf_dt$region[i] <- NA
    }
    ctcf_dt[i,]  
  }
  outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
  save(ctcf2tad_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile,"\n"))
} else {
  outFile <- file.path(outFolder, "ctcf2tad_dt.Rdata")
  ctcf2tad_dt <- get(load(outFile))
}
# load("CTCF_AND_DA/ENCSR444WCZ_A549_40kb/TCGAluad_mutKRAS_mutEGFR/ctcf2tad_dt.Rdata")
cat(paste0("# of CTCF BS:\t",nrow(ctcf2tad_dt), "\n"))
cat(paste0("# of CTCF BS in TADs:\t",sum(!is.na(ctcf2tad_dt$region)), "\n"))
cat(paste0("# of CTCF BS out of TADs:\t",sum(is.na(ctcf2tad_dt$region)), "\n"))

merged_dt <- merge(ds_final_dt[,c("hicds", "exprds", "region", "meanLogFC", "meanCorr", "adjPvalComb")],
                   ctcf2tad_dt,
                   by="region", all=FALSE)
tmp <- merged_dt$region
tmp <- gsub("(.+)_.+", "\\1", tmp)
stopifnot(tmp == merged_dt$chr)


aggByOrientation_dt <- aggregate(chr ~ orientation + region + meanCorr + adjPvalComb, FUN=length, data=merged_dt)
colnames(aggByOrientation_dt)[colnames(aggByOrientation_dt) == "chr"] <- "CTCF_count"

aggByOrientation_dt$orientation_lab <- ifelse(aggByOrientation_dt$orientation == ">", "forward", 
                                          ifelse(aggByOrientation_dt$orientation == "<", "reverse", NA))
stopifnot(!is.na(aggByOrientation_dt$orientation_lab))

wide_aggByOrientation_dt <- reshape(aggByOrientation_dt[,c("region", "orientation_lab","CTCF_count")], 
                                    idvar="region", direction="wide", timevar = "orientation_lab")



agg_dt <- aggregate(chr ~ region + meanCorr + adjPvalComb, FUN=length, data=merged_dt)
colnames(agg_dt)[colnames(agg_dt) == "chr"] <- "CTCF_totCount"

agg_merged_dt <- merge(wide_aggByOrientation_dt, agg_dt, by="region", all=TRUE)
stopifnot(!is.na(agg_merged_dt))

agg_merged_dt$CTCF_count.forward[is.na(agg_merged_dt$CTCF_count.forward)] <- 0
agg_merged_dt$CTCF_count.reverse[is.na(agg_merged_dt$CTCF_count.reverse)] <- 0

stopifnot(agg_merged_dt$CTCF_count.reverse + agg_merged_dt$CTCF_count.forward == agg_merged_dt$CTCF_totCount)

plot(
  agg_merged_dt$CTCF_totCount ~ -log10(agg_merged_dt$adjPvalComb)
)
