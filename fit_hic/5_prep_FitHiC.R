
require(data.table)

options(scipen=100)

script_name <- "5_prep_FitHiC.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

# INPUT_AGG/KARPAS_DMSO/KARPAS_DMSO_chr1_noDS_merged_agg.txt 
# 0       1       1
# 0       14293   1
# 1       1       156
# 1       2       5


########################################################################################################
########################################################################################################
########################################################################################################
# FRAGSFILE This file stores the information about midpoints (or start indices) of the fragments. 
# 5 columns: 
#   1) chromosome name; 
#   2) will be ignored -> set to 0
#   3) midPoint
#   4) hitCount -> set to 1
#   5) will be ignored -> set to 0
# hitCount only used for filtering; set bins that need to be filtered to “0”  [filtering can be combined with mappability threshold]
# FRAGSFILE Chromosome.Name 	Column.2 	Mid.Point 	Hit.Count 	Column.5
# chr1 	0 	20000 	1 	0
# chr1 	0 	60000 	1 	0


# INTERSFILE This file stores the information about interactions between fragment pairs. 
# 5 columns: 
  # 1) and 3) chromosome name
  # 2) and 4) midPoints of the fragment pair
  # 5) contact count between the two bins. 
# ! midpoints in this file match those in fragsfile
# ! use RAW contact counts and NOT the normalized counts 
# INTERSFILE Chromosome1.Name 	Mid.Point.1 	Chromosome2.Name 	Mid.Point.2 	Hit.Count
# chr1 	100020000 	chr1 	100100000 	201
# chr1 	100020000 	chr1 	100140000 	232

# Rscript 5_prep_FitHiC.R INPUT_AGG/KARPAS_DMSO KARPAS_DMSO noDS_merged_agg.txt chr1 10000 PREP_FITHIC


input_dir <- file.path("INPUT_AGG", "KARPAS_DMSO")
file_prefix <- "KARPAS_DMSO"
file_suffix <- "noDS_merged_agg.txt"
chromo <- "chr1"
binResol <- 10*10^3
ouptut_dir <- "PREP_FITHIC"

args <- commandArgs(trailingOnly = TRUE)
cat(args, "\n")

stopifnot(length(args) == 6)
input_dir <- args[1]
file_prefix <- args[2]
file_suffix <- args[3]
chromo <- args[4]
binResol <- args[5]
output_dir <- args[6]

binResol <- as.numeric(binResol)
stopifnot(!is.na(binResol))

cat("> start with: \n")
cat(paste0(".... chromo\t:\t", chromo, "\n"))
cat(paste0(".... binResol\t:\t", binResol, "\n"))

#################################

input_file <- file.path(input_dir, paste0(file_prefix, "_", chromo, "_",  file_suffix))
cat(input_file, "\n")
stopifnot(file.exists(input_file))

dir.create(output_dir, recursive = TRUE)

aggDT <- read.table(input_file, header=FALSE, col.names=c("binA", "binB", "count"))
stopifnot(!is.na(aggDT))

aggDT$binA_bp <- aggDT$binA * binResol
aggDT$binB_bp <- aggDT$binB * binResol

aggDT$binA_midpoint <- aggDT$binA_bp + binResol/2
aggDT$binB_midpoint <- aggDT$binB_bp + binResol/2

max_midpoint <- max(aggDT$binB_midpoint)
stopifnot(max_midpoint == aggDT$binB_midpoint[nrow(aggDT)])


frags_DT <- data.frame(
  chr = chromo,
  col2 = 0,
  midpoint = seq(from=binResol/2, to=max_midpoint, by=binResol),
  hitcount = 1,
  col5 = 0,
  stringsAsFactors = FALSE
)

inters_DT <- data.frame(
  chrA = chromo,
  midpointA = aggDT$binA_midpoint,
  chrB = chromo,
  midpointB = aggDT$binB_midpoint,
  hitcount = aggDT$count,
  stringsAsFactors = FALSE
)

stopifnot(inters_DT$midpointA %in% frags_DT$midpoint)
stopifnot(inters_DT$midpointB %in% frags_DT$midpoint)


output_file_frags <- file.path(output_dir, paste0(file_prefix, "_", chromo, "_", gsub(".txt", "", file_suffix), "_FitHiC_fragsfile.txt"))
write.table(frags_DT, file = output_file_frags, col.names = FALSE, row.names=FALSE, quote=FALSE, append=FALSE, sep="\t")
cat(paste0("... written: ", output_file_frags, "\n"))


output_file_inters <- file.path(output_dir, paste0(file_prefix, "_", chromo, "_", gsub(".txt", "", file_suffix), "_FitHiC_intersfile.txt"))
write.table(inters_DT, file = output_file_inters, col.names = FALSE, row.names=FALSE, quote=FALSE, append=FALSE, sep="\t")
cat(paste0("... written: ", output_file_inters, "\n"))


########################################################################################################
########################################################################################################
########################################################################################################
cat("*** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))